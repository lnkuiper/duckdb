#include "duckdb/reoptimizer/reoptimizer.hpp"

#include "duckdb/catalog/catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/unordered_map.hpp"
#include "duckdb/common/enums/catalog_type.hpp"
#include "duckdb/common/enums/expression_type.hpp"
#include "duckdb/common/enums/logical_operator_type.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/optimizer/join_order_optimizer.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_join.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"
#include "duckdb/storage/data_table.hpp"

using namespace duckdb;
using namespace std;

// TODO: write tests that figure out if this shit actually works
ReOptimizer::ReOptimizer(ClientContext &context, Binder &binder) : context(context), binder(binder) {
}

unique_ptr<LogicalOperator> ReOptimizer::ReOptimize(unique_ptr<LogicalOperator> plan, const string query) {
	const string tablename_prefix = "_reopt_temp_" + to_string(hash<string>{}(query));
	// re-optimization loop
	for (int iter = 0; true; iter++) {
		// profiling for each iteration
		context.profiler.StartPhase(to_string(iter));
		// create and execute subquery, adjust plan accordingly
		const string temp_table_name = tablename_prefix + "_" + to_string(iter);
		plan = PerformPartialPlan(move(plan), temp_table_name);
		if (done) {
			context.profiler.EndPhase();
			break;
		}
		// TODO: if q-error is low we could get away without calling optimizer (or some other criterion)
		context.profiler.StartPhase("optimizer");
		plan = CallOptimizer(move(plan));
		context.profiler.EndPhase();
		// if (ExtractJoinOperators(*plan).size() <= 1)
		// 	break;
		// end iteration profiling phase
		context.profiler.EndPhase();
	}
	context.profiler.StartPhase("clear_left_proj_maps");
	plan = ClearLeftProjectionMaps(move(plan));
	context.profiler.EndPhase();
	return plan;
}

unique_ptr<LogicalOperator> ReOptimizer::PerformPartialPlan(unique_ptr<LogicalOperator> plan,
                                                            const string temporary_table_name) {
	context.profiler.StartPhase("decide_subquery");
	auto *subquery_plan = DecideSubQueryPlan(*plan);
	context.profiler.EndPhase();

	if (done)
		return plan;
	auto *join_subquery_plan = static_cast<LogicalComparisonJoin *>(subquery_plan);

	context.profiler.StartPhase("generate_projection_maps");
	plan = GenerateProjectionMaps(move(plan));
	context.profiler.EndPhase();

	plan->Print();

	Printer::Print("----------------------------- before");
	plan->children[0]->GetColumnBindings();
	string bcres0 = "PROJECTION ";
	for (idx_t i = 0; i < plan->children[0]->expressions.size(); i++) {
		auto &e = *plan->children[0]->expressions[i];
		if (e.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
			BoundColumnRefExpression bcre = (BoundColumnRefExpression &) e;
			bcres0 += bcre.ToString();
		}
	}
	Printer::Print(bcres0);
	plan->children[0]->children[0]->GetColumnBindings();
	plan->children[0]->children[0]->children[0]->GetColumnBindings();
	Printer::Print("-----------------------------");

	context.profiler.StartPhase("map_binding_names");
	CreateMaps(*plan);
	context.profiler.EndPhase();

	context.profiler.StartPhase("create_subquery");
	vector<string> queried_tables;
	vector<string> where_conditions;
	string subquery = CreateSubQuery(*join_subquery_plan, temporary_table_name, queried_tables, where_conditions);
	context.profiler.EndPhase();

	Printer::Print(subquery);

	context.profiler.StartPhase("execute_subquery");
	ExecuteSubQuery(subquery);
	context.profiler.EndPhase();

	context.profiler.StartPhase("inject_cardinalities");
	InjectCardinalities(*join_subquery_plan, temporary_table_name);
	context.profiler.EndPhase();

	context.profiler.StartPhase("adjust_plan");
	plan = AdjustPlan(move(plan), *join_subquery_plan, temporary_table_name);
	context.profiler.EndPhase();

	context.profiler.StartPhase("fix_bindings");
	FixColumnBindings(*plan);
	context.profiler.EndPhase();

	// TODO: implement
	// context.profiler.StartPhase("create_next_plan");
	// plan = CreateNextPlan(move(plan));
	// context.profiler.EndPhase();

	plan->Print();

	Printer::Print("----------------------------- after");
	plan->children[0]->GetColumnBindings();
	string bcres1 = "PROJECTION ";
	for (idx_t i = 0; i < plan->children[0]->expressions.size(); i++) {
		auto &e = *plan->children[0]->expressions[i];
		if (e.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
			BoundColumnRefExpression bcre = (BoundColumnRefExpression &) e;
			bcres1 += bcre.ToString();
		}
	}
	Printer::Print(bcres1);
	plan->children[0]->children[0]->GetColumnBindings();
	plan->children[0]->children[0]->children[0]->GetColumnBindings();
	Printer::Print("-----------------------------");

	return plan;
}

LogicalOperator *ReOptimizer::DecideSubQueryPlan(LogicalOperator &plan) {
	vector<LogicalOperator *> joins = ExtractJoinOperators(plan);
	if (joins.size() <= 1) {
		done = true;
		return &plan;
	}

	// TODO: implement selection procedure
	// select last one for now - deciding this is basically the whole re-optimization strategy
	// the selection procedure might also decide that we are done, instead of selecting an index
	idx_t selected_join_i = joins.size() - 1;
	return joins[selected_join_i];
}

vector<LogicalOperator *> ReOptimizer::ExtractJoinOperators(LogicalOperator &plan) {
	// the deepest join operators appear at the end of the vector - DuckDB only uses right-heavy
	vector<LogicalOperator *> joins;
	if (plan.children.empty())
		return joins;
	switch (plan.type) {
	// case LogicalOperatorType::JOIN:
	// case LogicalOperatorType::ANY_JOIN:
	// case LogicalOperatorType::DELIM_JOIN:
	// case LogicalOperatorType::CROSS_PRODUCT:
	case LogicalOperatorType::COMPARISON_JOIN: // Pretty sure we only get COMPARISON_JOIN in JOB
		joins.push_back(&plan);
		break;
	default:
		// nothing to do
		break;
	}
	for (auto &child : plan.children) {
		vector<LogicalOperator *> child_joins = ExtractJoinOperators(*child);
		joins.insert(joins.end(), child_joins.begin(), child_joins.end());
	}
	return joins;
}

//! Gets column bindings from a join, but assumes projection maps are actually filled
static vector<ColumnBinding> GetColumnBindings(LogicalComparisonJoin &join) {
	vector<ColumnBinding> cbs;
	if (!join.left_projection_map.empty()) {
		vector<ColumnBinding> left_cbs =
		    LogicalOperator::MapBindings(join.children[0]->GetColumnBindings(), join.left_projection_map);
		for (auto cb : left_cbs)
			cbs.push_back(cb);
	}
	if (!join.right_projection_map.empty()) {
		vector<ColumnBinding> right_cbs =
		    LogicalOperator::MapBindings(join.children[1]->GetColumnBindings(), join.right_projection_map);
		for (auto cb : right_cbs)
			cbs.push_back(cb);
	}
	return cbs;
}

unique_ptr<LogicalOperator> ReOptimizer::GenerateProjectionMaps(unique_ptr<LogicalOperator> plan) {
	bool done = true;
	for (auto &child : plan->children) {
		if (child->type == LogicalOperatorType::COMPARISON_JOIN || child->type == LogicalOperatorType::PROJECTION) {
			done = false;
			break;
		}
	}
	if (done)
		return plan;
	// store the column bindings that are needed by the parents of the child join
	vector<ColumnBinding> column_bindings;
	if (plan->type == LogicalOperatorType::PROJECTION) {
		// from expressions
		for (idx_t expr_i = 0; expr_i < plan->expressions.size(); expr_i++) {
			auto &expr = *plan->expressions[expr_i];
			if (expr.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
				continue;
			auto bcre = (BoundColumnRefExpression &)expr;
			column_bindings.push_back(bcre.binding);
		}
		// get the child join, modify its left projection map
		LogicalComparisonJoin *child_join;
		if (plan->children[0]->type == LogicalOperatorType::COMPARISON_JOIN)
			child_join = static_cast<LogicalComparisonJoin *>(plan->children[0].get());
		// FIXME: #1 PRIORITY use LogicalAggregate::GetColumnBindings to fix this shit, or some other method i dont know anymore...
		// LogicalOperator::expressions is called select_list for LogicalAggregate
		// All of these have class/type BoundAggregateExpression
		// We can't do anything with this because we get some deleted pointer error
		// How the fuck are we gonna get the bindings then???
		else if (plan->children[0]->type != LogicalOperatorType::COMPARISON_JOIN)
			child_join = static_cast<LogicalComparisonJoin *>(plan->children[0]->children[0].get());
		auto left_cbs = child_join->children[0]->GetColumnBindings();
		// add index of needed bindings
		for (column_t cb_i = 0; cb_i < left_cbs.size(); cb_i++) {
			for (auto cb : column_bindings) {
				if (left_cbs[cb_i] == cb) {
					child_join->left_projection_map.push_back(cb_i);
				}
			}
		}
	} else if (plan->type == LogicalOperatorType::COMPARISON_JOIN) {
		auto *join = static_cast<LogicalComparisonJoin *>(plan.get());
		for (idx_t cond_i = 0; cond_i < join->conditions.size(); cond_i++) {
			// from conditions
			JoinCondition &jc = join->conditions[cond_i];
			auto l = (BoundColumnRefExpression &)*jc.left.get();
			auto r = (BoundColumnRefExpression &)*jc.right.get();
			column_bindings.push_back(l.binding);
			column_bindings.push_back(r.binding);
		}
		// from bindings of this plan
		for (auto new_cb : GetColumnBindings(*join)) {
			bool found = false;
			for (auto existing_cb : column_bindings) {
				if (existing_cb == new_cb) {
					found = true;
					break;
				}
			}
			if (!found)
				column_bindings.push_back(new_cb);
		}
		// FIXME: something small is going wrong here
		Printer::Print(to_string(column_bindings.size()));

		// add index of needed bindings, and keep track of removed columns
		auto *child_join = static_cast<LogicalComparisonJoin *>(plan->children[1].get());
		vector<ColumnBinding> left_cbs =
		    LogicalOperator::MapBindings(child_join->children[0]->GetColumnBindings(), child_join->left_projection_map);
		vector<column_t> removed_columns;
		for (column_t cb_i = 0; cb_i < left_cbs.size(); cb_i++) {
			bool keep = false;
			for (auto cb : column_bindings) {
				if (left_cbs[cb_i] == cb) {
					keep = true;
					break;
				}
			}
			if (keep) {
				child_join->left_projection_map.push_back(cb_i);
			} else {
				removed_columns.push_back(cb_i);
			}
		}
		// update right projection map accordingly
		for (idx_t rpi = 0; rpi < join->right_projection_map.size(); rpi++) {
			column_t decr = 0;
			for (column_t removed_col : removed_columns) {
				if (removed_col < join->right_projection_map[rpi]) {
					decr++;
				}
			}
			join->right_projection_map[rpi] -= decr;
		}
	}
	// recursively propagate operator tree
	for (idx_t child_i = 0; child_i < plan->children.size(); child_i++) {
		plan->children[child_i] = GenerateProjectionMaps(move(plan->children[child_i]));
	}
	return plan;
}

void ReOptimizer::CreateMaps(LogicalOperator &plan) {
	if (plan.type == LogicalOperatorType::GET) {
		auto *get = static_cast<LogicalGet *>(&plan);
		for (idx_t i = 0; i < get->column_ids.size(); i++) {
			bta["#[" + std::to_string(get->table_index) + "." + std::to_string(i) + "]"] = get->table->columns[get->column_ids[i]].name;
		}
	}
	// recursively propagate operator tree
	for (auto &child : plan.children) {
		CreateMaps(*child);
	}
}

// TODO: perhaps instead of sending a subquery, we add project + create table onto the plan and hack it in
string ReOptimizer::CreateSubQuery(LogicalComparisonJoin &join, const string temporary_table_name,
                                   vector<string> &queried_tables, vector<string> &where_conditions) {
	// join conditions
	for (idx_t cond_i = 0; cond_i < join.conditions.size(); cond_i++) {
		JoinCondition &join_condition = join.conditions[cond_i];
		auto l_bind = ((BoundColumnRefExpression &)*join_condition.left.get()).binding;
		auto r_bind = ((BoundColumnRefExpression &)*join_condition.right.get()).binding;
		where_conditions.push_back("t" + to_string(l_bind.table_index) + "." + bta[l_bind.ToString()] + " = " + "t" +
		                           to_string(r_bind.table_index) + "." + bta[r_bind.ToString()]);
	}
	for (idx_t child_i = 0; child_i < join.children.size(); child_i++) {
		auto &child = join.children[child_i];
		if (child->type == LogicalOperatorType::COMPARISON_JOIN) {
			// recursively propagate operator tree
			auto child_join = static_cast<LogicalComparisonJoin *>(child.get());
			CreateSubQuery(*child_join, temporary_table_name, queried_tables, where_conditions);
		} else {
			LogicalGet *logical_get;
			switch (child->type) {
			case LogicalOperatorType::GET: {
				logical_get = static_cast<LogicalGet *>(child.get());
				break;
			}
			case LogicalOperatorType::FILTER: {
				// extract filter conditions
				assert(child->children.size() == 1);
				auto *logical_filter = static_cast<LogicalFilter *>(child.get());
				for (auto &expr : logical_filter->expressions) {
					if (expr->GetExpressionClass() != ExpressionClass::BOUND_COMPARISON)
						throw new ReOptimizerException("Expected filter class to be BOUND_COMPARISON, got '%s'",
						                               expr->GetExpressionClass());
					auto *bce = static_cast<BoundComparisonExpression *>(expr.get());
					// filter conditions -table reference is 'always' placed on the left, constant on the right
					if (bce->GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
						throw new ReOptimizerException("Expected LHS of comparison to be BOUND_COLUMN_REF, got '%s'",
						                               bce->GetExpressionClass());
					auto binding = static_cast<BoundColumnRefExpression *>(bce->left.get())->binding;
					string lhs = "t" + to_string(binding.table_index) + "." + bta[binding.ToString()];
					if (bce->right->return_type == TypeId::VARCHAR) {
						// strings need to be escaped with quotes
						where_conditions.push_back(lhs + ExpressionTypeToOperator(bce->type) + "'" +
						                           bce->right->ToString() + "'");
					} else {
						where_conditions.push_back(lhs + ExpressionTypeToOperator(bce->type) + bce->right->ToString());
					}
				}
				logical_get = static_cast<LogicalGet *>(logical_filter->children[0].get());
				break;
			}
			default:
				throw new ReOptimizerException("Expected child of JOIN to be JOIN, GET or FILTER, got '%s'",
				                               child->type);
			}
			// extract queried tables
			queried_tables.push_back(logical_get->table->schema->name + "." + logical_get->table->name + " AS t" +
			                         to_string(logical_get->table_index));
		}
	}
	// selected columns only for the top-level join (can't use join.GetColumnBindings() because empty proj map -> all)
	vector<string> selected_columns;
	for (auto cb : GetColumnBindings(join))
		selected_columns.push_back("t" + to_string(cb.table_index) + "." + bta[cb.ToString()]);
	// create query and return
	return "CREATE TEMPORARY TABLE main." + temporary_table_name + " AS (" + "SELECT " +
	       JoinStrings(selected_columns, ", ") + " " + "FROM " + JoinStrings(queried_tables, ", ") + " " + "WHERE " +
	       JoinStrings(where_conditions, " AND ") + ");";
}

//! Extracts table names from GET operators
static vector<string> GetRelationSet(LogicalOperator &plan) {
	vector<string> relations;
	if (plan.type == LogicalOperatorType::GET) {
		auto *get = static_cast<LogicalGet *>(&plan);
		relations.push_back(get->table->schema->name + "." + get->table->name);
	} else {
		for (auto &child : plan.children) {
			vector<string> child_indices = GetRelationSet(*child);
			relations.insert(relations.end(), child_indices.begin(), child_indices.end());
		}
	}
	sort(relations.begin(), relations.end());
	return relations;
}

void ReOptimizer::InjectCardinalities(LogicalOperator &plan, string temp_table_name) {
	// TODO: derive information about the cardinality of the intermediate results, implement
	cardinalities[JoinStrings(GetRelationSet(plan), ",")] = GetTable("main", temp_table_name)->storage->cardinality;
}

unique_ptr<LogicalOperator> ReOptimizer::CreateNextPlan(unique_ptr<LogicalOperator> plan) {
	idx_t lowest_cost = 0;
	unique_ptr<LogicalOperator> best_plan;
	for (auto subset : TempTablePowerset()) {
		// TODO: create plan using leaf nodes
		unique_ptr<LogicalOperator> tentative_plan;
		// TODO: figure out if this works (bindings might get fucked, but would rather this than call Optimizer)
		JoinOrderOptimizer optimizer(cardinalities);
		tentative_plan = optimizer.Optimize(move(plan));
		// TODO: create cost function that uses injected cardinalities
		idx_t cost = tentative_plan->ComputeCost();
		if (cost < lowest_cost) {
			lowest_cost = cost;
			best_plan = move(tentative_plan);
		}
	}
	// TODO: perhaps call optimizer here?
	return best_plan;
}

vector<vector<string>> ReOptimizer::TempTablePowerset() {
	idx_t powerset_size = pow(2, temp_tables.size());
	vector<vector<string>> powerset(powerset_size);
	for (idx_t i = 0; i < powerset_size; i++) {
		vector<string> subset;
		for (idx_t j = 0; j < temp_tables.size(); j++) {
			if (i & (1 << j)) {
				subset.push_back(temp_tables[j]);
			}
		}
		powerset.push_back(subset);
	}
	return powerset;
}

unique_ptr<LogicalOperator> ReOptimizer::AdjustPlan(unique_ptr<LogicalOperator> plan, LogicalComparisonJoin &old_op,
                                                    const string temporary_table_name) {
	TableCatalogEntry *table = GetTable("main", temporary_table_name);
	// Create a LogicalGet for the newly made temporary table (empty column_ids - filled in "ReplaceLogicalOperator")
	unique_ptr<LogicalGet> temp_table_get = make_unique<LogicalGet>(table, 0, vector<column_t>());
	// replace 'subplan' with 'temp_table_get' in 'plan'
	rebind_mapping = {}; // reset
	ReplaceLogicalOperator(*plan, old_op, table);
	return plan;
}

TableCatalogEntry *ReOptimizer::GetTable(string schema, string table_name) {
	// Catalog::GetTable can only be called if there is an active transaction - else segfault
	CatalogEntry * entry;
	TableCatalogEntry *table;
	if (!context.transaction.HasActiveTransaction()) {
		context.transaction.BeginTransaction();
		entry = context.catalog.GetEntry(context, CatalogType::TABLE, schema, table_name);
		context.transaction.Commit();
	} else {
		entry = context.catalog.GetEntry(context, CatalogType::TABLE, schema, table_name);
	}
	table = static_cast<TableCatalogEntry *>(entry);
	return table;
}

void ReOptimizer::ReplaceLogicalOperator(LogicalOperator &plan, LogicalComparisonJoin &old_op, TableCatalogEntry *table,
                                         idx_t depth) {
	if (plan.children.empty())
		return;
	// search children
	for (idx_t child_i = 0; child_i < plan.children.size(); child_i++) {
		auto &child = plan.children[child_i];
		if (child->type != LogicalOperatorType::COMPARISON_JOIN)
			continue;
		auto *join = static_cast<LogicalComparisonJoin *>(child.get());

		if (join->ParamsToString() == old_op.ParamsToString()) {
			vector<ColumnBinding> old_op_cbs = GetColumnBindings(*join);
			for (column_t new_index = 0; new_index < old_op_cbs.size(); new_index++)
				rebind_mapping[old_op_cbs[new_index].ToString()] = ColumnBinding(depth, new_index);

			vector<column_t> column_ids;
			for (column_t column_id = 0; column_id < table->columns.size(); column_id++)
				column_ids.push_back(column_id);

			unique_ptr<LogicalGet> replacement_operator = make_unique<LogicalGet>(table, depth, column_ids);
			plan.children[child_i] = move(replacement_operator);
			return;
		}
	}
	// enter recursion until the operator is found
	for (auto &child : plan.children)
		ReplaceLogicalOperator(*child.get(), old_op, table, depth + 1);
}

void ReOptimizer::FixColumnBindings(LogicalOperator &plan) {
	// fix BoundColumRefs in JoinConditions
	if (plan.type == LogicalOperatorType::COMPARISON_JOIN) {
		LogicalComparisonJoin *join = static_cast<LogicalComparisonJoin *>(&plan);
		for (idx_t condition_i = 0; condition_i < join->conditions.size(); condition_i++) {
			JoinCondition &jc = join->conditions[condition_i];
			auto l = (BoundColumnRefExpression &)*jc.left.get();
			auto r = (BoundColumnRefExpression &)*jc.right.get();
			if (rebind_mapping.find(l.ToString()) != rebind_mapping.end()) {
				unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
				    l.alias, l.return_type, rebind_mapping[l.ToString()], l.depth);
				join->conditions[condition_i].left = move(fixed_bcre);
			}
			if (rebind_mapping.find(r.ToString()) != rebind_mapping.end()) {
				unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
				    r.alias, r.return_type, rebind_mapping[r.ToString()], r.depth);
				join->conditions[condition_i].right = move(fixed_bcre);
			}
		}
	}
	// fix other BoundColumnRefExpressions found in expressions (for e.g. LogicalProjection)
	for (idx_t expr_i = 0; expr_i < plan.expressions.size(); expr_i++) {
		auto &expr = *plan.expressions[expr_i];
		if (expr.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
			continue;
		auto e = (BoundColumnRefExpression &)expr;
		if (rebind_mapping.find(e.ToString()) != rebind_mapping.end()) {
			unique_ptr<BoundColumnRefExpression> fixed_bcre =
			    make_unique<BoundColumnRefExpression>(e.alias, e.return_type, rebind_mapping[e.ToString()], e.depth);
			plan.expressions[expr_i] = move(fixed_bcre);
		}
	}
	// recursively propagate operator tree
	for (auto &child : plan.children) {
		FixColumnBindings(*child);
	}
}

// FIXME: disabling optimizer indeed causes less overhead (ONLY TESTED FOR SMALL QUERIES), but it's not very clean
// we already have an optimized subquery plan at this point
// slapping a create table and project on top of it should be doable... right?
// also don't know if not having filter pushdown sucks - join order should be fixed tho (by how query is made)

// FIXME: getting cardinality from data_storage only works on autocommit, assume no transactions for now
// this might be because of the fact that cardinality is atomic?
void ReOptimizer::ExecuteSubQuery(const string subquery) {
	// store state of autocommit, profiler and optimizer
	// bool auto_commit = context.transaction.IsAutoCommit();
	context.transaction.Commit();
	bool profiler = context.profiler.IsEnabled();
	// bool optimizer = context.enable_optimizer;
	// disable them
	// context.transaction.SetAutoCommit(false);
	context.profiler.Disable();
	// context.enable_optimizer = false;
	// hack in a query - without autocommit, profiling, optimizer or lock
	// this prevents segfault / assertion fail / unnecessary overhead / deadlock respectively
	context.QueryWithoutLock(subquery, false);
	// restore the state of autocommit, profiler and optimizer
	// context.transaction.SetAutoCommit(auto_commit);
	context.transaction.BeginTransaction();
	if (profiler)
		context.profiler.Enable();
	// context.enable_optimizer = optimizer;
}

unique_ptr<LogicalOperator> ReOptimizer::CallOptimizer(unique_ptr<LogicalOperator> plan) {
	// not interested in optimizer sub-phases
	bool enabled = context.profiler.IsEnabled();
	context.profiler.Disable();
	Optimizer optimizer(binder, context);
	plan = optimizer.Optimize(move(plan));
	if (enabled)
		context.profiler.Enable();
	return plan;
}

unique_ptr<LogicalOperator> ReOptimizer::ClearLeftProjectionMaps(unique_ptr<LogicalOperator> plan) {
	if (plan->children.empty())
		return plan;
	if (plan->type == LogicalOperatorType::COMPARISON_JOIN) {
		auto *join = static_cast<LogicalComparisonJoin *>(plan.get());
		join->left_projection_map.clear();
	}
	for (idx_t i = 0; i < plan->children.size(); i++) {
		plan->children[i] = ClearLeftProjectionMaps(move(plan->children[i]));
	}
	return plan;
}

string ReOptimizer::JoinStrings(vector<string> strings, string delimiter) {
	string joined_strings = "";
	for (idx_t i = 0; i < strings.size(); i++) {
		joined_strings += strings[i];
		if (i < strings.size() - 1)
			joined_strings += delimiter;
	}
	return joined_strings;
}
