#include "duckdb/reoptimizer/reoptimizer.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/unordered_map.hpp"
#include "duckdb/common/enums/expression_type.hpp"
#include "duckdb/common/enums/logical_operator_type.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_join.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"

using namespace std;
using namespace duckdb;

// TODO: write tests that figure out if this shit actually works
ReOptimizer::ReOptimizer(ClientContext &context, Binder &binder) : context(context), binder(binder) {
}

unique_ptr<LogicalOperator> ReOptimizer::ReOptimize(unique_ptr<LogicalOperator> plan, const string query) {
	const string tablename_prefix = "_reopt_temp_" + to_string(hash<string>{}(query));
	for (int i = 0; true; i++) {
		const string temp_table_name = tablename_prefix + "_" + to_string(i);
		plan = SubQuery(move(plan), temp_table_name);
		// if 0 or 1 joins remain the plan can be finished normally
		if (remaining_joins <= 1)
			break;
		context.profiler.StartPhase("optimizer");
		bool enabled = context.profiler.IsEnabled();
		context.profiler.Disable(); // not interested in optimizer sub-phases
		Optimizer optimizer(binder, context);
		plan = optimizer.Optimize(move(plan));
		if (enabled)
			context.profiler.Enable();
		context.profiler.EndPhase();
	}
	context.profiler.StartPhase("clear_left_proj_maps");
	plan = ClearLeftProjectionMaps(move(plan));
	context.profiler.EndPhase();
	return plan;
}

unique_ptr<LogicalOperator> ReOptimizer::SubQuery(unique_ptr<LogicalOperator> plan,
                                                        const string temporary_table_name) {
	vector<LogicalOperator *> joins = ExtractJoinOperators(*plan);
	remaining_joins = joins.size();
	if (remaining_joins <= 1)
		return plan;

	auto *first_join = static_cast<LogicalComparisonJoin *>(joins.back());

	plan->Print();
	// Printer::Print("----------------------------- before");
	// plan->children[0]->GetColumnBindings();
	// string bcres0 = "PROJECTION ";
	// for (size_t i = 0; i < plan->children[0]->expressions.size(); i++) {
	// 	auto &e = *plan->children[0]->expressions[i];
	// 	if (e.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
	// 		BoundColumnRefExpression bcre = (BoundColumnRefExpression &) e;
	// 		bcres0 += bcre.ToString();
	// 	}
	// }
	// Printer::Print(bcres0);
	// plan->children[0]->children[0]->GetColumnBindings();
	// Printer::Print("-----------------------------");

	context.profiler.StartPhase("generate_projection_maps");
	plan = GenerateProjectionMaps(move(plan));
	context.profiler.EndPhase();

	context.profiler.StartPhase("map_binding_names");
	CreateBindingNameMapping(*plan);
	context.profiler.EndPhase();

	context.profiler.StartPhase("create_subquery");
	string subquery = CreateStepQuery(*first_join, temporary_table_name);
	context.profiler.EndPhase();

	Printer::Print(subquery);

	context.profiler.StartPhase("execute_subquery");
	// store state of autocommit and profiler
	bool auto_commit = context.transaction.IsAutoCommit();
	bool profiler_enabled = context.profiler.IsEnabled();
	// disable both
	context.transaction.SetAutoCommit(false);
	context.profiler.Disable();
	// hack in a query - without autocommit, profiling or lock (else segfault / assertion fail / deadlock)
	context.QueryWithoutLock(subquery, false);
	// restore the state of autocommit and profiler
	if (profiler_enabled)
		context.profiler.Enable();
	context.transaction.SetAutoCommit(auto_commit);
	context.profiler.EndPhase();

	context.profiler.StartPhase("adjust_plan");
	plan = AdjustPlan(move(plan), *first_join, temporary_table_name);
	context.profiler.EndPhase();

	context.profiler.StartPhase("fix_bindings");
	FixColumnBindings(*plan);
	context.profiler.EndPhase();

	plan->Print();

	return plan;
}

vector<LogicalOperator *> ReOptimizer::ExtractJoinOperators(LogicalOperator &plan) {
	// this selects the deepest rightmost join operator - DuckDB only uses right-heavy
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
		for (size_t i = 0; i < plan->expressions.size(); i++) {
			auto &expr = *plan->expressions[i];
			if (expr.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
				continue;
			auto bcre = (BoundColumnRefExpression &)expr;
			column_bindings.push_back(bcre.binding);
		}
		// get the child join, modify its left projection map
		auto *child_join = static_cast<LogicalComparisonJoin *>(plan->children[0].get());
		auto left_cbs = child_join->children[0]->GetColumnBindings();
		// add index of needed bindings
		for (column_t cb_index = 0; cb_index < left_cbs.size(); cb_index++) {
			for (auto cb : column_bindings) {
				if (left_cbs[cb_index] == cb) {
					child_join->left_projection_map.push_back(cb_index);
				}
			}
		}
	} else if (plan->type == LogicalOperatorType::COMPARISON_JOIN) {
		auto *join = static_cast<LogicalComparisonJoin *>(plan.get());
		for (size_t i = 0; i < join->conditions.size(); i++) {
			// from conditions
			JoinCondition &jc = join->conditions[i];
			auto l = (BoundColumnRefExpression &)*jc.left.get();
			auto r = (BoundColumnRefExpression &)*jc.right.get();
			column_bindings.push_back(l.binding);
			column_bindings.push_back(r.binding);
			// from bindings of this plan
			for (auto new_cb : plan->GetColumnBindings()) {
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
		}
		// add index of needed bindings, and keep track of removed columns
		auto *child_join = static_cast<LogicalComparisonJoin *>(plan->children[1].get());
		vector<ColumnBinding> left_cbs = LogicalOperator::MapBindings(child_join->children[0]->GetColumnBindings(), child_join->left_projection_map);
		vector<column_t> removed_columns;
		for (column_t cb_index = 0; cb_index < left_cbs.size(); cb_index++) {
			bool keep = false;
			for (auto cb : column_bindings) {
				if (left_cbs[cb_index] == cb) {
					keep = true;
					break;
				}
			}
			if (keep) {
				child_join->left_projection_map.push_back(cb_index);
			} else {
				removed_columns.push_back(cb_index);
			}
		}
		// update right projection map accordingly
		for (size_t rpi = 0; rpi < join->right_projection_map.size(); rpi++) {
			column_t decr = 0;
			for (column_t removed_col : removed_columns) {
				if (removed_col < join->right_projection_map[rpi]) {
					decr ++;
				}
			}
			join->right_projection_map[rpi] -= decr;
		}
	}
	// recursively propagate operator tree
	for (size_t child_i = 0; child_i < plan->children.size(); child_i++) {
		plan->children[child_i] = GenerateProjectionMaps(move(plan->children[child_i]));
	}
	return plan;
}

void ReOptimizer::CreateBindingNameMapping(LogicalOperator &plan) {
	if (plan.children.empty())
		return;
	// find bindings in JoinConditions
	if (plan.type == LogicalOperatorType::COMPARISON_JOIN) {
		auto *join = static_cast<LogicalComparisonJoin *>(&plan);
		for (size_t condition_index = 0; condition_index < join->conditions.size(); condition_index++) {
			JoinCondition &join_condition = join->conditions[condition_index];
			auto l = (BoundColumnRefExpression &)*join_condition.left.get();
			auto r = (BoundColumnRefExpression &)*join_condition.right.get();
			binding_name_mapping[l.binding.ToString()] = l.alias;
			binding_name_mapping[r.binding.ToString()] = r.alias;
		}
	} else {
		// find in expressions (e.g. for LogicalProjection)
		for (size_t expr_index = 0; expr_index < plan.expressions.size(); expr_index++) {
			auto &expr = *plan.expressions[expr_index];
			if (expr.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
				continue;
			auto bcre = (BoundColumnRefExpression &)expr;
			binding_name_mapping[bcre.binding.ToString()] = bcre.alias;
		}
	}
	// recursively propagate operator tree
	for (auto &child : plan.children) {
		CreateBindingNameMapping(*child);
	}
}

string ReOptimizer::CreateStepQuery(LogicalComparisonJoin &join, const string temporary_table_name) {
	vector<string> queried_tables;
	vector<string> queried_columns;
	vector<string> where_conditions;
	for (size_t child_index = 0; child_index < join.children.size(); child_index++) {
		auto &child = join.children.at(child_index);
		string table_alias = "t" + to_string(child_index);
		LogicalGet *logical_get;
		switch (child->type) {
		case LogicalOperatorType::GET: {
			logical_get = static_cast<LogicalGet *>(child.get());
			break;
		}
		case LogicalOperatorType::FILTER: { // if filter we need its child GET
			assert(child->children.size() == 1);
			auto *logical_filter = static_cast<LogicalFilter *>(child.get());
			for (auto &expr : logical_filter->expressions) {
				if (expr->GetExpressionClass() != ExpressionClass::BOUND_COMPARISON)
					throw new ReOptimizerException("Expected filter class to be BOUND_COMPARISON, got '%s' instead",
					                               expr->GetExpressionClass());
				auto *bce = static_cast<BoundComparisonExpression *>(expr.get());
				// filter conditions
				// table reference is 'always' placed on the left, constant on the right
				if (bce->right->return_type == TypeId::VARCHAR) // strings need to be escaped with quotes
					where_conditions.push_back(table_alias + "." + bce->left->alias +
					                           ExpressionTypeToOperator(bce->type) + "'" + bce->right->ToString() +
					                           "'");
				else
					where_conditions.push_back(table_alias + "." + expr->GetName());
			}
			logical_get = static_cast<LogicalGet *>(logical_filter->children.at(0).get());
			break;
		}
		default:
			throw new ReOptimizerException("Expected child of join to be GET or FILTER, got '%s' instead", child->type);
		}

		string schema = logical_get->table->schema->name;
		string table_name = logical_get->table->name;
		queried_tables.push_back(schema + "." + table_name + " AS " + table_alias);

		vector<ColumnBinding> child_bindings = logical_get->GetColumnBindings();
		vector<column_t> projection_map = child_index == 0 ? join.left_projection_map : join.right_projection_map;
		// take all bindings if empty
		if (projection_map.empty()) {
			for (size_t proj_index = 0; proj_index < child_bindings.size(); proj_index++) {
				projection_map.push_back(proj_index);
			}
		}

		// projection map stores indices of bindings, the bindings have column indices
		for (column_t binding_index : projection_map) {
			string col_name = binding_name_mapping[child_bindings[binding_index].ToString()];
			queried_columns.push_back(table_alias + "." + col_name + " AS " + col_name);
		}
	}
	// join conditions
	for (auto &cond : join.conditions) {
		where_conditions.push_back("t0." + cond.left->GetName() + " = t1." + cond.right->GetName());
	}
	// create query and return
	return "CREATE TEMPORARY TABLE main." + temporary_table_name + " AS (" + "SELECT " +
	       JoinStrings(queried_columns, ", ") + " " + "FROM " + JoinStrings(queried_tables, ", ") + " " + "WHERE " +
	       JoinStrings(where_conditions, " AND ") + ");";
}

unique_ptr<LogicalOperator> ReOptimizer::AdjustPlan(unique_ptr<LogicalOperator> plan, LogicalComparisonJoin &step,
                                                    const string temporary_table_name) {
	TableCatalogEntry *table = GetTable("main", temporary_table_name);

	// Create a LogicalGet for the newly made temporary table (empty column_ids - filled in "ReplaceLogicalOperator")
	unique_ptr<LogicalGet> temp_table_get = make_unique<LogicalGet>(table, 0, vector<column_t>());

	// replace 'step' with 'temp_table_get' in 'plan'
	new_bindings_mapping = {}; // reset
	ReplaceLogicalOperator(*plan, step, table);
	remaining_joins--;
	return plan;
}

TableCatalogEntry *ReOptimizer::GetTable(string schema, string table_name) {
	// Catalog::GetTable can only be called if there is an active transaction - else segfault
	TableCatalogEntry *table;
	if (!context.transaction.HasActiveTransaction()) {
		context.transaction.BeginTransaction();
		table = context.catalog.GetTable(context, schema, table_name);
		context.transaction.Commit();
	} else {
		table = context.catalog.GetTable(context, schema, table_name);
	}
	return table;
}

void ReOptimizer::ReplaceLogicalOperator(LogicalOperator &plan, LogicalComparisonJoin &old_op, TableCatalogEntry *table,
                                         index_t depth) {
	if (plan.children.empty())
		return;
	// search children
	for (size_t child_index = 0; child_index < plan.children.size(); child_index++) {
		auto &child = plan.children.at(child_index);
		if (child->type != LogicalOperatorType::COMPARISON_JOIN)
			continue;
		auto *join = static_cast<LogicalComparisonJoin *>(child.get());

		if (join->ParamsToString() == old_op.ParamsToString()) {
			for (auto &child : old_op.children) {
				for (auto child_cb : child->GetColumnBindings()) {
					new_bindings_mapping[child_cb.ToString()] = ColumnBinding(depth, child_cb.column_index);
				}
			}
			vector<column_t> column_ids;
			for (column_t column_id = 0; column_id < old_op.GetColumnBindings().size(); column_id++)
				column_ids.push_back(column_id);

			unique_ptr<LogicalGet> replacement_operator = make_unique<LogicalGet>(table, depth, column_ids);
			plan.children.at(child_index) = move(replacement_operator);
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
		for (size_t condition_index = 0; condition_index < join->conditions.size(); condition_index++) {
			JoinCondition &jc = join->conditions[condition_index];
			auto l = (BoundColumnRefExpression &)*jc.left.get();
			auto r = (BoundColumnRefExpression &)*jc.right.get();
			if (new_bindings_mapping.find(l.ToString()) != new_bindings_mapping.end()) {
				unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
				    l.alias, l.return_type, new_bindings_mapping[l.ToString()], l.depth);
				join->conditions.at(condition_index).left = move(fixed_bcre);
			}
			if (new_bindings_mapping.find(r.ToString()) != new_bindings_mapping.end()) {
				unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
				    r.alias, r.return_type, new_bindings_mapping[r.ToString()], r.depth);
				join->conditions.at(condition_index).right = move(fixed_bcre);
			}
		}
	}
	// fix other BoundColumnRefExpressions found in expressions (for e.g. LogicalProjection)
	for (size_t expr_index = 0; expr_index < plan.expressions.size(); expr_index++) {
		auto &expr = *plan.expressions[expr_index];
		if (expr.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
			continue;
		auto e = (BoundColumnRefExpression &)expr;
		if (new_bindings_mapping.find(e.ToString()) != new_bindings_mapping.end()) {
			unique_ptr<BoundColumnRefExpression> fixed_bcre =
			    make_unique<BoundColumnRefExpression>(e.alias, e.return_type, new_bindings_mapping[e.ToString()], e.depth);
			plan.expressions.at(expr_index) = move(fixed_bcre);
		}
	}
	// recursively propagate operator tree
	for (auto &child : plan.children) {
		FixColumnBindings(*child);
	}
}

unique_ptr<LogicalOperator> ReOptimizer::ClearLeftProjectionMaps(unique_ptr<LogicalOperator> plan) {
	if (plan->children.empty())
		return plan;
	if (plan->type == LogicalOperatorType::COMPARISON_JOIN) {
		auto *join = static_cast<LogicalComparisonJoin *>(plan.get());
		join->left_projection_map.clear();
	}
	for (size_t i = 0; i < plan->children.size(); i++) {
		plan->children[i] = ClearLeftProjectionMaps(move(plan->children[i]));
	}
	return plan;
}

string ReOptimizer::JoinStrings(vector<string> strings, string delimiter) {
	string joined_strings = "";
	for (size_t i = 0; i < strings.size(); i++) {
		joined_strings += strings.at(i);
		if (i < strings.size() - 1)
			joined_strings += delimiter;
	}
	return joined_strings;
}
