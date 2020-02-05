#include "duckdb/reoptimizer/reoptimizer.hpp"

#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/unordered_map.hpp"
#include "duckdb/common/enums/expression_type.hpp"
#include "duckdb/common/enums/logical_operator_type.hpp"
#include "duckdb/main/table_description.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/parser/expression/columnref_expression.hpp"
#include "duckdb/parser/statement/prepare_statement.hpp"
#include "duckdb/planner/joinside.hpp"
#include "duckdb/planner/planner.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_join.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"
#include "duckdb/transaction/transaction.hpp"

using namespace std;
using namespace duckdb;

// TODO: write tests that figure out if this shit actually works
// FIXME: avoid use of unique_ptr::get - this copies the whole thang
/**
 * One big fucking TODO:
 * Need to actually set the value of left_projection_map in every LogicalComparisonJoin in the plan
 * Already know how to compute this, just need to write the code for it
 * 
 * Editing the left_projection_map will screw with the projection maps of the parents
 * - those will need to be changed as well
 * 
 * How the fuck do we begin? - might be too difficult to pull off
 * 
 * 
 * 
 * 
 * Another approach could be to re-edit the plan after calling Optimize to re-add the column that it optimizes away
 * This will lead to the same issues as with the other approach - projection maps will be fucked, AND it will be slower
 *  */
ReOptimizer::ReOptimizer(ClientContext &context, Binder &binder) : context(context), binder(binder) {
}

unique_ptr<LogicalOperator> ReOptimizer::CreateStepPlan(unique_ptr<LogicalOperator> plan,
                                                        const string temporary_table_name) {
	vector<LogicalOperator *> joins = ExtractJoinOperators(*plan);
	remaining_joins = joins.size();
	if (remaining_joins <= 1)
		return plan;

	plan->Print();
	Printer::Print("----------------------------- before");
	plan->children[0]->GetColumnBindings();
	string bcres = "PROJECTION ";
	for (size_t i = 0; i < plan->children[0]->expressions.size(); i++) {
		auto &e = *plan->children[0]->expressions[i];
		if (e.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
			BoundColumnRefExpression bcre = (BoundColumnRefExpression &) e;
			bcres += bcre.ToString();
		}
	}
	Printer::Print(bcres);
	plan->children[0]->children[0]->GetColumnBindings();
	Printer::Print("-----------------------------");
	// binder.bind_context.Print();

	LogicalComparisonJoin *first_join = static_cast<LogicalComparisonJoin *>(joins.back());

	context.profiler.StartPhase("create_step");
	CreateBindingNameMapping(*plan);
	string step_query = CreateStepQuery(*first_join, temporary_table_name);
	context.profiler.EndPhase();

	Printer::Print(step_query);

	context.profiler.StartPhase("execute_step");
	// ClientContext::Query will commit if we are on autocommit - not desirable since we are technically mid-query
	// store whether we are on autocommit so we can restore it after executing our subquery
	bool auto_commit = context.transaction.IsAutoCommit();
	context.transaction.SetAutoCommit(false);
	// the profiler will be confused by the double query - disable and restore just like autocommit
	bool profiler_enabled = context.profiler.IsEnabled();
	context.profiler.Disable();
	// we now execute the step: this is tricky since there are safeguards preventing concurrent queries
	context.QueryWithoutLock(step_query, false);
	if (profiler_enabled)
		context.profiler.Enable();
	context.transaction.SetAutoCommit(auto_commit);
	context.profiler.EndPhase();

	context.profiler.StartPhase("edit_plan");
	plan = AdjustPlan(move(plan), *first_join, temporary_table_name);
	assert(plan);
	remaining_joins--;
	context.profiler.EndPhase();

	// Printer::Print("##########");
	// for (unordered_map<string, ColumnBinding>::iterator it = bindings_mapping.begin(); it != bindings_mapping.end();
	//      it++) {
	// 	Printer::Print(it->first + " -> " + it->second.ToString());
	// }
	// Printer::Print("##########");

	context.profiler.StartPhase("fix_bindings");
	FixColumnBindings(*plan);
	context.profiler.EndPhase();

	plan->Print();
	Printer::Print("----------------------------- after");
	plan->children[0]->GetColumnBindings();
	string bcres2 = "PROJECTION ";
	for (size_t i = 0; i < plan->children[0]->expressions.size(); i++) {
		auto &e = *plan->children[0]->expressions[i];
		if (e.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
			BoundColumnRefExpression bcre = (BoundColumnRefExpression &) e;
			bcres2 += bcre.ToString();
		}
	}
	Printer::Print(bcres2);
	plan->children[0]->children[0]->GetColumnBindings();
	Printer::Print("-----------------------------");
	// binder.bind_context.Print();

	return plan;
}

void ReOptimizer::CreateBindingNameMapping(LogicalOperator &plan) {
	if (plan.children.empty())
		return;
	// find bindings in JoinConditions
	if (plan.type == LogicalOperatorType::COMPARISON_JOIN) {
		LogicalComparisonJoin *join = static_cast<LogicalComparisonJoin *>(&plan);
		for (size_t condition_index = 0; condition_index < join->conditions.size(); condition_index++) {
			JoinCondition &join_condition = join->conditions[condition_index];
			BoundColumnRefExpression l = (BoundColumnRefExpression &)*join_condition.left.get();
			BoundColumnRefExpression r = (BoundColumnRefExpression &)*join_condition.right.get();
			binding_name_mapping[l.binding.ToString()] = l.alias;
			binding_name_mapping[r.binding.ToString()] = r.alias;
		}
	} else {
		// find in expressions (e.g. for LogicalProjection)
		for (size_t expr_index = 0; expr_index < plan.expressions.size(); expr_index++) {
			auto &expr = *plan.expressions[expr_index];
			if (expr.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
				continue;
			BoundColumnRefExpression e = (BoundColumnRefExpression &)expr;
			binding_name_mapping[e.binding.ToString()] = e.alias;
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
			LogicalFilter *logical_filter = static_cast<LogicalFilter *>(child.get());
			for (auto &expr : logical_filter->expressions) {
				if (expr->GetExpressionClass() != ExpressionClass::BOUND_COMPARISON)
					throw new ReOptimizerException("Expected filter class to be BOUND_COMPARISON, got '%s' instead",
					                               expr->GetExpressionClass());
				BoundComparisonExpression *bce = static_cast<BoundComparisonExpression *>(expr.get());
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
		// Set queried table name values for use later (FIXME: might not be needed - see TODO in ReplaceLogicalOperator)
		// if (i == 0) left_table_name = logical_get->table->name;
		// if (i == 1) right_table_name = logical_get->table->name;

		string schema = logical_get->table->schema->name;
		string table_name = logical_get->table->name;
		queried_tables.push_back(schema + "." + table_name + " AS " + table_alias);

		// TODO: working on this - big shit
		// TableCatalogEntry *table = GetTable(schema, table_name);

		// FIXME: left_projection_map is empty - projecting all columns - some of which are removed later
		vector<ColumnBinding> child_bindings = logical_get->GetColumnBindings();
		vector<column_t> projection_map = child_index == 0 ? join.left_projection_map : join.right_projection_map;
		// take all bindings if empty FIXME: probably not needed anymore
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
	return plan;
}

// FIXME: we already know which operator to replace, but we need to replace it using the parent
// e.g.: parent.children.at(i) = move(new_op). However, we don't have the parent
void ReOptimizer::ReplaceLogicalOperator(LogicalOperator &plan, LogicalComparisonJoin &old_op, TableCatalogEntry *table,
                                         index_t depth) {
	if (plan.children.empty())
		return;
	// search children
	for (size_t child_index = 0; child_index < plan.children.size(); child_index++) {
		auto &child = plan.children.at(child_index);
		if (child->type != LogicalOperatorType::COMPARISON_JOIN)
			continue;
		LogicalComparisonJoin *join = static_cast<LogicalComparisonJoin *>(child.get());

		if (join->ParamsToString() == old_op.ParamsToString()) {
			// TODO: find out if these two lines are needed (seems like NO, maybe when Optimize again?)
			// binder.bind_context.ReplaceBindingIndex(binder.bind_context.GetBindingAlias(left_table_name),
			// depth);
			// binder.bind_context.ReplaceBindingIndex(binder.bind_context.GetBindingAlias(right_table_name),
			// depth);

			// create mapping to account for new depth of bindings - FIXME: creates too many column bindings because no left map
			// FIXME: already fixed this in query creation, but we did not set the map there (perhaps do this anyway? - need unique_ptr)
			// also need to keep own index instead of taking child column index, since we shove everything to the left when one gets removed
			// -----------------------------------------------------------------------------------------------
			// vector<column_t> column_ids;
			// int column_index = 0;
			// for (size_t i = 0; i < old_op.children.size(); i++) {
			// 	auto &child = old_op.children[i];
			// 	vector<column_t> projection_map = i == 0 ? left_projection_map : old_op.right_projection_map;
			// 	vector<ColumnBinding> child_bindings = child->GetColumnBindings();
			// 	for (column_t col : projection_map) {
			// 		new_bindings_mapping[child_bindings[col].ToString()] = ColumnBinding(depth, column_index);
			// 		column_ids.push_back(column_index);
			// 		column_index++;
			// 	}
			// }
			// -----------------------------------------------------------------------------------------------

			// int column_index = 0; ????? FIXME:
			for (auto &child : old_op.children) {
				for (ColumnBinding child_cb : child->GetColumnBindings()) {
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

// FIXME:123 the projection map must be changed to reflect the changes in bindings as well
// max(projection_map) must be < than length of bindings
void ReOptimizer::FixColumnBindings(LogicalOperator &plan) {
	// fix BoundColumRefs in JoinConditions
	if (plan.type == LogicalOperatorType::COMPARISON_JOIN) {
		LogicalComparisonJoin *join = static_cast<LogicalComparisonJoin *>(&plan);
		for (size_t condition_index = 0; condition_index < join->conditions.size(); condition_index++) {
			JoinCondition &jc = join->conditions[condition_index];
			BoundColumnRefExpression l = (BoundColumnRefExpression &)*jc.left.get();
			BoundColumnRefExpression r = (BoundColumnRefExpression &)*jc.right.get();
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
		// TODO: Might want to use LogicalProjection::table_index ??
		auto &expr = *plan.expressions[expr_index];
		if (expr.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
			continue;
		BoundColumnRefExpression e = (BoundColumnRefExpression &)expr;
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

vector<LogicalOperator *> ReOptimizer::ExtractJoinOperators(LogicalOperator &plan) {
	// TODO: this selects the deepest rightmost join operator - prefer a better selection method
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

std::set<TableCatalogEntry *> ReOptimizer::ExtractGetTables(LogicalOperator &plan) {
	std::set<TableCatalogEntry *> gets;
	if (plan.type == LogicalOperatorType::GET) {
		LogicalGet *lg = static_cast<LogicalGet *>(&plan);
		gets.insert(lg->table);
	}
	for (auto &child : plan.children) {
		std::set<TableCatalogEntry *> child_gets = ExtractGetTables(*child);
		for (std::set<TableCatalogEntry *>::iterator it = child_gets.begin(); it != child_gets.end(); it++)
			gets.insert(*it);
	}
	return gets;
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

string ReOptimizer::JoinStrings(vector<string> strings, string delimiter) {
	string joined_strings = "";
	for (size_t i = 0; i < strings.size(); i++) {
		joined_strings += strings.at(i);
		if (i < strings.size() - 1)
			joined_strings += delimiter;
	}
	return joined_strings;
}
