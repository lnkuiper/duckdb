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
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_join.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"

using namespace std;
using namespace duckdb;

ReOptimizer::ReOptimizer(ClientContext &context, LogicalOperator &plan, Binder &binder) : context(context), binder(binder) {
}

unique_ptr<LogicalOperator> ReOptimizer::CreateStepPlan(unique_ptr<LogicalOperator> plan,
                                                        const string temporary_table_name) {
	vector<LogicalOperator *> joins = ExtractJoinOperators(*plan);
	if (joins.empty())
		return plan;

	plan->Print();

	used_columns_per_table = {}; // reset
	ExtractUsedColumns(*plan);

	LogicalComparisonJoin *first_join = static_cast<LogicalComparisonJoin *>(joins.back());

	context.profiler.StartPhase("create_step");
	string step_query = CreateStepQuery(*first_join, temporary_table_name);
	context.profiler.EndPhase();

	Printer::Print(step_query);

	context.profiler.StartPhase("execute_step");
	context.Query(step_query, false);
	context.profiler.EndPhase();

	/* TODO: edit 'plan', inserting the temporary table for the join node
	 * this will probably be quite hard, because of the binder
	 * might need to Skype with Mark/Hannes
	 * 
	 * UPDATE: it seems to work! needs more testing */

	context.profiler.StartPhase("edit_plan");
	plan = AdjustPlan(move(plan), *first_join, temporary_table_name);
	FixColumnBindings(*plan);
	context.profiler.EndPhase();

	plan->Print();

	return plan;
}

string ReOptimizer::CreateStepQuery(LogicalComparisonJoin &plan, const string temporary_table_name) {
	vector<string> queried_tables;
	vector<string> queried_columns;
	vector<string> where_conditions;
	assert(plan.children.size() == 2);
	for (int i = 0; i < plan.children.size(); i++) {
		auto &child = plan.children.at(i);
		string table_alias = "t" + to_string(i);
		LogicalGet *logical_get;
		switch (child->type) {
		case LogicalOperatorType::GET: {
			logical_get = static_cast<LogicalGet *>(child.get());
			break;
		}
		case LogicalOperatorType::FILTER: { // with filter we need the child GET
			assert(child->children.size() == 1);
			LogicalFilter *logical_filter = static_cast<LogicalFilter *>(child.get());
			for (auto &expr : logical_filter->expressions) {
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

		queried_tables.push_back(logical_get->table->schema->name + "." + logical_get->table->name + " AS " + table_alias);
		for (string col : used_columns_per_table[logical_get->table->name]) {
			queried_columns.push_back(table_alias + "." + col + " AS " + col);
		}
	}
	for (auto &cond : plan.conditions) {
		where_conditions.push_back("t0." + cond.left->GetName() + " = t1." + cond.right->GetName());
	}
	
	return "CREATE TEMPORARY TABLE " + temporary_table_name + " AS (" + "SELECT " + JoinStrings(queried_columns, ", ") + " " +
	       "FROM " + JoinStrings(queried_tables, ", ") + " " + "WHERE " + JoinStrings(where_conditions, " AND ") + ");";
}

unique_ptr<LogicalOperator> ReOptimizer::AdjustPlan(unique_ptr<LogicalOperator> plan, LogicalComparisonJoin &step,
                                                    const string temporary_table_name) {
	// FIXME: derive schema name instead of using "main" by default?
	TableCatalogEntry *table = context.catalog.GetTable(context, "main", temporary_table_name);	

	// Create a LogicalGet for the newly made temporary table (empty column_ids - filled in "ReplaceLogicalOperator")
	unique_ptr<LogicalGet> temp_table_get = make_unique<LogicalGet>(table, 0, vector<column_t>());

	// replace 'step' with 'temp_table_get' in 'plan'
	bindings_mapping = {}; // reset
	ReplaceLogicalOperator(*plan, step, move(temp_table_get));
	return plan;
}

// FIXME: check if we can simply store a pointer to &old_op (e.g. plan.children.at(i) or something and replace without recursion)
void ReOptimizer::ReplaceLogicalOperator(LogicalOperator &plan, LogicalComparisonJoin &old_op,
                                     unique_ptr<LogicalGet> new_op) {
	if (plan.children.empty())
		return;

	// recursively find the step to replace
	bool found = false;
	for (int i = 0; i < plan.children.size(); i++) {
        new_op->table_index++; // increase by depth
		auto &child = plan.children.at(i);
		if (child->type != LogicalOperatorType::COMPARISON_JOIN)
			continue;
		LogicalComparisonJoin *join = static_cast<LogicalComparisonJoin *>(child.get());
		bool equal = true;
		for (int j = 0; j < join->conditions.size(); j++) {
			if (join->conditions.at(j).left->GetName() != old_op.conditions.at(j).left->GetName()) {
				equal = false;
				break;
			}
		}
		if (equal) {
			found = true;
			// TODO: find out if the following two lines are needed (seems like NO, perhaps when calling Optimize again?)
			// binder.bind_context.ReplaceBindingIndex(binder.bind_context.GetBindingAlias(left_table_name), new_op->table_index);
			// binder.bind_context.ReplaceBindingIndex(binder.bind_context.GetBindingAlias(right_table_name), new_op->table_index);
			for (ColumnBinding cb : old_op.GetColumnBindings()) {
				// Fill with mapping from old bindings to new bindings (visit children too)
				bindings_mapping[cb.ToString()] = ColumnBinding(new_op->table_index, cb.column_index);
				for (auto &child : old_op.children) {
					for (ColumnBinding cbc : child->GetColumnBindings()) {
						bindings_mapping[cbc.ToString()] = ColumnBinding(new_op->table_index, cbc.column_index);
					}
				}
				// Add column id to LogicalGet column_ids (no duplicate column ids)
				if (find(new_op->column_ids.begin(), new_op->column_ids.end(), cb.column_index) == new_op->column_ids.end())
					new_op->column_ids.push_back(cb.column_index);
			}
			plan.children.at(i) = move(new_op);
			break;
		}
	}
	// Enter recursion until the operator in question is found
	if (!found) {
		for (auto &child : plan.children)
			ReplaceLogicalOperator(*child.get(), old_op, move(new_op));
	}
}

void ReOptimizer::FixColumnBindings(LogicalOperator &plan) {
	// Replace references to old bindings with the new bindings
	for (int i = 0; i < plan.expressions.size(); i++) {
		auto &expr = *plan.expressions[i];
		if (expr.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
			continue;
		BoundColumnRefExpression e = (BoundColumnRefExpression &)expr;
		if (bindings_mapping.find(e.ToString()) != bindings_mapping.end()) {
			unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(e.alias, e.return_type, bindings_mapping[e.ToString()], e.depth);
			plan.expressions.at(i) = move(fixed_bcre);
		}
	}
	// Recursively propagate operator tree
	for (auto &child : plan.children) {
		FixColumnBindings(*child);
	}
}

void ReOptimizer::ExtractUsedColumns(LogicalOperator &plan) {
	if (plan.children.empty())
		return;

	// Extract set of column names that are used by 'plan' or its children
	std::set<string> used_columns;
	switch (plan.type) {
	case LogicalOperatorType::PROJECTION: {
		LogicalProjection *lp = static_cast<LogicalProjection *>(&plan);
		for (int i = 0; i < lp->expressions.size(); i++)
			used_columns.insert(lp->expressions.at(i)->GetName());
		break;
	}
	case LogicalOperatorType::COMPARISON_JOIN: {
		LogicalComparisonJoin *lcj = static_cast<LogicalComparisonJoin *>(&plan);
		for (auto &cond : lcj->conditions) {
			used_columns.insert(cond.left->GetName());
			used_columns.insert(cond.right->GetName());
		}
		break;
	}
	}
	// Match column names in the referenced tables
	if (!used_columns.empty()) {
		std::set<TableCatalogEntry *> tables = ExtractGetTables(plan);
		for (std::set<TableCatalogEntry *>::iterator it = tables.begin(); it != tables.end(); it++) {
			TableCatalogEntry *table = *it;
			for (ColumnDefinition &cd : table->columns) {
				if (used_columns.find(cd.name) != used_columns.end()) {
					used_columns_per_table[table->name].insert(cd.name);
				}
			}
		}
	}
	// Enter recursion
	for (auto &child : plan.children) {
		ExtractUsedColumns(*child);
	}
}

vector<LogicalOperator *> ReOptimizer::ExtractJoinOperators(LogicalOperator &plan) {
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

string ReOptimizer::JoinStrings(vector<string> strings, string delimiter) {
	string joined_strings = "";
	for (int i = 0; i < strings.size(); i++) {
		joined_strings += strings.at(i);
		if (i < strings.size() - 1)
			joined_strings += delimiter;
	}
	return joined_strings;
}
