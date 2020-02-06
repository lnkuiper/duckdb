//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/reoptimizer/reoptimizer.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/main/client_context.hpp"
#include "duckdb/planner/binder.hpp"
#include "duckdb/planner/column_binding.hpp"
#include "duckdb/planner/logical_operator.hpp"

namespace duckdb {
class Binder;

class ReOptimizer {
public:
	ReOptimizer(ClientContext &context, Binder &binder);

	//! Reoptimization loop until only 1 join remains
	unique_ptr<LogicalOperator> ReOptimize(unique_ptr<LogicalOperator> plan, const string query);

private:
	//! Executes the first join, then adapts the plan accordingly FIXME
	unique_ptr<LogicalOperator> SubQuery(unique_ptr<LogicalOperator> plan, const string temporary_table_name);
	//! Returns all join operators in the plan - the last element is the first one to be executed
	vector<LogicalOperator *> ExtractJoinOperators(LogicalOperator &plan);
	//! Generate left projection map for joins in the plan, and change right map accordingly
	unique_ptr<LogicalOperator> GenerateProjectionMaps(unique_ptr<LogicalOperator> plan);
	//! Fill binding_name_mapping with (binding -> alias)
	void CreateBindingNameMapping(LogicalOperator &plan);
	//! Creates a CREATE TEMPORARY TABLE query string for the first join to be executed in 'plan'
	string CreateSubQuery(LogicalComparisonJoin &plan, const string temporary_table_name);
	//! Adjusts the original plan by replacing the join with a LogicalGet on the temporary table
	unique_ptr<LogicalOperator> AdjustPlan(unique_ptr<LogicalOperator> plan, LogicalComparisonJoin &old_op,
	                                       string temporary_table_name);
	//! Get a (works around autocommit stuff)
	TableCatalogEntry *GetTable(string schema, string table_name);
	//! Replaces the join 'old_op' in 'plan' with the given operator 'new_op'
	void ReplaceLogicalOperator(LogicalOperator &plan, LogicalComparisonJoin &old_op, TableCatalogEntry *table,
	                            index_t depth = 3);
	//! Fixes column bindings after replacing JOIN with GET
	void FixColumnBindings(LogicalOperator &plan);
	//! Empty all left projection maps again (required by PhysicalPlanGenerator - PhysicalHashJoin assert)
	unique_ptr<LogicalOperator> ClearLeftProjectionMaps(unique_ptr<LogicalOperator> plan);

	//! Utility to join strings like in Java, Python
	string JoinStrings(vector<string> strings, string delimiter);

	//! The client context
	ClientContext &context;
	//! The binder
	Binder &binder;

	//! Binding to name mapping (created by CountBindingReferences)
	unordered_map<string, string> binding_to_alias;
	//! The new column bindings (after replacing JOIN with GET)
	unordered_map<string, ColumnBinding> rebind_mapping;

	//! The amount of remaining joins in the plan given to the last call to CreateSubQuery
	int remaining_joins = 0;
};

} // namespace duckdb
