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
	//! First half of the re-optimization iteration procedure: perform subquery and adjust plan
	unique_ptr<LogicalOperator> PerformPartialPlan(unique_ptr<LogicalOperator> plan, const string temporary_table_name);
	//! Decides which join in the plan to execute as subquery
	LogicalOperator *DecideSubQueryPlan(LogicalOperator &plan);
	//! Returns all join operators in the plan - the last element is the first one to be executed
	vector<LogicalOperator *> ExtractJoinOperators(LogicalOperator &plan);
	//! Generate left projection map for joins in the plan, and change right map accordingly
	unique_ptr<LogicalOperator> GenerateProjectionMaps(unique_ptr<LogicalOperator> plan);
	//! Fill binding_name_mapping with (binding -> alias)
	void CreateMaps(LogicalOperator &plan);
	//! Creates a CREATE TEMPORARY TABLE query string for the first join to be executed in 'plan'
	string CreateSubQuery(LogicalComparisonJoin &plan, const string temporary_table_name,
	                      vector<string> &queried_tables, vector<string> &where_conditions);
	//! Adjusts the original plan by replacing the join with a LogicalGet on the temporary table
	unique_ptr<LogicalOperator> AdjustPlan(unique_ptr<LogicalOperator> plan, LogicalComparisonJoin &old_op,
	                                       string temporary_table_name);
	//! Call Catalog::GetTable (works around autocommit stuff)
	TableCatalogEntry *GetTable(string schema, string table_name);
	//! Replaces the join 'old_op' in 'plan' with the given operator 'new_op'
	void ReplaceLogicalOperator(LogicalOperator &plan, LogicalComparisonJoin &old_op, TableCatalogEntry *table,
	                            index_t depth = 3);
	//! Fixes column bindings after replacing JOIN with GET
	void FixColumnBindings(LogicalOperator &plan);
	//! Executes a query in the middle of the re-optimization process
	void ExecuteSubQuery(const string subquery);
	//! Second half of the re-optimization iteration procedure: call Optimizer::Optimize on the adjusted plan
	unique_ptr<LogicalOperator> CallOptimizer(unique_ptr<LogicalOperator> plan);
	//! Empty all left projection maps again (required by PhysicalPlanGenerator - PhysicalHashJoin assert)
	unique_ptr<LogicalOperator> ClearLeftProjectionMaps(unique_ptr<LogicalOperator> plan);

	//! Utility to join strings like in Java, Python
	string JoinStrings(vector<string> strings, string delimiter);

	//! The client context
	ClientContext &context;
	//! The binder
	Binder &binder;

	//! Binding to column alias mapping (created by CreateMaps)
	unordered_map<string, string> bta;
	//! The new column bindings (after replacing JOIN with GET)
	unordered_map<string, ColumnBinding> rebind_mapping;

	//! Whether we are done re-optimizing the plan
	bool done = false;
};

} // namespace duckdb
