//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/reoptimizer/reoptimizer.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include <set>

#include "duckdb/main/client_context.hpp"
#include "duckdb/planner/binder.hpp"
#include "duckdb/planner/column_binding.hpp"
#include "duckdb/planner/logical_operator.hpp"

namespace duckdb {
class Binder;

class ReOptimizer {
public:
	ReOptimizer(ClientContext &context, LogicalOperator &plan, Binder &binder);

	//! Executes the first join, then adapts the plan accordingly FIXME: rename
	unique_ptr<LogicalOperator> CreateStepPlan(unique_ptr<LogicalOperator> plan, const string temporary_table_name);

	//! Names of the left and right tables of the first join
	string left_table_name;
	string right_table_name;

	//! The client context
	ClientContext &context;
    //! The binder
    Binder &binder;

private:
	//! Creates a CREATE TEMPORARY TABLE query string for the first join to be executed in 'plan'
	string CreateStepQuery(LogicalComparisonJoin &plan, const string temporary_table_name);
	//! Adjusts the original plan by replacing the join with a LogicalGet on the temporary table
	unique_ptr<LogicalOperator> AdjustPlan(unique_ptr<LogicalOperator> plan, LogicalComparisonJoin &step,
	                                       string temporary_table_name);
	//! Replaces the join 'old_op' in 'plan' with the given operator 'new_op'
	void ReplaceLogicalOperator(LogicalOperator &plan, LogicalComparisonJoin &old_op, unique_ptr<LogicalGet> new_op);
	//! Fixes column bindings after replacing JOIN with GET
	void FixColumnBindings(LogicalOperator &plan);

	//! sets mapping of tablename -> set of columns used in 'plan' in 'used_columns_per_table'
	void ExtractUsedColumns(LogicalOperator &plan);
	//! returns all join operators in the plan - the last element is the first one to be executed
	vector<LogicalOperator *> ExtractJoinOperators(LogicalOperator &plan);
	//! Extracts pointers to TableCatalogEntry for each LogicalGet in 'plan'
	std::set<TableCatalogEntry *> ExtractGetTables(LogicalOperator &plan);

	//! Utility to join strings like in Java, Python
	string JoinStrings(vector<string> strings, string delimiter);

	//! mapping of tablename -> set of columns FIXME: rename
	unordered_map<string, std::set<string>> used_columns_per_table;
	//! The new column bindings (after replacing JOIN with GET)
	unordered_map<string, ColumnBinding> bindings_mapping;
};

} // namespace duckdb
