//===----------------------------------------------------------------------===//
//                         DuckDB
//
// reoptimizer/reoptimizer.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/planner/logical_operator.hpp"

namespace duckdb {

class ReOptimizer {
public:
	ReOptimizer(ClientContext &context);

    unique_ptr<LogicalOperator> CreateFirstStepPlan(unique_ptr<LogicalOperator> plan, string table_name);
    unique_ptr<LogicalOperator> remaining_plan;

    ClientContext &context;

    string step_query;

private:
    void SetTemporaryTableQuery(LogicalOperator &plan, string table_name);
    vector<LogicalOperator> GetJoinOperators(LogicalOperator &plan);
    vector<string> ColumnNamesFromLogicalGet(LogicalGet &logical_get);
    string JoinStrings(vector<string> strings, string delimiter);
    string JoinConditionFromLogicalPlan(LogicalOperator &plan, string schema);
};

} // namespace duckdb
