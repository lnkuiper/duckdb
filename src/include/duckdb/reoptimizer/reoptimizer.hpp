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

    string CreateFirstStepQuery(unique_ptr<LogicalOperator> plan, string table_name);
    unique_ptr<LogicalOperator> remaining_plan;

    ClientContext &context;

private:
    string CreateTemporaryTableQuery(unique_ptr<LogicalOperator> plan, string table_name);
    vector<unique_ptr<LogicalOperator>> GetJoinOperators(unique_ptr<LogicalOperator> plan);
    vector<string> ColumnNamesFromLogicalGet(LogicalGet* logical_get);
    string JoinStrings(vector<string> strings, string delimiter);
    string JoinConditionFromLogicalPlan(unique_ptr<LogicalOperator> plan);
};

} // namespace duckdb
