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
	ReOptimizer();

    unique_ptr<LogicalOperator> ReOptimizer::FirstStepWithTempTable(unique_ptr<LogicalOperator> plan, string temp_table_name);

private:
    vector<unique_ptr<LogicalOperator>> ReOptimizer::GetJoinOperators(unique_ptr<LogicalOperator> plan);
};

} // namespace duckdb
