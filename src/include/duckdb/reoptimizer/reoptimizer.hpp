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

    unique_ptr<LogicalOperator> ReOptimizer::FirstStepAsTempTable(unique_ptr<LogicalOperator> plan, string table_name);

private:
    vector<unique_ptr<LogicalOperator>> ReOptimizer::GetJoinOperators(unique_ptr<LogicalOperator> plan);

    string ReOptimizer::CreateTemporaryTableQuery(unique_ptr<LogicalOperator> plan, string table_name);
};

} // namespace duckdb
