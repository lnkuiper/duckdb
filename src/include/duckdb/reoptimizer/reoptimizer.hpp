//===----------------------------------------------------------------------===//
//                         DuckDB
//
// reoptimizer/reoptimizer.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

// #include <map>

#include "duckdb/main/client_context.hpp"
#include "duckdb/planner/binder.hpp"
#include "duckdb/planner/logical_operator.hpp"

namespace duckdb {
class Binder;

class ReOptimizer {
public:
	ReOptimizer(ClientContext &context, const string &query);

    unique_ptr<LogicalOperator> CreateStepPlan(unique_ptr<LogicalOperator> plan, string temporary_table_name);
    unique_ptr<LogicalOperator> remaining_plan;

    ClientContext &context;

    string remaining_query;

    string step_query;

private:
    void SetTemporaryTableQuery(LogicalComparisonJoin &plan, string table_name);
    void ExtractUsedColumns(const string &query);
    vector<LogicalOperator*> GetJoinOperators(LogicalOperator &plan);
    string JoinStrings(vector<string> strings, string delimiter);
    
    unordered_map<string, vector<string>> columns_per_table = {};
};

} // namespace duckdb
