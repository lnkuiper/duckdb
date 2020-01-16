//===----------------------------------------------------------------------===//
//                         DuckDB
//
// reoptimizer/reoptimizer.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include <set>

#include "duckdb/main/client_context.hpp"
#include "duckdb/planner/binder.hpp"
#include "duckdb/planner/logical_operator.hpp"

namespace duckdb {
class Binder;

class ReOptimizer {
public:
	ReOptimizer(ClientContext &context, const string &query, LogicalOperator &plan);

    unique_ptr<LogicalOperator> CreateStepPlan(unique_ptr<LogicalOperator> plan, string temporary_table_name);
    unique_ptr<LogicalOperator> remaining_plan;

    ClientContext &context;

    string remaining_query;

    string step_query;

private:
    void SetTemporaryTableQuery(LogicalComparisonJoin &plan, string table_name);
    // void ExtractUsedColumns(const string &query);
    void ExtractUsedColumns(LogicalOperator &plan);
    vector<LogicalOperator*> ExtractJoinOperators(LogicalOperator &plan);
    std::set<TableCatalogEntry*> ExtractGetTables(LogicalOperator &plan);
    string JoinStrings(vector<string> strings, string delimiter);
    unordered_map<string, std::set<string>> used_columns_per_table;
};

} // namespace duckdb
