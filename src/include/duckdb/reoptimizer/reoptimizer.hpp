//===----------------------------------------------------------------------===//
//                         DuckDB
//
// reoptimizer/reoptimizer.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/main/client_context.hpp"
#include "duckdb/planner/binder.hpp"
#include "duckdb/planner/logical_operator.hpp"

namespace duckdb {
class Binder;

class ReOptimizer {
public:
	ReOptimizer(Binder &binder, ClientContext &context);

    unique_ptr<LogicalOperator> CreateFirstStepPlan(unique_ptr<LogicalOperator> plan, string table_name);
    unique_ptr<LogicalOperator> remaining_plan;

    ClientContext &context;
    Binder &binder;

    string step_query;

private:
    void SetTemporaryTableQuery(LogicalComparisonJoin &plan, string table_name);
    vector<LogicalOperator*> GetJoinOperators(LogicalOperator &plan);
    string JoinStrings(vector<string> strings, string delimiter);
};

} // namespace duckdb
