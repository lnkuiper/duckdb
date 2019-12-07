#include "duckdb/reoptimizer/reoptimizer.hpp"

#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/optimizer/optimizer.hpp"

using namespace std;
using namespace duckdb;

ReOptimizer::ReOptimizer() {
}

unique_ptr<LogicalOperator> ReOptimizer::FirstStepWithTempTable(unique_ptr<LogicalOperator> plan, string temp_table_name) {
    unique_ptr<LogicalOperator> first_step = GetJoinOperators(move(plan)).back;
    // FIXME: change this first step to be a "CREATE TEMPORARY TABLE" plan with temp_table_name
    return plan;
}

vector<unique_ptr<LogicalOperator>> ReOptimizer::GetJoinOperators(unique_ptr<LogicalOperator> plan) {
    // returns a vector of all join operators in the plan
    // the last element is a "join leaf" - has no joins in children
    vector<unique_ptr<LogicalOperator>> join_nodes;
    if (plan->children.empty) {
        return join_nodes;
    }
    switch(plan->GetOperatorType) {
    case LogicalOperatorType::ANY_JOIN:
    case LogicalOperatorType::DELIM_JOIN:
    case LogicalOperatorType::COMPARISON_JOIN:
    case LogicalOperatorType::CROSS_PRODUCT:
        join_nodes.push_back(plan);
    default:
        for (auto &child : plan->children) {
            vector<unique_ptr<LogicalOperator>> children_join_nodes = GetJoinOperators(move(child));
            join_nodes.insert(join_nodes.end(), children_join_nodes.begin(), children_join_nodes.end());
        }
        return join_nodes;
    }
}