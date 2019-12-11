#include "duckdb/reoptimizer/reoptimizer.hpp"

#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/planner/operator/logical_any_join.hpp"
#include "duckdb/planner/operator/logical_delim_join.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_cross_product.hpp"

using namespace std;
using namespace duckdb;

ReOptimizer::ReOptimizer() {
}

unique_ptr<LogicalOperator> ReOptimizer::FirstStepAsTempTable(unique_ptr<LogicalOperator> plan, string table_name) {
    vector<unique_ptr<LogicalOperator>> join_nodes = GetJoinOperators(move(plan));
    if (join_nodes.empty)
        return plan;
    unique_ptr<LogicalOperator> first_join = join_nodes.back;
    /* TODO: change this first step to be a "CREATE TEMPORARY TABLE" plan with temp_table_name
     * It needs a projection and create table on top of it
     * Needs SchemaCatalogInformation and BoundCreateTableInfo
     * Check how the plan is created in Planner for examples
     * Planner calls binder... the original query is bound to a (temporary?) table, even for select statements
     * Hopefully not a problem for this approach
     * 
     * Optimizer needs the planner.binder - planner needs "unique_ptr<SQLStatement> statement"
     * We need to go back from LogicalOperator to SQLStatement
     * The only easy way to do this is by going back to "string query" and giving it to the parser
     * Plan: get the table names and columns from first_join and make a new query
     */
    string temp_table_query = CreateTemporaryTableQuery(move(first_join), table_name);
    
    // Needs to return createtable
    return first_join;
}

/* returns a vector of all join operators in the plan
 * the last element has no joins in children */
vector<unique_ptr<LogicalOperator>> ReOptimizer::GetJoinOperators(unique_ptr<LogicalOperator> plan) {
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

string ReOptimizer::CreateTemporaryTableQuery(unique_ptr<LogicalOperator> plan, string table_name) {
    string query = "CREATE TEMPORARY TABLE " + table_name + " AS ";
    // TODO: get children here
    switch(plan->GetOperatorType) {
    case LogicalOperatorType::ANY_JOIN: {
        LogicalAnyJoin* join = (LogicalAnyJoin *) plan.get();
        string condition = join->condition->ToString;
        break;
    }
    case LogicalOperatorType::DELIM_JOIN: {
        LogicalDelimJoin* join = (LogicalDelimJoin *) plan.get();
        vector<unique_ptr<Expression>> columns = join->duplicate_eliminated_columns;
        break;
    }
    case LogicalOperatorType::COMPARISON_JOIN: {
        LogicalComparisonJoin* join = (LogicalComparisonJoin *) plan.get();
        vector<JoinCondition> conditions = join->conditions;
        break;
    }
    case LogicalOperatorType::CROSS_PRODUCT: {
        LogicalCrossProduct* join = (LogicalCrossProduct *) plan.get();
        // easy, no conditions
        break;
    }
    default:
        throw ReOptimizerException("Expected a join type, got '%s' instead", plan->GetOperatorType.ToString);
    }
    return "";
}
