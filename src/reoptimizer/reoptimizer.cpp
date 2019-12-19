#include "duckdb/reoptimizer/reoptimizer.hpp"

#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/enums/expression_type.hpp"
#include "duckdb/common/enums/logical_operator_type.hpp"
#include "duckdb/main/table_description.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/parser/expression/columnref_expression.hpp"
// #include "duckdb/parser/query_node.hpp"
#include "duckdb/planner/joinside.hpp"
#include "duckdb/planner/planner.hpp"
// #include "duckdb/planner/operator/logical_any_join.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
// #include "duckdb/planner/operator/logical_cross_product.hpp"
// #include "duckdb/planner/operator/logical_delim_join.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"

using namespace std;
using namespace duckdb;

ReOptimizer::ReOptimizer(ClientContext &context, Binder &binder) : context(context), binder(binder) {
}

unique_ptr<LogicalOperator> ReOptimizer::CreateFirstStepPlan(unique_ptr<LogicalOperator> plan, string temporary_table_name) {
    vector<LogicalOperator*> join_nodes = GetJoinOperators(*plan);
    if (join_nodes.empty())
        return plan;
    LogicalComparisonJoin* first_join = static_cast<LogicalComparisonJoin*>(join_nodes.back());
    context.profiler.StartPhase("create_step");
    SetTemporaryTableQuery(*first_join, temporary_table_name);

    context.profiler.EndPhase();
    return plan;
}

/* returns a vector of all join operators in the plan
 * the last element should have two children, both of which have LogicalGet as children */
vector<LogicalOperator*> ReOptimizer::GetJoinOperators(LogicalOperator &plan) {
    vector<LogicalOperator*> join_nodes;
    if (plan.children.empty())
        return join_nodes;
    switch(plan.GetOperatorType()) {
    // Pretty sure we only get COMPARISON_JOIN in JOB
    // case LogicalOperatorType::JOIN:
    // case LogicalOperatorType::ANY_JOIN:
    // case LogicalOperatorType::DELIM_JOIN:
    // case LogicalOperatorType::CROSS_PRODUCT:
    case LogicalOperatorType::COMPARISON_JOIN:
        join_nodes.push_back(&plan);
    default:
        for (auto &child : plan.children) {
            vector<LogicalOperator*> children_join_nodes = GetJoinOperators(*child);
            join_nodes.insert(join_nodes.end(), children_join_nodes.begin(), children_join_nodes.end());
        }
        return join_nodes;
    }
}

void ReOptimizer::SetTemporaryTableQuery(LogicalComparisonJoin &plan, string temporary_table_name) {
    vector<string> queried_tables;
    vector<string> queried_columns;
    vector<string> where_conditions;
    assert(plan.children.size() == 2);
    for (int i = 0; i < plan.children.size(); i++) {
        auto &child = plan.children.at(i);
        string table_alias = "t" + to_string(i);
        switch(child->GetOperatorType()) {
        case LogicalOperatorType::GET: {
            LogicalGet* logical_get = static_cast<LogicalGet*>(child.get());
            queried_tables.push_back(logical_get->table->schema->name + "." + logical_get->table->name + " AS " + table_alias);
            for (ColumnDefinition &col : context.TableInfo(logical_get->table->schema->name, logical_get->table->name)->columns) {
                queried_columns.push_back(table_alias + "." + col.name);
            }
            break;
        }
        case LogicalOperatorType::FILTER: {
            assert(child->children.size() == 1);
            LogicalFilter* logical_filter = static_cast<LogicalFilter*>(child.get());
            LogicalGet* logical_get = static_cast<LogicalGet*>(logical_filter->children.at(0).get());
            queried_tables.push_back(logical_get->table->schema->name + "." + logical_get->table->name + " AS " + table_alias);
            for (auto &expr : logical_filter->expressions) {
                // TODO: push_back where_conditions
                
            }
            break;
        }
        default:
            throw new ReOptimizerException("Expected child of join to be GET or FILTER, got '%s' instead", child->GetOperatorType());
        }
    }
    for (auto &condition : plan.conditions) {
        where_conditions.push_back("t1." + condition.left->GetName() + " = t2." + condition.right->GetName());
    }
    step_query = "CREATE TEMPORARY TABLE " + temporary_table_name + " AS ("
               + "SELECT " + JoinStrings(queried_columns, ", ") + " "
               + "FROM " + JoinStrings(queried_tables, ", ") + " "
               + "WHERE " + JoinStrings(where_conditions, " AND ") + ");";
}

string ReOptimizer::JoinStrings(vector<string> strings, string delimiter) {
    string joined_strings = "";
    for (int i = 0; i < strings.size(); i++) {
        joined_strings += strings.at(i);
        if (i < strings.size() - 1)
            joined_strings += delimiter;
    }
    return joined_strings;
}
