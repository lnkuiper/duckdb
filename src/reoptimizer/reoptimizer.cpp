#include "duckdb/reoptimizer/reoptimizer.hpp"

#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/enums/expression_type.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/table_description.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/parser/expression/columnref_expression.hpp"
#include "duckdb/planner/joinside.hpp"
#include "duckdb/planner/planner.hpp"
#include "duckdb/planner/operator/logical_any_join.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_cross_product.hpp"
#include "duckdb/planner/operator/logical_delim_join.hpp"
#include "duckdb/planner/operator/logical_get.hpp"

using namespace std;
using namespace duckdb;

ReOptimizer::ReOptimizer(ClientContext &context) : context(context) {
}

string ReOptimizer::CreateFirstStepQuery(unique_ptr<LogicalOperator> plan, string temporary_table_name) {
    vector<unique_ptr<LogicalOperator>> join_nodes = GetJoinOperators(move(plan));
    if (join_nodes.empty)
        return "";
    context.profiler.StartPhase("create_step");
    unique_ptr<LogicalOperator> first_join = join_nodes.back;
    string query = CreateTemporaryTableQuery(move(first_join), temporary_table_name);

    // TODO: need to remove this first step from the plan now (replace with GET on temp table)
    // Do this by setting remaining_plan in this method - can take it later
    // does "plan" change if "first join" is changed to a GET on the temp table?
    // it is a unique_ptr after all
    // apparently the object is disposed of when the unique_ptr goes out of scope
    // that would mean that this class does not even need "remaining_plan"
    // it does mean that this method must call move(plan), probably
    // OK it is probably easier to set remaining_plan = move(plan) after all

    // This method needs to return the plan
    // set a string variable variable with the step_plan
    // that way we don't have to fuck around with pointers
    
    // Needs to return createtable
    context.profiler.EndPhase();
    return query;
}

/* returns a vector of all join operators in the plan
 * the last element should have two children, both of which have LogicalGet as children */
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

vector<string> ReOptimizer::ColumnNamesFromLogicalGet(LogicalGet* logical_get) {
    vector<column_t> queried_column_ids = logical_get->column_ids;
    vector<string> column_names;
    for (column_t column_id : queried_column_ids) {
        for (ColumnDefinition &cd : logical_get->table->columns) {
            if (column_id == cd.oid) {
                column_names.push_back(cd.name);
                break;
            }
        }
    }
    return column_names;
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

string ReOptimizer::JoinConditionFromLogicalPlan(unique_ptr<LogicalOperator> plan) {
    switch(plan->GetOperatorType) {
    case LogicalOperatorType::ANY_JOIN: {
        LogicalAnyJoin* join = (LogicalAnyJoin *) plan.get();
        // this join has a unique_ptr<Expression> condition, instead of vector<JoinCondition> conditions (like LogicalComparisonJoin),
        // which makes it much harder to parse (many switches for expressiontype and such)
        break;
    }
    case LogicalOperatorType::DELIM_JOIN: {
        // I don't fully understand  what this join type does
        LogicalDelimJoin* join = (LogicalDelimJoin *) plan.get();
        vector<unique_ptr<Expression>> columns = join->duplicate_eliminated_columns;
        break;
    }
    case LogicalOperatorType::COMPARISON_JOIN: {
        LogicalComparisonJoin* join = (LogicalComparisonJoin *) plan.get();
        vector<string> condition_strings;
        for (JoinCondition &condition : join->conditions) {
            string comparison_str;
            switch (condition.comparison) {
            case ExpressionType::COMPARE_EQUAL:
                comparison_str = " = ";
                break;
            case ExpressionType::COMPARE_NOTEQUAL:
                comparison_str = " != ";
                break;
            default:
                throw new ReOptimizerException("Expected ExpressionType::COMPARE_EQUAL, got 'ExpressionType%s' instead", condition.comparison);
            }
            ColumnRefExpression* l_expr = (ColumnRefExpression *) condition.left.get();
            ColumnRefExpression* r_expr = (ColumnRefExpression *) condition.right.get();
            condition_strings.push_back(l_expr->table_name + "." + l_expr->column_name + comparison_str +
                                        r_expr->table_name + "." + r_expr->column_name);
        }
        return " WHERE " + JoinStrings(condition_strings, " AND ");
    }
    case LogicalOperatorType::CROSS_PRODUCT: {
        // cross product has empty join condition
        return "";
    }
    default:
        throw ReOptimizerException("Expected a join type, got '%s' instead", plan->GetOperatorType.ToString);
    }
}

string ReOptimizer::CreateTemporaryTableQuery(unique_ptr<LogicalOperator> plan, string temporary_table_name) {
    vector<string> queried_columns;
    vector<string> queried_tables;
    for (int i = 0; i < plan->children.size(); i++) {
        auto &child = plan->children.at(i);
        switch(child->GetOperatorType) {
        case LogicalOperatorType::GET: {
            LogicalGet* logical_get = (LogicalGet *) child.get();
            for (string column_name : ColumnNamesFromLogicalGet(logical_get))
                queried_columns.push_back("t" + to_string(i) + "." + column_name);
            queried_tables.push_back(logical_get->table->schema->name + "." + logical_get->table->name);
        }
        case LogicalOperatorType::FILTER: {
            // TODO: FILTER is also possible, I forgot - filtering on multiple can be done in one LogicalFilter
            break;
        }
        default:
            throw new ReOptimizerException("Expected child of join to be GET or FILTER, got '%s' instead", child->GetOperatorType);
        }
    }
    return "CREATE TEMPORARY TABLE " + temporary_table_name + " AS ("
             + "SELECT " + JoinStrings(queried_columns, ", ") + " "
             + "FROM " + JoinStrings(queried_tables, ", ") + " "
             + JoinConditionFromLogicalPlan(move(plan)) + ");";
}
