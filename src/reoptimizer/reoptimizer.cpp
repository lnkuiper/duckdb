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

unique_ptr<LogicalOperator> ReOptimizer::CreateFirstStepPlan(unique_ptr<LogicalOperator> plan, string temporary_table_name) {
    vector<LogicalOperator> join_nodes = GetJoinOperators(*plan);
    if (join_nodes.empty)
        return plan;
    LogicalOperator &first_join = join_nodes.back;
    context.profiler.StartPhase("create_step");
    SetTemporaryTableQuery(first_join, temporary_table_name);

    context.profiler.EndPhase();
    return plan;
}

/* returns a vector of all join operators in the plan
 * the last element should have two children, both of which have LogicalGet as children */
vector<LogicalOperator> ReOptimizer::GetJoinOperators(LogicalOperator &plan) {
    vector<LogicalOperator> join_nodes;
    if (plan.children.empty)
        return join_nodes;
    switch(plan.GetOperatorType) {
    case LogicalOperatorType::JOIN:
    case LogicalOperatorType::ANY_JOIN:
    case LogicalOperatorType::DELIM_JOIN:
    case LogicalOperatorType::COMPARISON_JOIN:
    case LogicalOperatorType::CROSS_PRODUCT:
        join_nodes.push_back(plan);
        break;
    default:
        for (auto &child : plan.children) {
            vector<LogicalOperator> children_join_nodes = GetJoinOperators(*child);
            join_nodes.insert(join_nodes.end(), children_join_nodes.begin(), children_join_nodes.end());
        }
        return join_nodes;
    }
}

void ReOptimizer::SetTemporaryTableQuery(LogicalOperator &plan, string temporary_table_name) {
    vector<string> queried_columns;
    vector<string> queried_tables;
    string schema_name;
    for (auto &child : plan.children) {
        switch(child->GetOperatorType) {
        case LogicalOperatorType::GET: {
            LogicalGet* logical_get = static_cast<LogicalGet*>(child.get());
            schema_name = logical_get->table->schema->name;
            string table_name = logical_get->table->name;
            for (string column_name : ColumnNamesFromLogicalGet(*logical_get))
                queried_columns.push_back(schema_name + "." + table_name + "." + column_name);
            queried_tables.push_back(schema_name + "." + table_name);
            break;
        }
        case LogicalOperatorType::FILTER: {
            // TODO: FILTER is also possible, I forgot - filtering on multiple can be done in one LogicalFilter
            break;
        }
        default:
            throw new ReOptimizerException("Expected child of join to be GET or FILTER, got '%s' instead", child->GetOperatorType);
        }
    }
    step_query = "CREATE TEMPORARY TABLE " + temporary_table_name + " AS ("
               + "SELECT " + JoinStrings(queried_columns, ", ") + " "
               + "FROM " + JoinStrings(queried_tables, ", ") + " "
               + JoinConditionFromLogicalPlan(plan, schema_name) + ");";
}

vector<string> ReOptimizer::ColumnNamesFromLogicalGet(LogicalGet &logical_get) {
    // TODO: this method will not do what it is supposed to
    // The column_ids found in the logical_get are only used within the context of the binder
    // We need the binder
    vector<column_t> queried_column_ids = logical_get.column_ids;
    vector<string> column_names;
    for (column_t column_id : queried_column_ids) {
        for (ColumnDefinition &cd : logical_get.table->columns) {
            if (column_id == cd.oid) {
                column_names.push_back(cd.name);
                break;
            }
        }
    }
    return column_names;
}

string ReOptimizer::JoinConditionFromLogicalPlan(LogicalOperator &plan, string schema) {
    switch(plan.GetOperatorType) {
    case LogicalOperatorType::ANY_JOIN: {
        LogicalAnyJoin &join = static_cast<LogicalAnyJoin&>(plan);
        // this join has a unique_ptr<Expression> condition, instead of vector<JoinCondition> conditions (like LogicalComparisonJoin),
        // which makes it much harder to parse (many switches for expressiontype and such)
        break;
    }
    case LogicalOperatorType::DELIM_JOIN: {
        // I don't fully understand  what this join type does
        LogicalDelimJoin &join = static_cast<LogicalDelimJoin&>(plan);
        vector<unique_ptr<Expression>> columns = join.duplicate_eliminated_columns;
        break;
    }
    case LogicalOperatorType::COMPARISON_JOIN: {
        LogicalComparisonJoin &join = static_cast<LogicalComparisonJoin&>(plan);
        vector<string> condition_strings;
        for (JoinCondition &condition : join.conditions) {
            string comparison_str;
            switch (condition.comparison) {
            case ExpressionType::COMPARE_EQUAL:
                comparison_str = " = ";
                break;
            case ExpressionType::COMPARE_NOTEQUAL:
                comparison_str = " != ";
                break;
            default:
                throw new ReOptimizerException("Expected ExpressionType::COMPARE_EQUAL or ExpressionType::COMPARE_NOTEQUAL, got 'ExpressionType%s' instead", condition.comparison);
            }
            // Pretty sure that the expression found in bottom-level joins have to be bound_columnref_expressions
            BoundColumnRefExpression* l_expr = static_cast<BoundColumnRefExpression*>(condition.left.get());            
            BoundColumnRefExpression* r_expr = static_cast<BoundColumnRefExpression*>(condition.left.get());

            // Documentation says: "Column index set by the binder, used to generate the final BoundExpression"
            // That means that these indices are only used internally within this plan
            // Need to find a way to go back to the names
            // We probably need the binder for this, however the binder only has a map<string, index_t>
            // We cannot simply make a map<index_t, string> here, since the same index_t will likely be re-used
            // We need the binder to have a map<pair<index_t,index_t>, string> so we can go from (table_idx, col_idx) to "table.column"
            // Keep looking in bind_insert.cpp to see if this is possible
            index_t l_table_idx = l_expr->binding.table_index;
            index_t l_col_idx = l_expr->binding.column_index;
            
            // Given a schema name, we can get a TableDescription from 
            // TableDescription td = context.TableInfo(schema, )

            // ColumnRefExpression* l_expr = static_cast<ColumnRefExpression*>(condition.left.get());
            // ColumnRefExpression* r_expr = static_cast<ColumnRefExpression*>(condition.right.get());
            // condition_strings.push_back(l_expr->table_name + "." + l_expr->column_name + comparison_str +
            //                             r_expr->table_name + "." + r_expr->column_name);
        }
        return " WHERE " + JoinStrings(condition_strings, " AND ");
    }
    case LogicalOperatorType::CROSS_PRODUCT: {
        // cross product has empty join condition
        return "";
    }
    default:
        throw ReOptimizerException("Expected a join type, got '%s' instead", plan.GetOperatorType.ToString);
    }
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
