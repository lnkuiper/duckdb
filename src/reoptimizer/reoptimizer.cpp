#include "duckdb/reoptimizer/reoptimizer.hpp"

#include <regex>
#include <string>

#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/unordered_map.hpp"
#include "duckdb/common/enums/expression_type.hpp"
#include "duckdb/common/enums/logical_operator_type.hpp"
#include "duckdb/main/table_description.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/parser/expression/columnref_expression.hpp"
#include "duckdb/parser/statement/prepare_statement.hpp"
#include "duckdb/planner/joinside.hpp"
#include "duckdb/planner/planner.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"

using namespace std;
using namespace duckdb;

ReOptimizer::ReOptimizer(ClientContext &context, const string &query, LogicalOperator &plan) : context(context) {   
    remaining_query = query;
}

unique_ptr<LogicalOperator> ReOptimizer::CreateStepPlan(unique_ptr<LogicalOperator> plan, string temporary_table_name) {
    vector<LogicalOperator*> join_nodes = ExtractJoinOperators(*plan);
    if (join_nodes.empty()) {
        step_query = "";
        return plan;
    }

    used_columns_per_table = {};
    ExtractUsedColumns(*plan);

    context.profiler.StartPhase("create_step");
    LogicalComparisonJoin* first_join = static_cast<LogicalComparisonJoin*>(join_nodes.back());
    SetTemporaryTableQuery(*first_join, temporary_table_name);
    context.profiler.EndPhase();

    /* TODO: edit 'plan', inserting the temporary table for the join node
     * this will probably be quite hard, because of the binder
     * might need to Skype with Mark/Hannes */

    return plan;
}

void ReOptimizer::SetTemporaryTableQuery(LogicalComparisonJoin &plan, string temporary_table_name) {
    // taken from the overridden LogicalOperator::Print()/ToString() methods
    vector<string> queried_tables;
    vector<string> queried_columns;
    vector<string> where_conditions;
    assert(plan.children.size() == 2);
    for (int i = 0; i < plan.children.size(); i++) {
        auto &child = plan.children.at(i);
        string table_alias = "t" + to_string(i);
        LogicalGet* logical_get;
        switch(child->type) {
        case LogicalOperatorType::GET: {
            logical_get = static_cast<LogicalGet*>(child.get());
            break;
        }
        case LogicalOperatorType::FILTER: {
            assert(child->children.size() == 1);
            LogicalFilter* logical_filter = static_cast<LogicalFilter*>(child.get());
            for (auto &expr : logical_filter->expressions) {
                where_conditions.push_back(table_alias + "." + expr->GetName());
            }
            logical_get = static_cast<LogicalGet*>(logical_filter->children.at(0).get());
            break;
        }
        default:
            throw new ReOptimizerException("Expected child of join to be GET or FILTER, got '%s' instead", child->type);
        }
        /* TODO: adding all the columns in each temporary table will cause unnecessary IO 
         * should traverse the plan to find only the columns that are needed for future joins or the result.
         * Also, the current approach is returning an empty string vector */
        queried_tables.push_back(logical_get->table->schema->name + "." + logical_get->table->name + " AS " + table_alias);
        for (string col : used_columns_per_table[logical_get->table->name]) {
            queried_columns.push_back(table_alias + "." + col);
        }
    }
    for (auto &cond : plan.conditions) {
        where_conditions.push_back("t0." + cond.left->GetName() + " = t1." + cond.right->GetName());
    }
    // create table within a transaction that is never committed for now (error with ActiveTransaction().is_invalidated))
    // step_query = "CREATE TEMPORARY TABLE " + temporary_table_name + " AS ("
    step_query = "CREATE TABLE " + temporary_table_name + " AS ("
               + "SELECT " + JoinStrings(queried_columns, ", ") + " "
               + "FROM " + JoinStrings(queried_tables, ", ") + " "
               + "WHERE " + JoinStrings(where_conditions, " AND ") + ");";
}

/* sets mapping of tablename -> set of columns used in the plan */
void ReOptimizer::ExtractUsedColumns(LogicalOperator &plan) {
    if (plan.children.empty())
        return;
    std::set<string> used_columns;
    switch (plan.type) {
    case LogicalOperatorType::PROJECTION: {
        LogicalProjection* lp = static_cast<LogicalProjection*>(&plan);
        for (int i = 0; i < lp->expressions.size(); i++)
            used_columns.insert(lp->expressions.at(i)->GetName());
        break;
    }
    case LogicalOperatorType::COMPARISON_JOIN: {
        LogicalComparisonJoin* lcj = static_cast<LogicalComparisonJoin*>(&plan);
        for (auto &cond : lcj->conditions) {
            used_columns.insert(cond.left->GetName());
            used_columns.insert(cond.right->GetName());
        }
        break;
    }
    }

    if (!used_columns.empty()) {
        std::set<TableCatalogEntry*> tables = ExtractGetTables(plan);
        for (std::set<TableCatalogEntry*>::iterator it = tables.begin(); it != tables.end(); it++) {
            TableCatalogEntry* table = *it;
            for (ColumnDefinition &cd : table->columns) {
                if (used_columns.find(cd.name) != used_columns.end()) {
                    used_columns_per_table[table->name].insert(cd.name);
                }
            }
        }
    }

    for (auto &child : plan.children) {
        ExtractUsedColumns(*child);
    }
}

/* returns a vector of all join operators in the plan
 * the last element should have two children, both of which have LogicalGet as children */
vector<LogicalOperator*> ReOptimizer::ExtractJoinOperators(LogicalOperator &plan) {
    vector<LogicalOperator*> joins;
    if (plan.children.empty())
        return joins;
    switch(plan.type) {
    // case LogicalOperatorType::JOIN:
    // case LogicalOperatorType::ANY_JOIN:
    // case LogicalOperatorType::DELIM_JOIN:
    // case LogicalOperatorType::CROSS_PRODUCT:
    // Pretty sure we only get COMPARISON_JOIN in JOB
    case LogicalOperatorType::COMPARISON_JOIN:
        joins.push_back(&plan);
    default:
        for (auto &child : plan.children) {
            vector<LogicalOperator*> child_joins = ExtractJoinOperators(*child);
            joins.insert(joins.end(), child_joins.begin(), child_joins.end());
        }
        return joins;
    }
}

std::set<TableCatalogEntry*> ReOptimizer::ExtractGetTables(LogicalOperator &plan) {
    std::set<TableCatalogEntry*> gets;
    if (plan.type == LogicalOperatorType::GET) {
        LogicalGet* lg = static_cast<LogicalGet*>(&plan);
        gets.insert(lg->table);
    }
    for (auto &child : plan.children) {
        std::set<TableCatalogEntry*> child_gets = ExtractGetTables(*child);
        for (std::set<TableCatalogEntry*>::iterator it = child_gets.begin(); it != child_gets.end(); it++)
            gets.insert(*it);
    }
    return gets;
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
