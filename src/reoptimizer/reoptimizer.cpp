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

using namespace std;
using namespace duckdb;

ReOptimizer::ReOptimizer(ClientContext &context, const string &query) : context(context) {   
    remaining_query = query;
    ExtractUsedColumns(query);
}

unique_ptr<LogicalOperator> ReOptimizer::CreateStepPlan(unique_ptr<LogicalOperator> plan, string temporary_table_name) {
    vector<LogicalOperator*> join_nodes = GetJoinOperators(*plan);
    if (join_nodes.empty()) {
        step_query = "";
        return plan;
    }
    
    context.profiler.StartPhase("create_step");
    LogicalComparisonJoin* first_join = static_cast<LogicalComparisonJoin*>(join_nodes.back());
    SetTemporaryTableQuery(*first_join, temporary_table_name);
    context.profiler.EndPhase();

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
        switch(child->GetOperatorType()) {
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
            throw new ReOptimizerException("Expected child of join to be GET or FILTER, got '%s' instead", child->GetOperatorType());
        }
        /* TODO: adding all the columns in each temporary table will cause unnecessary IO 
         * should traverse the plan to find only the columns that are needed for future joins or the result.
         * Also, the current approach is returning an empty string vector */
        queried_tables.push_back(logical_get->table->schema->name + "." + logical_get->table->name + " AS " + table_alias);
        for (string col : columns_per_table[logical_get->table->name]) {
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

/* sets mapping of tablename -> vector of columns of the columns used in the plan */
void ReOptimizer::ExtractUsedColumns(const string &query) {
    // map from aliases to tablename
    unordered_map<string, string> aliases;
    
    /* FIXME: regex implementation for now. Not happy with this,
     * but approaches using SQLStatement and LogicalOperator did not work */

    // Extract FROM block
    smatch from_smatch;
	regex from_re("FROM(.*?)WHERE");
	if (regex_search(query, from_smatch, from_re)) {
		string from_clauses = from_smatch[1];

        // Match aliases in FROM, add to map
		smatch alias_smatch;
		regex alias_re("([^ ]+ AS [^ ,]+)");
		string::const_iterator search_start(from_clauses.cbegin());
		while (regex_search(search_start, from_clauses.cend(), alias_smatch, alias_re)) {
			string alias_def = alias_smatch[0];
            aliases[alias_def.substr(alias_def.find(" AS ") + 4, alias_def.length())] =
                                                alias_def.substr(0, alias_def.find(" AS "));
			search_start = alias_smatch.suffix().first;
		}
	}

    // Map from full table names to queried columns
    for (unordered_map<string, string>::iterator it = aliases.begin(); it != aliases.end(); it++) {
        columns_per_table[it->second] = vector<string>();
    }

    // Extract SELECT block
    smatch select_smatch;
    regex select_re("SELECT(.*?)FROM");
    if (regex_search(query, select_smatch, select_re)) {
    	string selects = select_smatch[1];

        // Extract queried columns (don't have to worry about 'SELECT *' in JOB)
    	smatch col_smatch;
    	regex col_re("([^ .\\(]+\\.[^ .\\)]+)");
    	string::const_iterator search_start(selects.cbegin());
		while (regex_search(search_start, selects.cend(), col_smatch, col_re)) {
			string col = col_smatch[0];
			string alias = col.substr(0, col.find("."));
            string column = col.substr(col.find("."), col.length());
            columns_per_table[aliases[alias]].push_back(column);
			search_start = col_smatch.suffix().first;
		}
    }

    // Extract WHERE block
    smatch where_smatch;
	regex where_re("WHERE(.*?);");
	if (regex_search(query, where_smatch, where_re)) {
		string where_clauses = where_smatch[1];

        // Match equality joins
		smatch col_smatch;
		regex col_re("([^ ']+ = [^ ',]+)");
		string::const_iterator search_start(where_clauses.cbegin());
		while (regex_search(search_start, where_clauses.cend(), col_smatch, col_re)) {
			string col_def = col_smatch[0];

            // Extract columns from both sides of equality
            string lhs = col_def.substr(0, col_def.find(" = "));
            string l_alias = lhs.substr(0, lhs.find("."));
            string l_column =  lhs.substr(lhs.find(".") + 1, lhs.length());
            columns_per_table[aliases[l_alias]].push_back(l_column);

            string rhs = col_def.substr(col_def.find(" = ") + 3, col_def.length());
            string r_alias = rhs.substr(0, rhs.find("."));
            string r_column = rhs.substr(rhs.find(".") + 1, rhs.length());
            columns_per_table[aliases[r_alias]].push_back(r_column);

			search_start = col_smatch.suffix().first;
		}
	}
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

string ReOptimizer::JoinStrings(vector<string> strings, string delimiter) {
    string joined_strings = "";
    for (int i = 0; i < strings.size(); i++) {
        joined_strings += strings.at(i);
        if (i < strings.size() - 1)
            joined_strings += delimiter;
    }
    return joined_strings;
}
