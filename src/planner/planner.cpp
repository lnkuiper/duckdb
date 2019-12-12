#include "duckdb/planner/planner.hpp"

#include "duckdb/common/serializer.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parser/statement/list.hpp"
#include "duckdb/planner/binder.hpp"
#include "duckdb/planner/bound_sql_statement.hpp"
#include "duckdb/planner/expression/bound_parameter_expression.hpp"
#include "duckdb/planner/logical_plan_generator.hpp"
#include "duckdb/planner/operator/logical_explain.hpp"
#include "duckdb/planner/operator/logical_prepare.hpp"
#include "duckdb/planner/statement/bound_select_statement.hpp"
#include "duckdb/planner/query_node/bound_select_node.hpp"
#include "duckdb/planner/query_node/bound_set_operation_node.hpp"

#include "duckdb/planner/pragma_handler.hpp"

using namespace duckdb;
using namespace std;

Planner::Planner(ClientContext &context) : binder(context), context(context) {
}

void Planner::CreatePlan(SQLStatement &statement) {
	vector<BoundParameterExpression *> bound_parameters;

	// first bind the tables and columns to the catalog
	context.profiler.StartPhase("binder");
	binder.parameters = &bound_parameters;
	auto bound_statement = binder.Bind(statement);
	context.profiler.EndPhase();

	VerifyQuery(*bound_statement);

	this->names = bound_statement->GetNames();
	this->sql_types = bound_statement->GetTypes();

	// now create a logical query plan from the query
	context.profiler.StartPhase("logical_planner");
	LogicalPlanGenerator logical_planner(binder, context);
	this->plan = logical_planner.CreatePlan(*bound_statement);
	context.profiler.EndPhase();

	// set up a map of parameter number -> value entries
	for (auto &expr : bound_parameters) {
		// check if the type of the parameter could be resolved
		if (expr->return_type == TypeId::INVALID) {
			throw BinderException("Could not determine type of parameters: try adding explicit type casts");
		}
		auto value = make_unique<Value>(expr->return_type);
		expr->value = value.get();
		// check if the parameter number has been used before
		if (value_map.find(expr->parameter_nr) != value_map.end()) {
			throw BinderException("Duplicate parameter index. Use $1, $2 etc. to differentiate.");
		}
		PreparedValueEntry entry;
		entry.value = move(value);
		entry.target_type = expr->sql_type;
		value_map[expr->parameter_nr] = move(entry);
	}
}

void Planner::CreatePlan(unique_ptr<SQLStatement> statement) {
	assert(statement);
	switch (statement->type) {
	case StatementType::SELECT:
	case StatementType::INSERT:
	case StatementType::COPY:
	case StatementType::DELETE:
	case StatementType::UPDATE:
	case StatementType::CREATE_INDEX:
	case StatementType::CREATE_TABLE:
	case StatementType::EXECUTE:
	case StatementType::CREATE_VIEW:
	case StatementType::CREATE_SCHEMA:
	case StatementType::CREATE_SEQUENCE:
	case StatementType::DROP:
	case StatementType::ALTER:
	case StatementType::TRANSACTION:
	case StatementType::EXPLAIN:
		CreatePlan(*statement);
		break;
	case StatementType::PRAGMA: {
		auto &stmt = *reinterpret_cast<PragmaStatement *>(statement.get());
		PragmaHandler handler(context);
		// some pragma statements have a "replacement" SQL statement that will be executed instead
		// use the PragmaHandler to get the (potential) replacement SQL statement
		auto new_stmt = handler.HandlePragma(*stmt.info);
		if (new_stmt) {
			CreatePlan(move(new_stmt));
		} else {
			CreatePlan(stmt);
		}
		break;
	}
	case StatementType::PREPARE: {
		auto &stmt = *reinterpret_cast<PrepareStatement *>(statement.get());
		auto statement_type = stmt.statement->type;
		// create a plan of the underlying statement
		CreatePlan(move(stmt.statement));
		// now create the logical prepare
		auto prepared_data = make_unique<PreparedStatementData>(statement_type);
		prepared_data->names = names;
		prepared_data->sql_types = sql_types;
		prepared_data->value_map = move(value_map);

		auto prepare = make_unique<LogicalPrepare>(stmt.name, move(prepared_data), move(plan));
		names = {"Success"};
		sql_types = {SQLType(SQLTypeId::BOOLEAN)};
		plan = move(prepare);
		break;
	}
	default:
		throw NotImplementedException("Cannot plan statement of type %s!", StatementTypeToString(statement->type).c_str());
	}
}

void Planner::VerifyQuery(BoundSQLStatement &statement) {
	if (!context.query_verification_enabled) {
		return;
	}
	if (statement.type != StatementType::SELECT) {
		return;
	}
	auto &select = (BoundSelectStatement &)statement;
	VerifyNode(*select.node);
}

void Planner::VerifyNode(BoundQueryNode &node) {
	if (node.type == QueryNodeType::SELECT_NODE) {
		auto &select_node = (BoundSelectNode &)node;
		vector<unique_ptr<Expression>> copies;
		for (auto &expr : select_node.select_list) {
			VerifyExpression(*expr, copies);
		}
		if (select_node.where_clause) {
			VerifyExpression(*select_node.where_clause, copies);
		}
		for (auto &expr : select_node.groups) {
			VerifyExpression(*expr, copies);
		}
		if (select_node.having) {
			VerifyExpression(*select_node.having, copies);
		}
		for (auto &aggr : select_node.aggregates) {
			VerifyExpression(*aggr, copies);
		}
		for (auto &window : select_node.windows) {
			VerifyExpression(*window, copies);
		}

		// double loop to verify that (in)equality of hashes
		for (index_t i = 0; i < copies.size(); i++) {
			auto outer_hash = copies[i]->Hash();
			for (index_t j = 0; j < copies.size(); j++) {
				auto inner_hash = copies[j]->Hash();
				if (outer_hash != inner_hash) {
					// if hashes are not equivalent the expressions should not be equivalent
					assert(!Expression::Equals(copies[i].get(), copies[j].get()));
				}
			}
		}
	} else {
		assert(node.type == QueryNodeType::SET_OPERATION_NODE);
		auto &setop_node = (BoundSetOperationNode &)node;
		VerifyNode(*setop_node.left);
		VerifyNode(*setop_node.right);
	}
}

void Planner::VerifyExpression(Expression &expr, vector<unique_ptr<Expression>> &copies) {
	if (expr.HasSubquery()) {
		// can't copy subqueries
		return;
	}
	// verify that the copy of expressions works
	auto copy = expr.Copy();
	// copy should have identical hash and identical equality function
	assert(copy->Hash() == expr.Hash());
	assert(Expression::Equals(copy.get(), &expr));
	copies.push_back(move(copy));
}
