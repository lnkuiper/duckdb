#include "duckdb/main/connection.hpp"

#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection_manager.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/appender.hpp"
#include "duckdb/parser/parser.hpp"

using namespace duckdb;
using namespace std;

Connection::Connection(DuckDB &database) : db(database), context(make_unique<ClientContext>(database)) {
	db.connection_manager->AddConnection(this);
#ifdef DEBUG
	EnableProfiling();
#endif
}

Connection::~Connection() {
	if (!context->is_invalidated) {
		context->Cleanup();
		db.connection_manager->RemoveConnection(this);
	}
}

string Connection::GetProfilingInformation(ProfilerPrintFormat format) {
	if (context->is_invalidated) {
		return "Context is invalidated.";
	}
	if (format == ProfilerPrintFormat::JSON) {
		return context->profiler.ToJSON();
	} else {
		return context->profiler.ToString();
	}
}

void Connection::Interrupt() {
	context->Interrupt();
}

void Connection::EnableProfiling() {
	context->EnableProfiling();
}

void Connection::DisableProfiling() {
	context->DisableProfiling();
}

void Connection::EnableQueryVerification() {
#ifdef DEBUG
	context->query_verification_enabled = true;
#endif
}

unique_ptr<QueryResult> Connection::SendQuery(string query) {
	return context->Query(query, true);
}

unique_ptr<MaterializedQueryResult> Connection::Query(string query) {
	auto result = context->Query(query, false);
	assert(result->type == QueryResultType::MATERIALIZED_RESULT);
	return unique_ptr_cast<QueryResult, MaterializedQueryResult>(move(result));
}

unique_ptr<PreparedStatement> Connection::Prepare(string query) {
	return context->Prepare(query);
}

unique_ptr<QueryResult> Connection::QueryParamsRecursive(string query, vector<Value> &values) {
	auto statement = Prepare(query);
	if (!statement->success) {
		return make_unique<MaterializedQueryResult>(statement->error);
	}
	return statement->Execute(values);
}

unique_ptr<TableDescription> Connection::TableInfo(string table_name) {
	return TableInfo(DEFAULT_SCHEMA, table_name);
}

unique_ptr<TableDescription> Connection::TableInfo(string schema_name, string table_name) {
	return context->TableInfo(schema_name, table_name);
}

vector<unique_ptr<SQLStatement>> Connection::ExtractStatements(string query) {
	Parser parser;
	parser.ParseQuery(query);
	return move(parser.statements);
}

void Connection::Append(TableDescription &description, DataChunk &chunk) {
	context->Append(description, chunk);
}
