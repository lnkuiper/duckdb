#include "catch.hpp"
#include "test_helpers.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/common/set.hpp"
#include "duckdb/common/thread.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/storage/storage_info.hpp"

using namespace duckdb;

namespace {

mutex ordinal_threads_lock;
duckdb::set<thread_id> ordinal_threads;

struct OrdinalProbeLocal : FunctionLocalState {
	bool first_chunk = true;
};

duckdb::unique_ptr<FunctionLocalState> InitOrdinalProbe(ExpressionState &, const BoundFunctionExpression &,
                                                        FunctionData *) {
	return make_uniq<OrdinalProbeLocal>();
}

void ExecuteOrdinalProbe(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &local = ExecuteFunctionState::GetFunctionState(state)->Cast<OrdinalProbeLocal>();
	if (local.first_chunk) {
		local.first_chunk = false;
		{
			duckdb::lock_guard<mutex> guard(ordinal_threads_lock);
			ordinal_threads.insert(ThreadUtil::GetThreadId());
		}
		ThreadUtil::SleepMs(1);
	}
	result.Reference(args.data[0]);
}

} // namespace

TEST_CASE("Source ordinals and random execute on parallel workers", "[input_ordinality]") {
	DuckDB db(nullptr);
	Connection con(db);
	ScalarFunction function("ordinal_probe", {LogicalType::BIGINT}, LogicalType::BIGINT, ExecuteOrdinalProbe, nullptr,
	                        nullptr, InitOrdinalProbe);
	function.SetVolatile();
	CreateScalarFunctionInfo info(function);
	con.context->RunFunctionInTransaction(
	    [&] { Catalog::GetSystemCatalog(*con.context).CreateFunction(*con.context, info); });
	const auto count = idx_t(20) * DEFAULT_ROW_GROUP_SIZE + 37;
	REQUIRE_NO_FAIL(con.Query("PREPARE q AS SELECT sum(ordinal_probe(ordinality))::BIGINT, max(random()) FROM range(" +
	                          to_string(count) + ") WITH ORDINALITY"));
	for (auto threads : {4, 1, 4}) {
		REQUIRE_NO_FAIL(con.Query("SET threads=" + to_string(threads)));
		ordinal_threads.clear();
		auto result = con.Query("EXECUTE q");
		REQUIRE_NO_FAIL(*result);
		auto chunk = result->Fetch();
		REQUIRE(chunk);
		REQUIRE(chunk->GetValue(0, 0).GetValue<int64_t>() == NumericCast<int64_t>(count * (count + 1) / 2));
		REQUIRE(chunk->GetValue(1, 0).GetValue<double>() < 1);
#ifndef DUCKDB_NO_THREADS
		REQUIRE(ordinal_threads.size() == (threads == 1 ? 1 : 4));
#endif
	}
}

TEST_CASE("A schema-local random function keeps its original arguments", "[input_ordinality]") {
	DuckDB db(nullptr);
	Connection con(db);
	con.BeginTransaction();
	REQUIRE_NO_FAIL(con.Query("CREATE SCHEMA custom"));
	ScalarFunction function("random", {}, LogicalType::DOUBLE, [](DataChunk &args, ExpressionState &, Vector &result) {
		REQUIRE(args.ColumnCount() == 0);
		result.Reference(Value::DOUBLE(0.25), count_t(args.size()));
	});
	function.SetVolatile();
	CreateScalarFunctionInfo info(function);
	info.internal = false;
	info.SetQualifiedName(QualifiedName("memory", "custom", "random"));
	Catalog::GetCatalog(*con.context, Identifier("memory")).CreateFunction(*con.context, info);
	auto result = con.Query("SELECT sum(custom.random()) FROM range(4097)");
	REQUIRE_NO_FAIL(*result);
	auto chunk = result->Fetch();
	REQUIRE(chunk);
	REQUIRE(chunk->GetValue(0, 0).GetValue<double>() == 4097 * 0.25);
	con.Rollback();
}
