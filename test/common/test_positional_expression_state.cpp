#include "catch.hpp"
#include "test_helpers.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/common/random_engine.hpp"
#include "duckdb/common/vector/vector_writer.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/function/table/generated_rows.hpp"
#include "duckdb/parser/parsed_data/create_aggregate_function_info.hpp"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/parser/parsed_data/create_table_function_info.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"

using namespace duckdb;

namespace {

mutex preparation_lock;
duckdb::vector<string> prepared_expressions;

struct PositionProbeGlobal : FunctionGlobalState {};

struct PositionProbeLocal : FunctionLocalState {};

duckdb::unique_ptr<FunctionGlobalState> PreparePositionProbe(ClientContext &, const BoundFunctionExpression &expr) {
	duckdb::lock_guard<mutex> guard(preparation_lock);
	prepared_expressions.push_back(expr.Function().GetName().GetIdentifierName());
	return make_uniq<PositionProbeGlobal>();
}

duckdb::unique_ptr<FunctionLocalState> InitializePositionProbe(ExpressionState &, const BoundFunctionExpression &,
                                                               FunctionData *) {
	return make_uniq<PositionProbeLocal>();
}

void PositionProbe(DataChunk &args, ExpressionState &state, Vector &result) {
	int64_t position = -1;
	const bool prepared = bool(state.Cast<ExecuteFunctionState>().global_state);
	if (prepared) {
		auto &metadata = state.root.executor->GetInputMetadata();
		position =
		    NumericCast<int64_t>(metadata.batch.first_row_index.GetIndex() + metadata.offset_in_batch.GetIndex());
	}
	auto writer = FlatVector::Writer<int64_t>(result, args.size());
	for (idx_t row = 0; row < args.size(); row++) {
		writer.WriteValue(prepared ? position + NumericCast<int64_t>(row) : -1);
	}
}

void RegisterPositionProbes(Connection &con) {
	for (auto name : {"position_probe_a", "position_probe_b", "position_probe_parent", "unknown_init_probe"}) {
		duckdb::vector<LogicalType> arguments;
		if (string(name) == "position_probe_parent") {
			arguments.push_back(LogicalType::BIGINT);
		}
		ScalarFunction function(name, arguments, LogicalType::BIGINT, PositionProbe, nullptr, nullptr,
		                        InitializePositionProbe);
		function.SetVolatile();
		function.SetNullHandling(FunctionNullHandling::SPECIAL_HANDLING);
		if (string(name) != "unknown_init_probe") {
			function.SetInitGlobalStateCallback(PreparePositionProbe);
		}
		CreateScalarFunctionInfo info(function);
		con.context->RunFunctionInTransaction(
		    [&] { Catalog::GetSystemCatalog(*con.context).CreateFunction(*con.context, info); });
	}
}

void RegisterPositionSources(Connection &con) {
	for (const bool independent : {false, true}) {
		TableFunction function(
		    independent ? "independent_position_source" : "unknown_position_source", {},
		    [](ClientContext &, TableFunctionInput &input, DataChunk &output) {
			    auto &global = input.global_state->Cast<GeneratedRowsGlobalState>();
			    auto &local = input.local_state->Cast<GeneratedRowsLocalState>();
			    if (!local.Next(global)) {
				    return;
			    }
			    input.metadata = local.GetMetadata(global);
			    output.data[0].Sequence(Hugeint::Cast<int64_t>(local.offset), 1, local.count);
		    },
		    [](ClientContext &, TableFunctionBindInput &, duckdb::vector<LogicalType> &types,
		       duckdb::vector<Identifier> &names) -> duckdb::unique_ptr<FunctionData> {
			    types.push_back(LogicalType::BIGINT);
			    names.emplace_back("i");
			    return nullptr;
		    },
		    [](ClientContext &context, TableFunctionInitInput &input) -> duckdb::unique_ptr<GlobalTableFunctionState> {
			    return make_uniq<GeneratedRowsGlobalState>(context, Hugeint::Convert(4097), input.metadata);
		    },
		    GeneratedRowsLocalState::Init);
		function.get_partition_data = GeneratedRowsLocalState::GetPartitionData;
		function.source_metadata = BatchMetadataInfo::RowPositions();
		if (independent) {
			function.state_initialization = StateInitialization::INDEPENDENT;
		}
		CreateTableFunctionInfo info(function);
		con.context->RunFunctionInTransaction(
		    [&] { Catalog::GetSystemCatalog(*con.context).CreateTableFunction(*con.context, info); });
	}
}

void ClearPositionProbes() {
	duckdb::lock_guard<mutex> guard(preparation_lock);
	prepared_expressions.clear();
}

struct AggregateSeedInfo : AggregateFunctionInfo {
	explicit AggregateSeedInfo(ClientContext &context_p) : context(context_p) {
	}
	ClientContext &context;
	duckdb::vector<uint64_t> seeds;
};

void RegisterAggregateSeedProbe(Connection &con, duckdb::shared_ptr<AggregateSeedInfo> seed_info) {
	AggregateFunction function(
	    "aggregate_seed_probe", {LogicalType::DOUBLE}, LogicalType::BIGINT,
	    [](AggregateStateInput &) -> idx_t { return sizeof(int64_t); },
	    [](AggregateStateInput &input, data_ptr_t *states, idx_t count) {
		    auto &info = input.function.GetExtraFunctionInfo().Cast<AggregateSeedInfo>();
		    auto &engine = RandomEngine::Get(info.context);
		    duckdb::lock_guard<mutex> guard(engine.lock);
		    for (idx_t row = 0; row < count; row++) {
			    info.seeds.push_back(engine.NextRandomInteger64());
			    *reinterpret_cast<int64_t *>(states[row]) = 0;
		    }
	    },
	    [](Vector[], AggregateInputData &, idx_t, Vector &states, idx_t count) {
		    auto entries = states.Values<int64_t *>();
		    for (idx_t row = 0; row < count; row++) {
			    ++*entries[row].GetValue();
		    }
	    },
	    [](Vector &source, Vector &target, AggregateInputData &, idx_t count) {
		    auto sources = source.Values<int64_t *>();
		    auto targets = target.Values<int64_t *>();
		    for (idx_t row = 0; row < count; row++) {
			    *targets[row].GetValue() += *sources[row].GetValue();
		    }
	    },
	    [](Vector &states, AggregateFinalizeInputData &, Vector &result, idx_t count, idx_t offset) {
		    auto entries = states.Values<int64_t *>();
		    auto writer = FlatVector::Writer<int64_t>(result, count, offset);
		    for (idx_t row = 0; row < count; row++) {
			    writer.WriteValue(*entries[row].GetValue());
		    }
	    },
	    FunctionNullHandling::DEFAULT_NULL_HANDLING);
	function.SetExtraFunctionInfo(std::move(seed_info));
	CreateAggregateFunctionInfo info(function);
	con.context->RunFunctionInTransaction(
	    [&] { Catalog::GetSystemCatalog(*con.context).CreateFunction(*con.context, info); });
}

} // namespace

TEST_CASE("Shared expression preparation follows projection order and query generations",
          "[streaming_batch_metadata]") {
	DuckDB db;
	Connection con(db);
	RegisterPositionProbes(con);
	REQUIRE_NO_FAIL(con.Query("PRAGMA verify_parallelism"));
	REQUIRE_NO_FAIL(con.Query("PREPARE positions AS SELECT sum(a), sum(b) FROM "
	                          "(SELECT i, a, position_probe_b() b FROM "
	                          "(SELECT i, position_probe_a() a FROM range(4097) t(i)))"));
	duckdb::vector<idx_t> thread_counts {1};
#ifndef DUCKDB_NO_THREADS
	thread_counts.push_back(4);
#endif
	for (auto threads : thread_counts) {
		REQUIRE_NO_FAIL(con.Query("SET threads=" + std::to_string(threads)));
		for (idx_t execution = 0; execution < 2; execution++) {
			CAPTURE(threads, execution);
			ClearPositionProbes();
			auto result = con.Query("EXECUTE positions");
			REQUIRE_NO_FAIL(*result);
			REQUIRE(result->Collection().GetValue(0, 0).GetValue<int64_t>() == 4097 * 2048);
			REQUIRE(result->Collection().GetValue(1, 0).GetValue<int64_t>() == 4097 * 2048);
			const duckdb::vector<string> expected {"position_probe_a", "position_probe_b"};
			REQUIRE(prepared_expressions == expected);
		}
	}
	ClearPositionProbes();
	auto nested = con.Query("SELECT sum(position_probe_parent(position_probe_a())) FROM range(4097)");
	REQUIRE_NO_FAIL(*nested);
	REQUIRE(nested->Collection().GetValue(0, 0).GetValue<int64_t>() == 4097 * 2048);
	const duckdb::vector<string> expected {"position_probe_a", "position_probe_parent"};
	REQUIRE(prepared_expressions == expected);
}

TEST_CASE("Unsupported expression and stream domains never partially prepare shared seeds",
          "[streaming_batch_metadata]") {
	DuckDB db;
	Connection con(db);
	RegisterPositionProbes(con);
	REQUIRE_NO_FAIL(con.Query("PRAGMA verify_parallelism"));
#ifndef DUCKDB_NO_THREADS
	REQUIRE_NO_FAIL(con.Query("SET threads=4"));
#endif
	const duckdb::vector<string> queries {
	    "SELECT sum(position_probe_a()) FROM range(4097) t(i) WHERE i % 2 = 0",
	    "SELECT sum(a) + sum(b) FROM (SELECT position_probe_a() a, "
	    "CASE WHEN i % 2 = 0 THEN position_probe_b() ELSE 0 END b FROM range(4097) t(i))",
	    "SELECT sum(a), count(u) FROM (SELECT position_probe_a() a, uuid() u FROM range(4097))",
	    "SELECT sum(a) + sum(b) FROM (SELECT position_probe_a() a, unknown_init_probe() b FROM range(4097))",
	    "SELECT sum(position_probe_a()) FROM (SELECT i FROM range(2048) t(i) UNION ALL SELECT i FROM range(2049) t(i))",
	    "SELECT sum(position_probe_a()) FROM range(2048) UNION ALL SELECT sum(position_probe_b()) FROM range(2049)",
	    "SELECT sum(position_probe_a()) FROM "
	    "(SELECT i, unknown_init_probe() b FROM range(4097) t(i)) WHERE b < 0",
	    "SELECT sum((CASE WHEN i < 2048 THEN NULL ELSE 1 END) + position_probe_a()) FROM range(4097) t(i)",
	    "WITH RECURSIVE t(i) AS (SELECT 0 UNION ALL SELECT i+1 FROM t WHERE i<2) "
	    "SELECT sum(position_probe_a()) FROM t"};
	for (auto &query : queries) {
		CAPTURE(query);
		ClearPositionProbes();
		auto result = con.Query(query);
		REQUIRE_NO_FAIL(*result);
		REQUIRE(prepared_expressions.empty());
	}
}

TEST_CASE("Shared expressions require declared source and sink initialization effects", "[streaming_batch_metadata]") {
	DuckDB db;
	Connection con(db);
	RegisterPositionProbes(con);
	RegisterPositionSources(con);
	REQUIRE_NO_FAIL(con.Query("PRAGMA verify_parallelism"));
	for (const bool independent : {false, true}) {
		ClearPositionProbes();
		auto source = independent ? "independent_position_source" : "unknown_position_source";
		auto result = con.Query("SELECT sum(position_probe_a()) FROM " + string(source) + "()");
		REQUIRE_NO_FAIL(*result);
		REQUIRE(result->Collection().GetValue(0, 0).GetValue<int64_t>() == (independent ? 4097 * 2048 : -4097));
		REQUIRE(prepared_expressions.size() == idx_t(independent ? 1 : 0));
	}
	REQUIRE_NO_FAIL(con.Query("CREATE TABLE constrained(i BIGINT CHECK (unknown_init_probe() < 0))"));
	ClearPositionProbes();
	REQUIRE_NO_FAIL(con.Query("INSERT INTO constrained SELECT position_probe_a() FROM range(4097)"));
	REQUIRE(prepared_expressions.empty());
	REQUIRE_NO_FAIL(con.Query("CREATE TABLE indexed(i BIGINT); CREATE INDEX idx ON indexed(i)"));
	ClearPositionProbes();
	REQUIRE_NO_FAIL(con.Query("INSERT INTO indexed SELECT position_probe_a() FROM range(4097)"));
	REQUIRE(prepared_expressions.empty());
}

TEST_CASE("Aggregate extra info retains client seed initialization order", "[streaming_batch_metadata]") {
	DuckDB db;
	Connection con(db);
	RegisterPositionProbes(con);
	auto info = make_shared_ptr<AggregateSeedInfo>(*con.context);
	RegisterAggregateSeedProbe(con, info);
	REQUIRE_NO_FAIL(con.Query("SET threads=1"));
	REQUIRE_NO_FAIL(con.Query("SELECT setseed(0.42)"));
	ClearPositionProbes();
	auto result = con.Query("SELECT sum(a), aggregate_seed_probe(a), sum(p) FROM "
	                        "(SELECT random() a, position_probe_a() p FROM range(4097))");
	REQUIRE_NO_FAIL(*result);
	REQUIRE(result->Collection().GetValue(1, 0).GetValue<int64_t>() == 4097);
	REQUIRE(result->Collection().GetValue(2, 0).GetValue<int64_t>() == -4097);
	REQUIRE(prepared_expressions.empty());
	REQUIRE(!info->seeds.empty());
	// The aggregate global state consumes the first client seed, before projection initialization.
	REQUIRE(info->seeds[0] == UINT64_C(17148587759093342360));
}
