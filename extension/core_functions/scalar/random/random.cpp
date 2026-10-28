#include "core_functions/scalar/random_functions.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/client_context_state.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/common/random_engine.hpp"
#include "duckdb/common/types/uuid.hpp"
#include "duckdb/common/types/timestamp.hpp"

namespace duckdb {

namespace {

struct ExtractVersionUuidOperator {
	template <typename INPUT_TYPE, typename RESULT_TYPE>
	static RESULT_TYPE Operation(INPUT_TYPE input) {
		char uuid[36]; // Intentionally no initialize.
		BaseUUID::ToString(input, uuid);
		// UUIDv4 and UUIDv7 stores version as the 15-th uint8_t.
		return uuid[14] - '0';
	}
};

struct ExtractTimestampUuidOperator {
	template <typename INPUT_TYPE, typename RESULT_TYPE>
	static RESULT_TYPE Operation(INPUT_TYPE input) {
		// Validate whether the given UUID is v7.
		const uint8_t version = (static_cast<uint8_t>((input.upper) >> 8) & 0xf0) >> 4;
		if (version != 7) {
			throw InvalidInputException("Given UUID is with version %u, not version 7.", version);
		}

		// UUID v7 begins with a 48 bit big-endian Unix Epoch timestamp with millisecond granularity.
		int64_t upper = input.upper;
		// flip the top byte
		upper ^= NumericLimits<int64_t>::Minimum();
		int64_t unix_ts_milli = upper;
		unix_ts_milli = unix_ts_milli >> 16;

		static constexpr int64_t kMilliToMicro = 1000;
		const int64_t unix_ts_ms = kMilliToMicro * unix_ts_milli;
		return timestamp_t {unix_ts_ms};
	}
};

template <typename INPUT, typename OP>
void ExtractVersionFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	D_ASSERT(args.ColumnCount() == 1);
	const auto &input = args.data[0];
	UnaryExecutor::Execute<INPUT, uint32_t, OP>(input, result);
}

template <typename INPUT, typename OP>
void ExtractTimestampFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	D_ASSERT(args.ColumnCount() == 1);
	const auto &input = args.data[0];
	UnaryExecutor::Execute<INPUT, timestamp_t, OP>(input, result);
}

struct RandomQueryState : public ClientContextState {
	static uint64_t GetSeed(ClientContext &context, idx_t query_id, bool per_worker = false) {
		auto state = context.registered_state->GetOrCreate<RandomQueryState>("random_query_seed");
		lock_guard<mutex> guard(state->lock);
		if (!state->query_id.IsValid() || state->query_id.GetIndex() != query_id) {
			auto &random_engine = RandomEngine::Get(context);
			lock_guard<mutex> random_guard(random_engine.lock);
			state->seed = random_engine.NextRandomInteger64();
			// Worker allocation must not advance the session RNG or reuse the ordinal streams.
			state->worker_seeds.SetSeed(~state->seed);
			state->query_id = query_id;
		}
		return per_worker ? state->worker_seeds.NextRandomInteger64() : state->seed;
	}

	mutex lock;
	optional_idx query_id;
	uint64_t seed = 0;
	RandomEngine worker_seeds {0};
};

struct RandomLocalState : public FunctionLocalState {
	explicit RandomLocalState(uint64_t seed) : random_engine(0) {
		random_engine.SetSeed(seed);
	}

	RandomEngine random_engine;
	uint64_t position = 0;
	uint64_t stream_index = 0;
	optional_idx query_id;

	void InitializeQuery(ClientContext &context) {
		auto current_query = context.transaction.GetActiveQuery();
		if (query_id.IsValid() && query_id.GetIndex() == current_query) {
			return;
		}
		random_engine.SetSeed(RandomQueryState::GetSeed(context, current_query));
		random_engine.Advance(stream_index * uint64_t(2));
		auto seed = random_engine.NextRandomInteger64();
		random_engine.SetSeed(seed);
		position = 0;
		query_id = current_query;
	}
};

void RandomFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = ExecuteFunctionState::GetFunctionState(state)->Cast<RandomLocalState>();
	result.SetVectorType(VectorType::FLAT_VECTOR);
	auto result_data = FlatVector::Writer<double>(result, args.size());
	if (args.ColumnCount() == 0) {
		for (idx_t i = 0; i < args.size(); i++) {
			result_data.WriteValue(lstate.random_engine.NextRandom());
		}
		return;
	}

	D_ASSERT(args.ColumnCount() == 2);
	lstate.InitializeQuery(state.GetContext());
	auto ordinals = args.data[0].Values<int64_t>();
	for (auto ordinal : ordinals) {
		auto position = NumericCast<uint64_t>(ordinal.GetValue() - 1);
		if (position != lstate.position) {
			// Unsigned subtraction and the PCG period wrap modulo 2^64. A double consumes two outputs.
			lstate.random_engine.Advance((position - lstate.position) * uint64_t(2));
		}
		result_data.WriteValue(lstate.random_engine.NextRandom());
		lstate.position = position + 1;
	}
}

unique_ptr<FunctionLocalState> RandomInitLocalState(ExpressionState &state, const BoundFunctionExpression &expr,
                                                    FunctionData *bind_data) {
	auto &context = state.GetContext();
	if (!expr.GetChildren().empty()) {
		D_ASSERT(expr.GetChildren().size() == 2);
		auto result = make_uniq<RandomLocalState>(0);
		result->stream_index = expr.GetChildren()[1]->Cast<BoundConstantExpression>().GetValue().GetValue<uint64_t>();
		return result;
	}
	if (context.transaction.HasActiveTransaction()) {
		return make_uniq<RandomLocalState>(
		    RandomQueryState::GetSeed(context, context.transaction.GetActiveQuery(), true));
	}
	auto &random_engine = RandomEngine::Get(context);
	lock_guard<mutex> guard(random_engine.lock);
	return make_uniq<RandomLocalState>(random_engine.NextRandomInteger64());
}

void GenerateUUIDv4Function(DataChunk &args, ExpressionState &state, Vector &result) {
	D_ASSERT(args.ColumnCount() == 0);
	auto &lstate = ExecuteFunctionState::GetFunctionState(state)->Cast<RandomLocalState>();

	result.SetVectorType(VectorType::FLAT_VECTOR);
	auto result_data = FlatVector::Writer<hugeint_t>(result, args.size());
	for (idx_t i = 0; i < args.size(); i++) {
		result_data.WriteValue(UUIDv4::GenerateRandomUUID(lstate.random_engine));
	}
}

void GenerateUUIDv7Function(DataChunk &args, ExpressionState &state, Vector &result) {
	D_ASSERT(args.ColumnCount() == 0);
	auto &lstate = ExecuteFunctionState::GetFunctionState(state)->Cast<RandomLocalState>();

	result.SetVectorType(VectorType::FLAT_VECTOR);
	auto result_data = FlatVector::Writer<hugeint_t>(result, args.size());
	for (idx_t i = 0; i < args.size(); i++) {
		result_data.WriteValue(UUIDv7::GenerateRandomUUID(lstate.random_engine));
	}
}

} // namespace

ScalarFunction RandomFun::GetFunction() {
	ScalarFunction random("random", {}, LogicalType::DOUBLE, RandomFunction, nullptr, nullptr, RandomInitLocalState);
	random.SetStability(FunctionStability::VOLATILE);
	return random;
}

ScalarFunction UUIDFun::GetFunction() {
	return UUIDv4Fun::GetFunction();
}

ScalarFunction UUIDv4Fun::GetFunction() {
	ScalarFunction uuid_v4_function({}, LogicalType::UUID, GenerateUUIDv4Function, nullptr, nullptr,
	                                RandomInitLocalState);
	// generate a random uuid v4
	uuid_v4_function.SetStability(FunctionStability::VOLATILE);
	return uuid_v4_function;
}

ScalarFunction UUIDv7Fun::GetFunction() {
	ScalarFunction uuid_v7_function({}, LogicalType::UUID, GenerateUUIDv7Function, nullptr, nullptr,
	                                RandomInitLocalState);
	// generate a random uuid v7
	uuid_v7_function.SetStability(FunctionStability::VOLATILE);
	return uuid_v7_function;
}

ScalarFunction UUIDExtractVersionFun::GetFunction() {
	ScalarFunction function({}, LogicalType::UINTEGER, ExtractVersionFunction<hugeint_t, ExtractVersionUuidOperator>);
	function.GetSignature().AddParameter("uuid", LogicalType::UUID);
	return function;
}

ScalarFunction UUIDExtractTimestampFun::GetFunction() {
	ScalarFunction function({}, LogicalType::TIMESTAMP_TZ,
	                        ExtractTimestampFunction<hugeint_t, ExtractTimestampUuidOperator>);
	function.GetSignature().AddParameter("uuid", LogicalType::UUID);
	// throws if the UUID is not a version 7 UUID
	function.SetFallible();
	return function;
}

} // namespace duckdb
