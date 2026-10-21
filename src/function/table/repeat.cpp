#include "duckdb/function/table/range.hpp"
#include "duckdb/function/table/generated_rows.hpp"

namespace duckdb {

struct RepeatFunctionData : public TableFunctionData {
	RepeatFunctionData(Value value, idx_t target_count) : value(std::move(value)), target_count(target_count) {
	}

	Value value;
	idx_t target_count;
};

static unique_ptr<FunctionData> RepeatBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<Identifier> &names) {
	// the repeat function returns the type of the first argument
	auto &inputs = input.inputs;
	return_types.push_back(inputs[0].type());
	names.emplace_back(inputs[0].ToString());
	if (inputs[1].IsNull()) {
		throw BinderException("Repeat second parameter cannot be NULL");
	}
	auto repeat_count = inputs[1].GetValue<int64_t>();
	if (repeat_count < 0) {
		throw BinderException("Repeat second parameter cannot be be less than 0");
	}
	return make_uniq<RepeatFunctionData>(inputs[0], NumericCast<idx_t>(repeat_count));
}

static unique_ptr<GlobalTableFunctionState> RepeatInit(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<RepeatFunctionData>();
	return make_uniq<GeneratedRowsGlobalState>(context, Hugeint::Convert(bind_data.target_count));
}

static void RepeatFunction(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<RepeatFunctionData>();
	auto &global_state = data_p.global_state->Cast<GeneratedRowsGlobalState>();
	auto &local_state = data_p.local_state->Cast<GeneratedRowsLocalState>();
	if (!local_state.Next(global_state)) {
		return;
	}
	data_p.row_offset = Hugeint::Cast<idx_t>(local_state.offset);
	output.data[0].Reference(bind_data.value, count_t(local_state.count));
	output.SetChildCardinality(local_state.count);
}

static unique_ptr<NodeStatistics> RepeatCardinality(ClientContext &context, const FunctionData *bind_data_p) {
	auto &bind_data = bind_data_p->Cast<RepeatFunctionData>();
	return make_uniq<NodeStatistics>(bind_data.target_count, bind_data.target_count);
}

void RepeatTableFunction::RegisterFunction(BuiltinFunctions &set) {
	TableFunction repeat("repeat",
	                     FunctionSignature()
	                         .AddPositionalOnly("value", LogicalType::ANY)
	                         .AddPositionalOnly("count", LogicalType::BIGINT),
	                     RepeatFunction, RepeatBind, RepeatInit, GeneratedRowsLocalState::Init);
	repeat.cardinality = RepeatCardinality;
	repeat.table_scan_progress = GeneratedRowsGlobalState::Progress;
	repeat.get_partition_data = GeneratedRowsLocalState::GetPartitionData;
	repeat.provides_row_offsets = true;
	set.AddFunction(repeat);
}

} // namespace duckdb
