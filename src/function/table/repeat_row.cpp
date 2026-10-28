#include "duckdb/function/table/range.hpp"
#include "duckdb/function/table/generated_rows.hpp"

namespace duckdb {

struct RepeatRowFunctionData : public TableFunctionData {
	RepeatRowFunctionData(vector<Value> values, idx_t target_count)
	    : values(std::move(values)), target_count(target_count) {
	}

	const vector<Value> values;
	idx_t target_count;
};

static unique_ptr<FunctionData> RepeatRowBind(ClientContext &context, TableFunctionBindInput &input,
                                              vector<LogicalType> &return_types, vector<Identifier> &names) {
	auto &inputs = input.inputs;
	for (idx_t input_idx = 0; input_idx < inputs.size(); input_idx++) {
		return_types.push_back(inputs[input_idx].type());
		names.emplace_back("column" + std::to_string(input_idx));
	}
	auto &num_rows = input.named_parameters.at("num_rows");
	if (inputs.empty()) {
		throw BinderException("repeat_rows requires at least one column to be specified");
	}
	if (num_rows.IsNull()) {
		throw BinderException("num_rows should be an integer value >= 0");
	}
	auto repeat_rows = num_rows.GetValue<int64_t>();
	if (repeat_rows < 0) {
		throw BinderException("num_rows cannot be less than zero");
	}
	return make_uniq<RepeatRowFunctionData>(inputs, NumericCast<idx_t>(repeat_rows));
}

static unique_ptr<GlobalTableFunctionState> RepeatRowInit(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<RepeatRowFunctionData>();
	return make_uniq<GeneratedRowsGlobalState>(context, Hugeint::Convert(bind_data.target_count), input.column_ids);
}

static void RepeatRowFunction(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<RepeatRowFunctionData>();
	auto &global_state = data_p.global_state->Cast<GeneratedRowsGlobalState>();
	auto &local_state = data_p.local_state->Cast<GeneratedRowsLocalState>();
	if (!local_state.Next(global_state)) {
		return;
	}
	for (idx_t col = 0; col < global_state.column_ids.size(); col++) {
		auto id = global_state.column_ids[col];
		if (id == COLUMN_IDENTIFIER_ROW_NUMBER) {
			local_state.RowNumber(output.data[col]);
		} else {
			output.data[col].Reference(bind_data.values[id], count_t(local_state.count));
		}
	}
	output.SetChildCardinality(local_state.count);
}

static unique_ptr<NodeStatistics> RepeatRowCardinality(ClientContext &context, const FunctionData *bind_data_p) {
	auto &bind_data = bind_data_p->Cast<RepeatRowFunctionData>();
	return make_uniq<NodeStatistics>(bind_data.target_count, bind_data.target_count);
}

void RepeatRowTableFunction::RegisterFunction(BuiltinFunctions &set) {
	TableFunction repeat_row("repeat_row", {}, RepeatRowFunction, RepeatRowBind, RepeatRowInit,
	                         GeneratedRowsLocalState::Init);
	// a positional list, so it declares "*args" and no "**kwargs" - an argument named after no parameter is an error
	repeat_row.GetSignature().AddArgs("args", LogicalType::ANY);
	repeat_row.GetSignature().AddKeywordOnly("num_rows", LogicalType::BIGINT);
	repeat_row.cardinality = RepeatRowCardinality;
	repeat_row.table_scan_progress = GeneratedRowsGlobalState::Progress;
	repeat_row.get_partition_data = GeneratedRowsLocalState::GetPartitionData;
	repeat_row.get_virtual_columns = GeneratedRowsLocalState::GetVirtualColumns;
	repeat_row.projection_pushdown = true;
	repeat_row.supports_row_number = true;
	set.AddFunction(repeat_row);
}

} // namespace duckdb
