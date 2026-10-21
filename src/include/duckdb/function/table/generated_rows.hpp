//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/function/table/generated_rows.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/atomic.hpp"
#include "duckdb/common/types/hugeint.hpp"
#include "duckdb/execution/partition_info.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_config.hpp"
#include "duckdb/storage/storage_info.hpp"

namespace duckdb {

//! Allocate ordered batches without materializing the generated rows.
struct GeneratedRowsGlobalState : public GlobalTableFunctionState {
	GeneratedRowsGlobalState(ClientContext &context, hugeint_t row_count_p)
	    : row_count(row_count_p),
	      batch_size(ClientConfig::GetConfig(context).verify_parallelism ? STANDARD_VECTOR_SIZE
	                                                                     : DEFAULT_ROW_GROUP_SIZE),
	      batch_count(Hugeint::Cast<idx_t>((row_count + batch_size - 1) / batch_size)) {
		D_ASSERT(row_count >= 0);
	}

	const hugeint_t row_count;
	const idx_t batch_size;
	const idx_t batch_count;
	atomic<idx_t> next_batch {0};
	atomic<idx_t> completed_batches {0};

	idx_t MaxThreads() const override {
		return MaxValue<idx_t>(batch_count, 1);
	}

	static double Progress(ClientContext &, const FunctionData *, const GlobalTableFunctionState *global_state) {
		if (!global_state) {
			return -1;
		}
		auto &state = global_state->Cast<GeneratedRowsGlobalState>();
		return state.batch_count == 0 ? 100
		                              : 100.0 * double(state.completed_batches.load(std::memory_order_relaxed)) /
		                                    double(state.batch_count);
	}
};

struct GeneratedRowsLocalState : public LocalTableFunctionState {
	optional_idx batch_index;
	hugeint_t offset = 0;
	idx_t count = 0;
	idx_t remaining = 0;

	bool Next(GeneratedRowsGlobalState &global_state) {
		D_ASSERT(remaining >= count);
		offset += count;
		remaining -= count;
		count = 0;
		if (remaining == 0) {
			if (batch_index.IsValid()) {
				// Report progress once per batch rather than synchronizing for every vector.
				global_state.completed_batches.fetch_add(1, std::memory_order_relaxed);
				batch_index = optional_idx();
			}
			auto next_batch = global_state.next_batch.fetch_add(1, std::memory_order_relaxed);
			if (next_batch >= global_state.batch_count) {
				return false;
			}
			batch_index = next_batch;
			offset = Hugeint::Convert(next_batch) * global_state.batch_size;
			remaining =
			    Hugeint::Cast<idx_t>(MinValue<hugeint_t>(global_state.row_count - offset, global_state.batch_size));
		}
		count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);
		return true;
	}

	static unique_ptr<LocalTableFunctionState> Init(ExecutionContext &, TableFunctionInitInput &,
	                                                GlobalTableFunctionState *) {
		return make_uniq<GeneratedRowsLocalState>();
	}

	static OperatorPartitionData GetPartitionData(ClientContext &, TableFunctionGetPartitionInput &input) {
		auto &state = input.local_state->Cast<GeneratedRowsLocalState>();
		return OperatorPartitionData(state.batch_index.GetIndex());
	}
};

} // namespace duckdb
