//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/function/table/generated_rows.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/atomic.hpp"
#include "duckdb/common/mutex.hpp"
#include "duckdb/common/types/hugeint.hpp"
#include "duckdb/execution/partition_info.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_config.hpp"
#include "duckdb/storage/storage_info.hpp"

namespace duckdb {

struct GeneratedRowsProgress {
	atomic<double> rows {0};
};

//! Allocate ordered batches without materializing the generated rows.
struct GeneratedRowsGlobalState : public GlobalTableFunctionState {
	GeneratedRowsGlobalState(ClientContext &context, hugeint_t row_count_p)
	    : row_count(row_count_p), batch_size(ClientConfig::GetConfig(context).verify_parallelism
	                                             ? MaxValue<idx_t>(DEFAULT_STANDARD_VECTOR_SIZE, STANDARD_VECTOR_SIZE)
	                                             : DEFAULT_ROW_GROUP_SIZE),
	      batch_count(
	          Hugeint::Cast<idx_t>((row_count + Hugeint::Convert(batch_size) - 1) / Hugeint::Convert(batch_size))) {
		D_ASSERT(row_count >= 0);
	}

	const hugeint_t row_count;
	const idx_t batch_size;
	const idx_t batch_count;
	atomic<idx_t> next_batch {0};
	mutable mutex progress_lock;
	vector<unique_ptr<GeneratedRowsProgress>> local_progress;

	idx_t MaxThreads() const override {
		return MaxValue<idx_t>(batch_count, 1);
	}

	static double Progress(ClientContext &, const FunctionData *, const GlobalTableFunctionState *global_state) {
		if (!global_state) {
			return -1;
		}
		auto &state = global_state->Cast<GeneratedRowsGlobalState>();
		if (state.row_count == 0) {
			return 100;
		}
		lock_guard<mutex> guard(state.progress_lock);
		double rows = 0;
		for (auto &progress : state.local_progress) {
			rows += progress->rows.load(std::memory_order_relaxed);
		}
		return MinValue(100.0, 100.0 * rows / Hugeint::Cast<double>(state.row_count));
	}
};

struct GeneratedRowsLocalState : public LocalTableFunctionState {
	explicit GeneratedRowsLocalState(GeneratedRowsProgress &progress_p) : progress(progress_p) {
	}

	GeneratedRowsProgress &progress;
	hugeint_t completed_rows = 0;
	optional_idx batch_index;
	hugeint_t offset = 0;
	idx_t count = 0;
	idx_t remaining = 0;

	bool Next(GeneratedRowsGlobalState &global_state) {
		D_ASSERT(remaining >= count);
		offset += Hugeint::Convert(count);
		completed_rows += Hugeint::Convert(count);
		progress.rows.store(Hugeint::Cast<double>(completed_rows), std::memory_order_relaxed);
		remaining -= count;
		count = 0;
		if (remaining == 0) {
			if (batch_index.IsValid()) {
				batch_index = optional_idx();
			}
			auto next_batch = global_state.next_batch.fetch_add(1, std::memory_order_relaxed);
			if (next_batch >= global_state.batch_count) {
				return false;
			}
			batch_index = next_batch;
			offset = Hugeint::Convert(next_batch) * Hugeint::Convert(global_state.batch_size);
			remaining = Hugeint::Cast<idx_t>(
			    MinValue<hugeint_t>(global_state.row_count - offset, Hugeint::Convert(global_state.batch_size)));
		}
		count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);
		return true;
	}

	static unique_ptr<LocalTableFunctionState> Init(ExecutionContext &, TableFunctionInitInput &,
	                                                GlobalTableFunctionState *global_state) {
		auto &state = global_state->Cast<GeneratedRowsGlobalState>();
		auto progress = make_uniq<GeneratedRowsProgress>();
		auto result = make_uniq<GeneratedRowsLocalState>(*progress);
		lock_guard<mutex> guard(state.progress_lock);
		state.local_progress.push_back(std::move(progress));
		return std::move(result);
	}

	static OperatorPartitionData GetPartitionData(ClientContext &, TableFunctionGetPartitionInput &input) {
		auto &state = input.local_state->Cast<GeneratedRowsLocalState>();
		return OperatorPartitionData(state.batch_index.GetIndex());
	}
};

} // namespace duckdb
