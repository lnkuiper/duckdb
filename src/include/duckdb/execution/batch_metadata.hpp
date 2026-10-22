//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/execution/batch_metadata.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/optional_idx.hpp"

namespace duckdb {

//! Requested or guaranteed properties of a stream.
struct BatchMetadataInfo {
	bool batch_index = false;
	bool row_positions = false;
	bool completion = false;

	static BatchMetadataInfo BatchIndex() {
		BatchMetadataInfo result;
		result.batch_index = true;
		return result;
	}
	static BatchMetadataInfo RowPositions() {
		auto result = BatchIndex();
		result.row_positions = true;
		return result;
	}
	static BatchMetadataInfo Completion() {
		auto result = BatchIndex();
		result.completion = true;
		return result;
	}
	bool Any() const {
		return batch_index || row_positions || completion;
	}
	bool Supports(const BatchMetadataInfo &requested) const {
		return (!requested.batch_index || batch_index) && (!requested.row_positions || row_positions) &&
		       (!requested.completion || completion);
	}
	BatchMetadataInfo Intersect(const BatchMetadataInfo &other) const {
		BatchMetadataInfo result;
		result.batch_index = batch_index && other.batch_index;
		result.row_positions = row_positions && other.row_positions;
		result.completion = completion && other.completion;
		return result;
	}
	void Merge(const BatchMetadataInfo &other) {
		batch_index = batch_index || other.batch_index;
		row_positions = row_positions || other.row_positions;
		completion = completion || other.completion;
	}
};

struct BatchMetadata {
	optional_idx batch_index;
	//! Identifies the execution stream in which row ordinals are defined.
	optional_idx stream_id;
	optional_idx row_count;
	optional_idx first_row_index;

	void Verify(const BatchMetadataInfo &properties) const;
};

//! Positions, when present, describe a contiguous interval at this operator boundary.
struct ChunkMetadata {
	BatchMetadata batch;
	optional_idx offset_in_batch;

	void ClearPositions() {
		batch.first_row_index.SetInvalid();
		offset_in_batch.SetInvalid();
	}
	void Restrict(const BatchMetadataInfo &properties) {
		if (!properties.batch_index) {
			*this = ChunkMetadata();
		} else if (!properties.row_positions) {
			ClearPositions();
		}
	}
	void Verify(const BatchMetadataInfo &properties, idx_t count) const;
	bool CanAppend(const ChunkMetadata &other, idx_t count) const;
};

} // namespace duckdb
