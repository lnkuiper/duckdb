//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/execution/operator_chunk.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/execution/batch_metadata.hpp"

namespace duckdb {

//! Execution metadata travels with the rows, including when execution is suspended.
struct OperatorChunk {
	DataChunk chunk;
	ChunkMetadata metadata;

	void Reset() {
		chunk.Reset();
		metadata = ChunkMetadata();
	}
	void Move(OperatorChunk &other) {
		Move(other.chunk, other.metadata);
	}
	void Move(DataChunk &other, ChunkMetadata &other_metadata) {
		chunk.Move(other);
		metadata = other_metadata;
		other_metadata = ChunkMetadata();
	}
	void MoveTo(DataChunk &target, ChunkMetadata &target_metadata) {
		target.Move(chunk);
		target_metadata = metadata;
		metadata = ChunkMetadata();
	}
	void Reference(OperatorChunk &other) {
		chunk.Reference(other.chunk);
		metadata = other.metadata;
	}
};

} // namespace duckdb
