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
	void Reference(OperatorChunk &other) {
		chunk.Reference(other.chunk);
		metadata = other.metadata;
	}
};

} // namespace duckdb
