#include "duckdb/execution/batch_metadata.hpp"
#include "duckdb/common/exception.hpp"

namespace duckdb {

bool ChunkMetadata::CanAppend(const ChunkMetadata &other, idx_t count) const {
	if (batch.batch_index != other.batch.batch_index || batch.stream_id != other.batch.stream_id ||
	    batch.row_count != other.batch.row_count || batch.first_row_index != other.batch.first_row_index ||
	    offset_in_batch.IsValid() != other.offset_in_batch.IsValid()) {
		return false;
	}
	return !offset_in_batch.IsValid() || (other.offset_in_batch.GetIndex() >= offset_in_batch.GetIndex() &&
	                                      other.offset_in_batch.GetIndex() - offset_in_batch.GetIndex() == count);
}

void BatchMetadata::Verify(const BatchMetadataInfo &properties) const {
	if (!properties.Any()) {
		return;
	}
	if (!batch_index.IsValid() || !stream_id.IsValid()) {
		throw InternalException("Operator did not provide promised batch identity");
	}
	if (properties.row_positions && !first_row_index.IsValid()) {
		throw InternalException("Operator did not provide promised batch row position");
	}
}

void ChunkMetadata::Verify(const BatchMetadataInfo &properties, idx_t count) const {
	if (!properties.Any() || count == 0) {
		return;
	}
	batch.Verify(properties);
	if (!properties.row_positions) {
		return;
	}
	if (!offset_in_batch.IsValid()) {
		throw InternalException("Operator did not provide promised row positions");
	}
	if (batch.row_count.IsValid() && (offset_in_batch.GetIndex() > batch.row_count.GetIndex() ||
	                                  count > batch.row_count.GetIndex() - offset_in_batch.GetIndex())) {
		throw InternalException("Operator chunk exceeds its declared batch cardinality");
	}
}

} // namespace duckdb
