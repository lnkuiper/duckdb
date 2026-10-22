#include "duckdb/execution/batch_metadata.hpp"
#include "duckdb/common/exception.hpp"

namespace duckdb {

void ChunkMetadata::Verify(const BatchMetadataInfo &properties, idx_t count) const {
	if (!properties.Any() || count == 0) {
		return;
	}
	if (!batch.batch_index.IsValid() || !batch.stream_id.IsValid()) {
		throw InternalException("Operator did not provide promised batch identity");
	}
	if (!properties.row_positions) {
		return;
	}
	if (!batch.first_row_index.IsValid() || !offset_in_batch.IsValid()) {
		throw InternalException("Operator did not provide promised row positions");
	}
	if (batch.row_count.IsValid() && (offset_in_batch.GetIndex() > batch.row_count.GetIndex() ||
	                                  count > batch.row_count.GetIndex() - offset_in_batch.GetIndex())) {
		throw InternalException("Operator chunk exceeds its declared batch cardinality");
	}
}

} // namespace duckdb
