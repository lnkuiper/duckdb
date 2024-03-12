#include "duckdb/common/container_allocator.hpp"

#include "duckdb/common/allocator.hpp"

namespace duckdb {

data_ptr_t AllocatorWrapper::Allocate(idx_t size) {
	return Allocator::DefaultAllocate(nullptr, size);
}

void AllocatorWrapper::Free(data_ptr_t pointer, idx_t size) {
	Allocator::DefaultFree(nullptr, pointer, size);
}

} // namespace duckdb
