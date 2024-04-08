#define DUCKDB_EXTENSION_MAIN
#include "mimalloc_extension.hpp"

#include "duckdb/common/allocator.hpp"
// #include "mimalloc-new-delete.h" FIXME: we would like to override all allocations but it's tricky
#include "mimalloc.h"

namespace duckdb {

void MimallocExtension::Load(DuckDB &db) {
	// NOTE: This extension can only be loaded statically

	// These options can help free up memory: https://github.com/microsoft/mimalloc/issues/351
	mi_option_enable(mi_option_t::mi_option_purge_decommits);
	mi_option_enable(mi_option_t::mi_option_abandoned_page_purge);

	// 100ms is the default
	mi_option_set(mi_option_t::mi_option_purge_delay, 100);
}

std::string MimallocExtension::Name() {
	return "mimalloc";
}

data_ptr_t MimallocExtension::Allocate(PrivateAllocatorData *private_data, idx_t size) {
	return data_ptr_cast(mi_malloc(size));
}

void MimallocExtension::Free(PrivateAllocatorData *private_data, data_ptr_t pointer, idx_t size) {
	mi_free(pointer);
}

data_ptr_t MimallocExtension::Reallocate(PrivateAllocatorData *private_data, data_ptr_t pointer, idx_t old_size,
                                         idx_t size) {
	return data_ptr_cast(mi_realloc(pointer, size));
}

bool MiSumFree(const mi_heap_t *heap, const mi_heap_area_t *area, void *block, size_t block_size, void *arg) {
	auto &total_free = *reinterpret_cast<idx_t *>(arg);
	total_free += area->committed - area->used;
	return true;
}

void MimallocExtension::ThreadFlush(idx_t threshold) {
	auto heap = mi_heap_get_backing();
	idx_t total_free = 0;
	mi_heap_visit_blocks(heap, false, MiSumFree, &total_free);
	if (total_free > threshold) {
		mi_heap_collect(heap, false);
		mi_heap_collect(heap, true);
	}
}

void MimallocExtension::FlushAll() {
	auto heap = mi_heap_get_default();
	mi_heap_collect(heap, false);
	mi_heap_collect(heap, true);
}

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void mimalloc_init(duckdb::DatabaseInstance &db) {
	duckdb::DuckDB db_wrapper(db);
	db_wrapper.LoadExtension<duckdb::MimallocExtension>();
}

DUCKDB_EXTENSION_API const char *mimalloc_version() {
	return duckdb::DuckDB::LibraryVersion();
}
}

#ifndef DUCKDB_EXTENSION_MAIN
#error DUCKDB_EXTENSION_MAIN not defined
#endif
