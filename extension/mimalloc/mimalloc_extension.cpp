#define DUCKDB_EXTENSION_MAIN
#include "mimalloc_extension.hpp"

#include "duckdb/common/allocator.hpp"
#include "mimalloc.h"

namespace duckdb {

static const idx_t DEFAULT_PURGE_DELAY = 1000;

void MimallocExtension::Load(DuckDB &db) {
	// NOTE: This extension can only be loaded statically

	// These options can help free up memory: https://github.com/microsoft/mimalloc/issues/351
	// However, memory retention seems fine without them
	// mi_option_set_enabled_default(mi_option_t::mi_option_purge_decommits, true);
	// mi_option_set_enabled_default(mi_option_t::mi_option_abandoned_page_purge, true);
}

std::string MimallocExtension::Name() {
	return "mimalloc";
}

data_ptr_t MimallocExtension::Allocate(PrivateAllocatorData *, idx_t size) {
	return data_ptr_cast(mi_malloc(size));
}

void MimallocExtension::Free(PrivateAllocatorData *, data_ptr_t pointer, idx_t) {
	mi_free(pointer);
}

data_ptr_t MimallocExtension::Reallocate(PrivateAllocatorData *, data_ptr_t pointer, idx_t, idx_t size) {
	return data_ptr_cast(mi_realloc(pointer, size));
}

void MimallocExtension::ThreadFlush(idx_t) {
	// mimalloc's cleanup assumes threads are short-lived, i.e., do one task, but our threads live forever
	// this gets us the behavior from mimalloc that we want
	FlushAll();
	mi_thread_done();
}

void FlushHeap(mi_heap_t *heap) {
	mi_heap_collect(heap, false);
	mi_heap_collect(heap, true);
}

void MimallocExtension::FlushAll() {
	mi_option_set(mi_option_t::mi_option_purge_delay, 0);
	FlushHeap(mi_heap_get_default());
	FlushHeap(mi_heap_get_backing());
	mi_option_set(mi_option_t::mi_option_purge_delay, DEFAULT_PURGE_DELAY);
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
