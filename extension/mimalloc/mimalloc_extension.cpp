#define DUCKDB_EXTENSION_MAIN
#include "mimalloc_extension.hpp"

#include "duckdb/common/allocator.hpp"
// #include "mimalloc-new-delete.h" FIXME: we would like to override all allocations but it's tricky
#include "duckdb/common/pair.hpp"
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

data_ptr_t MimallocExtension::Allocate(PrivateAllocatorData *, idx_t size) {
	return data_ptr_cast(mi_malloc(size));
}

void MimallocExtension::Free(PrivateAllocatorData *, data_ptr_t pointer, idx_t) {
	mi_free(pointer);
}

data_ptr_t MimallocExtension::Reallocate(PrivateAllocatorData *, data_ptr_t pointer, idx_t, idx_t size) {
	return data_ptr_cast(mi_realloc(pointer, size));
}

bool SumCommitedAndUsed(const mi_heap_t *, const mi_heap_area_t *area, void *, size_t, void *arg) {
	auto &committed_and_used = *reinterpret_cast<pair<idx_t, idx_t> *>(arg);
	committed_and_used.first += area->committed;
	committed_and_used.second += area->used;
	return true;
}

void MimallocExtension::ThreadFlush(idx_t threshold) {
	auto heap = mi_heap_get_backing();

	auto committed_and_used = make_pair<idx_t, idx_t>(0, 0);
	mi_heap_visit_blocks(heap, false, SumCommitedAndUsed, &committed_and_used);
	auto &committed = committed_and_used.first;
	auto &used = committed_and_used.second;

	if (used > threshold) {
		// This thread has more than threshold in use after finishing task, most likely buffer-managed blocks
		// Delete the heap so that the main thread's heap takes ownership over them
		mi_heap_delete(heap);
	} else if (committed - used > threshold) {
		// This thread has more than theshold outstanding unused allocations, clean them up
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
