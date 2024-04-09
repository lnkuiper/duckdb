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

struct MiHeapStats {
	idx_t used = 0;
	idx_t committed = 0;
};

static inline bool SumHeapStats(const mi_heap_t *, const mi_heap_area_t *area, void *, size_t, void *arg) {
	auto &stats = *reinterpret_cast<MiHeapStats *>(arg);
	stats.used += area->used * area->block_size;
	stats.committed += area->committed;
	return true;
}

static inline void FlushHeap(mi_heap_t *heap, idx_t threshold) {
	MiHeapStats stats;
	mi_heap_visit_blocks(heap, false, SumHeapStats, &stats);
	Printer::PrintF("%llu before: %llu/%llu", heap, stats.used, stats.committed);

	if (stats.used > threshold) {
		// This thread has more than threshold in use after finishing task, most likely buffer-managed blocks
		// Delete the heap so that the main thread's heap takes ownership over them
		mi_thread_done();
		mi_thread_init();
//		mi_heap_destroy(heap);
	} else if (stats.committed - stats.used > threshold) {
		// This thread has more than theshold outstanding unused allocations, clean them up
		mi_heap_collect(heap, false);
		mi_heap_collect(heap, true);
	}

	stats.used = 0;
	stats.committed = 0;
	mi_heap_visit_blocks(heap, false, SumHeapStats, &stats);
	Printer::PrintF("%llu after: %llu/%llu", heap, stats.used, stats.committed);
}

void MimallocExtension::ThreadFlush(idx_t threshold) {
	FlushHeap(mi_heap_get_default(), threshold);
	FlushHeap(mi_heap_get_backing(), threshold);
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
