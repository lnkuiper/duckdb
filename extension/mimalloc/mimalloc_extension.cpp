#define DUCKDB_EXTENSION_MAIN
#include "mimalloc_extension.hpp"

#include "duckdb/common/allocator.hpp"
#include "mimalloc-new-delete.h"
#include "mimalloc.h"

namespace duckdb {

void MimallocExtension::Load(DuckDB &db) {
	// NOTE: This extension can only be loaded statically

	// These options should help free up memory: https://github.com/microsoft/mimalloc/issues/351
	mi_option_enable(mi_option_t::mi_option_purge_decommits);
	mi_option_enable(mi_option_t::mi_option_abandoned_page_purge);
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

void MimallocExtension::FlushAll() {
	mi_collect(false);
	mi_collect(true);
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
