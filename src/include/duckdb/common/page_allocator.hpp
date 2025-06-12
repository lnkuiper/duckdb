//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/page_allocator.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/mutex.hpp"
#include "duckdb/storage/storage_info.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "duckdb/common/atomic.hpp"

namespace duckdb {

class PageArena {
public:
	explicit PageArena(const idx_t &arena_idx);

public:
	data_ptr_t Allocate();

private:
	static data_ptr_t NewPage();

private:
	const idx_t arena_idx;

	mutex lock;
	data_ptr_t allocation;
	idx_t count;
};

class PageAllocator {
public:
	explicit PageAllocator(const idx_t &num_arenas);

public:
	data_ptr_t Allocate();
	void Free(const data_ptr_t &ptr);

private:
	PageArena &GetArena();

public:
	static constexpr idx_t PAGE_SIZE = DEFAULT_BLOCK_ALLOC_SIZE;
	static constexpr idx_t PAGES_PER_ALLOCATION = 8;
	static constexpr idx_t ALLOCATION_SIZE = PAGE_SIZE * PAGES_PER_ALLOCATION;

private:
	const FastMod<idx_t> fast_mod;
	unsafe_vector<unsafe_unique_ptr<PageArena>> arenas;
	atomic<idx_t> next_thread_idx;
};

} // namespace duckdb
