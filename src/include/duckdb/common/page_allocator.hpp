//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/page_allocator.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/atomic.hpp"
#include "duckdb/common/mutex.hpp"
#include "duckdb/storage/storage_info.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "duckdb/common/bit_utils.hpp"

namespace duckdb {

class PageAllocatorArena;
class PageAllocatorPool;

class PageAllocator {
	friend class PageAllocatorArena;
	friend class PageAllocatorPool;

public:
	explicit PageAllocator(const idx_t &num_arenas);

public:
	data_ptr_t Allocate();
	void Free(const data_ptr_t &ptr);

private:
	uint8_t GetThreadArenaIndex();
	static idx_t PointerArenaLookupTableIndex(const data_ptr_t &ptr);
	uint8_t PointerToArenaIndex(const data_ptr_t &ptr) const;

public:
	static constexpr idx_t PAGE_SIZE = DEFAULT_BLOCK_ALLOC_SIZE;
	static constexpr idx_t PAGES_PER_CHUNK = 8;
	static constexpr idx_t CHUNK_SIZE = PAGE_SIZE * PAGES_PER_CHUNK;

private:
	mutex lock;
	const FastMod<idx_t> fast_mod;
	atomic<idx_t> next_arena_index;

	static constexpr idx_t MAX_ARENAS = 256;
	mutex arena_locks[MAX_ARENAS];
	unsafe_vector<unsafe_unique_ptr<PageAllocatorArena>> arenas;

	static constexpr idx_t MAX_ADDRESS = 1ULL << 47;
	static constexpr idx_t POINTER_TO_ARENA_TABLE_SIZE = MAX_ADDRESS / CHUNK_SIZE;
	uint8_t pointer_arena_lookup_table[POINTER_TO_ARENA_TABLE_SIZE];
};

class PageAllocatorArena {
public:
	explicit PageAllocatorArena(const uint8_t &arena_idx);

public:
	data_ptr_t Allocate(PageAllocator &page_allocator);
	void Free(PageAllocator &page_allocator, const data_ptr_t &ptr);

private:
	const uint8_t arena_idx;
	vector<unsafe_unique_ptr<PageAllocatorPool>> pools;
};

class PageAllocatorPool {
public:
	PageAllocatorPool();

public:
	data_ptr_t Allocate(PageAllocator &page_allocator, const uint8_t &arena_idx);
	bool Free(PageAllocator &page_allocator, const uint8_t &arena_idx, const data_ptr_t &ptr);
	idx_t GetMappedChunks() const;

private:
	struct PoolHTEntry {
	public:
		PoolHTEntry() : occupied(false), chunk_idx(0) {
		}

	public:
		bool IsOccupied() const {
			return occupied;
		}

		void SetOccupied(const bool &occupied_p) {
			occupied = occupied_p;
		}

		uint16_t GetChunkIndex() const {
			return chunk_idx;
		}

		void SetChunkIndex(const uint16_t &chunk_idx_p) {
			D_ASSERT(chunk_idx_p < CHUNKS_PER_POOL);
			chunk_idx = chunk_idx_p;
		}

	private:
		bool occupied : 1;
		uint16_t chunk_idx : 15;
	};

	struct ChunkOccupancy {
	public:
		ChunkOccupancy() : bitmap(0) {
		}

	public:
		bool IsOccupied(const idx_t &page_idx) const {
			D_ASSERT(page_idx < PageAllocator::PAGES_PER_CHUNK);
			return bitmap & static_cast<uint8_t>(1) << page_idx;
		}

		void SetPageOccupied(const idx_t &page_idx) {
			D_ASSERT(!IsOccupied(page_idx));
			bitmap |= static_cast<uint8_t>(1) << page_idx;
			D_ASSERT(IsOccupied(page_idx));
		}

		void SetPageFree(const idx_t &page_idx) {
			D_ASSERT(IsOccupied(page_idx));
			bitmap &= ~(static_cast<uint8_t>(1) << page_idx);
			D_ASSERT(!IsOccupied(page_idx));
		}

		idx_t GetFirstFreePage() const {
			D_ASSERT(!IsFull()); // Should be handled before
			const auto page_index = CountZeros<uint8_t>::Trailing(~bitmap);
			D_ASSERT(page_index < PageAllocator::PAGES_PER_CHUNK);
			return page_index;
		}

		idx_t OccupiedPageCount() const {
			// Kernighan's algorithm
			idx_t occupied = 0;
			auto bitmap_copy = bitmap;
			while (bitmap_copy) {
				bitmap_copy &= (bitmap_copy - 1);
				occupied++;
			}
			return occupied;
		}

		bool IsEmpty() const {
			return bitmap == EMPTY;
		}

		bool IsFull() const {
			return bitmap == FULL;
		}

	private:
		static constexpr uint8_t EMPTY = 0;
		static constexpr uint8_t FULL = 255;

		uint8_t bitmap;
	};

private:
	void MapChunk(const uint16_t &chunk_idx);
	void UnmapChunk(const uint16_t &chunk_idx);

	idx_t GetHTIndex(const data_ptr_t &ptr, const bool &find_unoccupied) const;
	void InsertHTEntry(const data_ptr_t &ptr, const uint16_t &chunk_idx);
	const PoolHTEntry &GetHTEntry(const data_ptr_t &ptr) const;

	void Verify(const PageAllocator &page_allocator, const uint16_t &arena_idx) const;

private:
	static constexpr idx_t CHUNKS_PER_POOL = 32768;

	data_ptr_t chunks[CHUNKS_PER_POOL];
	ChunkOccupancy chunk_occupancies[CHUNKS_PER_POOL];

	//! Index into "ordered_chunk_idxs" with the first free chunk
	idx_t first_free_chunk;
	//! Chunk indices ordered by occupancy (fullest first)
	uint16_t ordered_chunk_idxs[CHUNKS_PER_POOL];
	//! Reverse mapping of "ordered_chunk_idxs"
	uint16_t chunk_idx_to_order[CHUNKS_PER_POOL];

	//! We use "chunk_occupancy" as a bitmap, where 0 means occupied, and 1 means free.
	//! We do this instead of the other way around so we can use CountZeros<uint8_t>::Leading,
	//! to find the first empty page slot in the chunk
	static constexpr uint16_t CHUNK_EMPTY = CHUNKS_PER_POOL - 1;
	static constexpr uint16_t CHUNK_FULL = 0;

	//! Mapping from pointer to chunk index
	static constexpr idx_t HASH_TABLE_SIZE = CHUNKS_PER_POOL * 2;
	PoolHTEntry hash_table[HASH_TABLE_SIZE];

	//! Statistics
	idx_t mapped_chunks;
	idx_t allocated_pages;
};

} // namespace duckdb
