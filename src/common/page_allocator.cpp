#include "duckdb/common/page_allocator.hpp"

#if defined(_WIN32)
#include "duckdb/common/windows.hpp"
#else
#include <sys/mman.h>
#endif

namespace duckdb {

PageAllocator::PageAllocator(const idx_t &num_arenas)
    : fast_mod(MinValue<idx_t>(num_arenas, MAX_ARENAS)), arenas(fast_mod.GetDivisor()) {
}

data_ptr_t PageAllocator::Allocate() {
	const auto arena_idx = GetThreadArenaIndex();
	lock_guard<mutex> guard(arena_locks[arena_idx]);
	auto &arena = arenas[arena_idx];
	if (!arena) {
		arena = make_unsafe_uniq<PageAllocatorArena>(arena_idx);
	}
	return arena.get()->Allocate(*this);
}

uint8_t PageAllocator::GetThreadArenaIndex() const {
	static constexpr idx_t ARENA_INDEX_ALLOCATION_REFRESH_INTERVAL = 8;
	thread_local idx_t allocation_count = ARENA_INDEX_ALLOCATION_REFRESH_INTERVAL - 1;
	thread_local idx_t cpu_id = DConstants::INVALID_INDEX;
	if (++allocation_count == ARENA_INDEX_ALLOCATION_REFRESH_INTERVAL) {
		cpu_id = TaskScheduler::GetEstimatedCPUId();
		allocation_count = 0;
	}
	return UnsafeNumericCast<uint8_t>(fast_mod.Mod(cpu_id));
}

void PageAllocator::Free(const data_ptr_t &ptr) {
	arenas[PointerToArenaIndex(ptr)]->Free(ptr);
}

idx_t PageAllocator::PointerArenaLookupTableIndex(const data_ptr_t &ptr) {
	return cast_pointer_to_uint64(ptr) / CHUNK_SIZE;
}

uint8_t PageAllocator::PointerToArenaIndex(const data_ptr_t &ptr) const {
	return pointer_arena_lookup_table[PointerArenaLookupTableIndex(ptr)];
}

PageAllocatorArena::PageAllocatorArena(const uint8_t &arena_idx_p) : arena_idx(arena_idx_p) {
}

data_ptr_t PageAllocatorArena::Allocate(PageAllocator &page_allocator) {
	for (auto &pool : pools) {
		const auto ptr = pool->Allocate(page_allocator, arena_idx);
		if (ptr) {
			return ptr;
		}
	}
	pools.emplace_back(make_unsafe_uniq<PageAllocatorPool>());
	return pools.back()->Allocate(page_allocator, arena_idx);
}

void PageAllocatorArena::Free(const data_ptr_t &ptr) {
	for (idx_t pool_idx = 0; pool_idx < pools.size(); pool_idx++) {
		auto &pool = pools[pool_idx];
		if (!pool->Free(ptr)) {
			continue; // Allocated in a different pool
		}
		// Move pool with most mapped chunks to the front to reduce fragmentation
		if (pool_idx > 0 && pool->GetMappedChunks() < pools[pool_idx - 1]->GetMappedChunks()) {
			std::swap(pools[pool_idx], pools[pool_idx - 1]);
		}
		return;
	}
	throw InternalException("Unable to free page!");
}

PageAllocatorPool::PageAllocatorPool() : first_free_chunk(0), mapped_chunks(0), allocated_pages(0) {
	for (idx_t i = 0; i < CHUNKS_PER_POOL; i++) {
		// NULL for sanity
		chunks[i] = nullptr;

		// Initialize all chunks as empty
		chunk_occupancies[i] = ChunkOccupancy();

		// Start off with normal order
		ordered_chunk_idxs[i] = NumericCast<uint8_t>(i);
		chunk_idx_to_order[i] = NumericCast<uint8_t>(i);
	}

	for (idx_t i = 0; i < HASH_TABLE_SIZE; i++) {
		hash_table[i].occupied = false; // Initialize HT as empty
	}
}

data_ptr_t PageAllocatorPool::Allocate(PageAllocator &page_allocator, const uint8_t &arena_idx) {
	if (first_free_chunk == CHUNKS_PER_POOL) {
		return nullptr; // Pool is full
	}

	const auto &chunk_idx = ordered_chunk_idxs[first_free_chunk];
	auto &chunk = chunks[chunk_idx];
	auto &chunk_occupancy = chunk_occupancies[chunk_idx];

	if (chunk_occupancy.IsEmpty()) {
		// Map a new chunk
		MapChunk(chunk_idx);

		// Add it to the global lookup table
		page_allocator.pointer_arena_lookup_table[PageAllocator::PointerArenaLookupTableIndex(chunk)] = arena_idx;
	}

	// Mark the first free page in this chunk as occupied
	const auto page_idx = chunk_occupancy.GetFirstFreePage();
	chunk_occupancy.SetPageOccupied(page_idx);
	allocated_pages++;

	if (chunk_occupancy.IsFull()) {
		first_free_chunk++; // Chunk is full, increment free index
	}

	Verify();
	return chunk + page_idx * PageAllocator::PAGE_SIZE;
}

bool PageAllocatorPool::Free(const data_ptr_t &ptr) {
	const auto &entry = GetHTEntry(ptr);
	if (!entry.occupied) {
		return false; // Pointer is not in this pool
	}

	const auto &chunk_idx = entry.chunk_idx;
	auto &chunk = chunks[chunk_idx];
	auto &chunk_occupancy = chunk_occupancies[chunk_idx];

	// Derive page index in this chunk, then mark as unoccupied
	const auto page_idx = NumericCast<idx_t>(ptr - chunk) / PageAllocator::PAGE_SIZE;
	chunk_occupancy.SetPageFree(page_idx);
	allocated_pages--;

	if (chunk_occupancy.IsEmpty()) {
		UnmapChunk(chunk_idx);
	}

	// We keep the chunk indices in sorted order to avoid fragmentation
	const auto occupied_page_count = chunk_occupancy.OccupiedPageCount();
	uint8_t idx_in_order_before = chunk_idx_to_order[chunk_idx];
	D_ASSERT(ordered_chunk_idxs[idx_in_order_before] == chunk_idx);

	// Find the first chunk that has less than or the same number of occupied pages
	idx_t idx_in_order_after = idx_in_order_before + 1;
	for (; idx_in_order_after < CHUNKS_PER_POOL; idx_in_order_after++) {
		const auto &next_chunk_idx = ordered_chunk_idxs[idx_in_order_after];
		if (chunk_occupancies[next_chunk_idx].OccupiedPageCount() <= occupied_page_count) {
			break;
		}
	}
	idx_in_order_after--;

	// Swap the indices
	std::swap(ordered_chunk_idxs[idx_in_order_before], ordered_chunk_idxs[idx_in_order_after]);
	std::swap(chunk_idx_to_order[idx_in_order_before], chunk_idx_to_order[idx_in_order_after]);

	if (first_free_chunk == CHUNKS_PER_POOL || idx_in_order_after < first_free_chunk) {
		first_free_chunk = idx_in_order_after; // Update which chunk
	}


	Verify();
	return true;
}

void PageAllocatorPool::MapChunk(const uint8_t &chunk_idx) {
	mapped_chunks++;
	auto &chunk = chunks[chunk_idx];

	if (chunk) {
		return; // We already mapped this before, but told the OS that we don't need it
	}

	// Windows/Linux support huge pages, so the allocation should be aligned to 2 MiB
	// MacOS does not support huge pages, so we have to align the allocation
	void *allocation;
#if defined(_WIN32)
	allocation =
	    VirtualAlloc(nullptr, PageAllocator::CHUNK_SIZE, MEM_RESERVE | MEM_COMMIT | MEM_LARGE_PAGES, PAGE_READWRITE);
#elif defined(__APPLE__)
	allocation =
	    mmap(nullptr, PageAllocator::CHUNK_SIZE * 2, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	allocation = AlignPointer<PageAllocator::CHUNK_SIZE>(allocation);
#else
	allocation = mmap(nullptr, PageAllocator::CHUNK_SIZE, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#endif
	D_ASSERT(allocation == AlignPointer<PageAllocator::CHUNK_SIZE>(allocation));

	chunk = data_ptr_cast(allocation);
	InsertHTEntry(chunk, chunk_idx);
}

void PageAllocatorPool::UnmapChunk(const uint8_t &chunk_idx) {
	mapped_chunks--;
	auto &chunk = chunks[chunk_idx];
	D_ASSERT(chunk);

	// Tell the OS that it can reclaim the pages without actually munmapping them
#if defined(_WIN32)
	VirtualAlloc(chunk, PageAllocator::CHUNK_SIZE, MEM_RESET, PAGE_READWRITE);
#else
	madvise(chunk, PageAllocator::CHUNK_SIZE, MADV_DONTNEED);
#endif
}

idx_t PageAllocatorPool::GetMappedChunks() const {
	return mapped_chunks;
}

idx_t PageAllocatorPool::GetHTIndex(const data_ptr_t &ptr, const bool &find_unoccupied) const {
	// Round pointer to chunk size so we get the chunk pointer, then lookup in HT
	const auto chunk_ptr_hashable = AlignValueFloor<uint64_t, PageAllocator::CHUNK_SIZE>(cast_pointer_to_uint64(ptr));
	idx_t index = Hash(chunk_ptr_hashable) % HASH_TABLE_SIZE;
	const auto chunk_ptr = cast_uint64_to_pointer(chunk_ptr_hashable);
	if (find_unoccupied) {
		// Return the first unoccupied entry
		for (; true; ++index %= HASH_TABLE_SIZE) { // Linear probing
			const auto &entry = hash_table[index];
			if (!entry.occupied) {
				break;
			}
		}
	} else {
		// Find the corresponding pointer
		for (; true; ++index %= HASH_TABLE_SIZE) { // Linear probing
			const auto &entry = hash_table[index];
			D_ASSERT(entry.occupied);
			if (chunk_ptr == chunks[entry.chunk_idx]) {
				break;
			}
		}
	}
	return index;
}

void PageAllocatorPool::InsertHTEntry(const data_ptr_t &ptr, const uint8_t &chunk_idx) {
	auto &entry = hash_table[GetHTIndex(ptr, true)];
	D_ASSERT(!entry.occupied);
	entry.occupied = true;
	entry.chunk_idx = chunk_idx;
}

const PageAllocatorPool::PoolHTEntry &PageAllocatorPool::GetHTEntry(const data_ptr_t &ptr) const {
	return hash_table[GetHTIndex(ptr, false)];
}

void PageAllocatorPool::Verify() const {
#ifdef DEBUG
	bool first_free_chunk_found = false;
	idx_t mapped_chunks_verification = 0;
	idx_t allocated_pages_verification = 0;
	for (idx_t i = 0; i < CHUNKS_PER_POOL; i++) {
		const auto &chunk_idx = ordered_chunk_idxs[i];
		const auto &chunk = chunks[chunk_idx];
		const auto &chunk_occupancy = chunk_occupancies[chunk_idx];

		// Count statistics
		mapped_chunks_verification += !chunk_occupancy.IsEmpty();
		allocated_pages_verification += chunk_occupancy.OccupiedPageCount();

		// Verify that the first free chunk index actually points to the first non-full chunk
		if (!first_free_chunk_found && !chunk_occupancy.IsFull()) {
			D_ASSERT(i == first_free_chunk);
			first_free_chunk_found = true;
		}

		// Verify that the reverse mapping is valid
		D_ASSERT(chunk_idx_to_order[chunk_idx] == i);

		// Verify that the chunk indices are ordered by page occupancy
		if (i < CHUNKS_PER_POOL - 1) {
			const auto next_chunk_idx = ordered_chunk_idxs[i + 1];
			const auto &next_chunk_occupancy = chunk_occupancies[next_chunk_idx];
			D_ASSERT(chunk_occupancy.OccupiedPageCount() >= next_chunk_occupancy.OccupiedPageCount());
		}

		// Verify the hash table
		if (chunk) {
			const auto &entry = GetHTEntry(chunk);
			D_ASSERT(entry.occupied && chunk_idx == entry.chunk_idx);
		}
	}

	// The statistics should be valid
	D_ASSERT(mapped_chunks_verification == mapped_chunks);
	D_ASSERT(allocated_pages_verification == allocated_pages);
#endif
}

} // namespace duckdb
