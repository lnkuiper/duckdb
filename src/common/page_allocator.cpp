#include "duckdb/common/page_allocator.hpp"

#if defined(_WIN32)
#include "duckdb/common/windows.hpp"
#else
#include <sys/mman.h>
#endif

namespace duckdb {

PageArena::PageArena(const idx_t &arena_idx_p)
    : arena_idx(arena_idx_p), allocation(NewPage()), count(PageAllocator::PAGES_PER_ALLOCATION) {
}

data_ptr_t PageArena::Allocate() {
	// lock_guard<mutex> guard(lock);
	// if (count == PageAllocator::PAGES_PER_ALLOCATION) {
	// 	allocation = NewPage();
	// 	count = 0;
	// }
	return allocation + count++ * PageAllocator::PAGE_SIZE;
}

data_ptr_t PageArena::NewPage() {
	void *allocation;
#if defined(_WIN32)
	allocation = VirtualAlloc(nullptr, largePageSize, MEM_RESERVE | MEM_COMMIT | MEM_LARGE_PAGES, PAGE_READWRITE);
#else
	allocation =
	    mmap(nullptr, PageAllocator::ALLOCATION_SIZE * 100, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#endif
	return data_ptr_cast(allocation);
}

PageAllocator::PageAllocator(const idx_t &num_arenas) : fast_mod(num_arenas), next_thread_idx(0) {
	arenas.reserve(num_arenas);
	for (idx_t arena_idx = 0; arena_idx < num_arenas; arena_idx++) {
		arenas.emplace_back(make_unsafe_uniq<PageArena>(arena_idx));
	}
}

data_ptr_t PageAllocator::Allocate() {
	return GetArena().Allocate();
}

void PageAllocator::Free(const data_ptr_t &ptr) {
	// TODO
}

PageArena &PageAllocator::GetArena() {
	thread_local idx_t thread_idx = DConstants::INVALID_INDEX;
	if (thread_idx == DConstants::INVALID_INDEX) {
		thread_idx = next_thread_idx++;
	}
	return *arenas[fast_mod.Mod(thread_idx)];
}

} // namespace duckdb
