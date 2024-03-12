//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/allocator.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/winapi.hpp"

#include <memory>

namespace duckdb {

//! Wrapper to prevent circular imports
struct AllocatorWrapper {
	DUCKDB_API static data_ptr_t Allocate(idx_t size);
	DUCKDB_API static void Free(data_ptr_t pointer);
};

template <typename _Tp>
struct container_allocator : public std::allocator<_Tp> {
public:
	using original = std::allocator<_Tp>;
	using original::original;

	typename original::pointer allocate(typename original::size_type n, const void * = 0) {
		return reinterpret_cast<typename original::pointer>(AllocatorWrapper::Allocate(n * sizeof(_Tp)));
	}

	void deallocate(typename original::pointer p, typename original::size_type n) {
		AllocatorWrapper::Free(data_ptr_cast(p));
	}

	template <typename _Up>
	struct rebind {
		typedef container_allocator<_Up> other;
	};
};

} // namespace duckdb
