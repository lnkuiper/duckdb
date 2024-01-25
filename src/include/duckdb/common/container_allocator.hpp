//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/allocator.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/allocator.hpp"

namespace duckdb {

template <typename _Tp>
struct container_allocator : public std::allocator<_Tp> {
public:
	using original = std::allocator<_Tp>;
	using original::original;

	typename original::pointer allocate(typename original::size_type n, const void * = 0) {
		return reinterpret_cast<typename original::pointer>(Allocator::DefaultAllocator().AllocateData(n));
	}

	void deallocate(typename original::pointer p, typename original::size_type n) {
		Allocator::DefaultAllocator().FreeData(data_ptr_cast(p), n);
	}

	template <typename _Up>
	struct rebind {
		typedef container_allocator<_Up> other;
	};
};

} // namespace duckdb
