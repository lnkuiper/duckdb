//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/list.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/container_allocator.hpp"

#include <list>

namespace duckdb {

template <class _Tp, class _Alloc = container_allocator<_Tp>>
using list = typename std::list<_Tp, _Alloc>;

}
