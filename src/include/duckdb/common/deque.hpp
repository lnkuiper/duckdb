//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/deque.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/container_allocator.hpp"

#include <deque>

namespace duckdb {

template <class _Tp, class _Allocator = container_allocator<_Tp>>
using deque = typename std::deque<_Tp, _Allocator>;

}
