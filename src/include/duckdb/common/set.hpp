//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/set.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/container_allocator.hpp"

#include <set>

namespace duckdb {

template <class _Key, class _Compare = std::less<_Key>, class _Allocator = container_allocator<_Key>>
using set = typename std::set<_Key, _Compare, _Allocator>;

template <class _Key, class _Compare = std::less<_Key>, class _Allocator = container_allocator<_Key>>
using multiset = typename std::multiset<_Key, _Compare, _Allocator>;

} // namespace duckdb
