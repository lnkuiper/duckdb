//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/map.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/container_allocator.hpp"
#include "duckdb/common/pair.hpp"

#include <map>

namespace duckdb {

template <class _Key, class _Tp, class _Compare = std::less<_Key>,
          class _Allocator = container_allocator<pair<const _Key, _Tp>>>
using map = typename std::map<_Key, _Tp, _Compare, _Allocator>;

template <class _Key, class _Tp, class _Compare = std::less<_Key>,
          class _Allocator = container_allocator<pair<const _Key, _Tp>>>
using multimap = typename std::multimap<_Key, _Tp, _Compare, _Allocator>;

} // namespace duckdb
