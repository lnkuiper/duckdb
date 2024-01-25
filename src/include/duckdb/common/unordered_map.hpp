//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/unordered_map.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/container_allocator.hpp"
#include "duckdb/common/pair.hpp"

#include <unordered_map>

namespace duckdb {

template <class _Key, class _Tp, class _Hash = std::hash<_Key>, class _Pred = std::equal_to<_Key>,
          class _Alloc = container_allocator<pair<const _Key, _Tp>>>
using unordered_map = typename std::unordered_map<_Key, _Tp, _Hash, _Pred, _Alloc>;

} // namespace duckdb
