//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/unordered_set.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/container_allocator.hpp"

#include <unordered_set>

namespace duckdb {

template <class _Value, class _Hash = std::hash<_Value>, class _Pred = std::equal_to<_Value>,
          class _Alloc = container_allocator<_Value>>
using unordered_set = typename std::unordered_set<_Value, _Hash, _Pred, _Alloc>;

} // namespace duckdb
