//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/stack.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/deque.hpp"

#include <stack>

namespace duckdb {

template <class _Tp, class _Container = deque<_Tp>>
using stack = typename std::stack<_Tp, _Container>;

}
