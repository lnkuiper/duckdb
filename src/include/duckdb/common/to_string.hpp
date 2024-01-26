//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/to_string.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/string.hpp"

namespace duckdb {

template <class T>
string to_string(const T &val) {
	return string(std::to_string(val));
}

} // namespace duckdb
