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

string to_string(int __val);
string to_string(unsigned __val);
string to_string(long __val);
string to_string(unsigned long __val);
string to_string(long long __val);
string to_string(unsigned long long __val);
string to_string(float __val);
string to_string(double __val);
string to_string(long double __val);

} // namespace duckdb
