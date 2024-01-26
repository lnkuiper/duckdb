//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/string.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/container_allocator.hpp"

#include <sstream>
#include <string>

namespace duckdb {

typedef std::basic_string<char, std::char_traits<char>, container_allocator<char>> string;

} // namespace duckdb
