//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/memory_safety_exceptions.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/winapi.hpp"

#include <stdexcept>

namespace duckdb {

DUCKDB_API std::exception VectorOutOfBoundsException(idx_t index, idx_t size);
DUCKDB_API std::exception VectorBackOnEmptyException();

DUCKDB_API std::exception UniquePtrNullException();
DUCKDB_API std::exception OptionalPtrNullException();

} // namespace duckdb
