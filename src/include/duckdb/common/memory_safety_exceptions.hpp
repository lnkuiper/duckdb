//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/memory_safety_exceptions.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/typedefs.hpp"

#include <stdexcept>

namespace duckdb {

std::exception VectorOutOfBoundsException(idx_t index, idx_t size);
std::exception VectorBackOnEmptyException();

std::exception UniquePtrNullException();
std::exception OptionalPtrNullException();

} // namespace duckdb
