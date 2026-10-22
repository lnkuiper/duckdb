//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/enums/state_initialization.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/constants.hpp"

namespace duckdb {

//! Whether state initialization consumes mutable client state, such as RNG seeds.
enum class StateInitialization : uint8_t { UNKNOWN, INDEPENDENT, SHARED };

} // namespace duckdb
