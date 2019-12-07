//===----------------------------------------------------------------------===//
//                         DuckDB
//
// reoptimizer/reoptimizer.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/planner/logical_operator.hpp"

namespace duckdb {
class ClientContext;

class ReOptimizer {
public:
	ReOptimizer(Binder &binder, ClientContext &context);

	ClientContext &context;
    Binder &binder;
};

} // namespace duckdb
