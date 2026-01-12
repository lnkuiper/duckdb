//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/optimizer/partitioned_execution.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/planner/logical_operator.hpp"

namespace duckdb {

class Optimizer;

class PartitionedExecution {
public:
	explicit PartitionedExecution(Optimizer &optimizer);

public:
	unique_ptr<LogicalOperator> Optimize(unique_ptr<LogicalOperator> op);

private:
	Optimizer &optimizer;
};

} // namespace duckdb
