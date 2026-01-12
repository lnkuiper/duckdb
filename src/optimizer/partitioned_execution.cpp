#include "duckdb/optimizer/partitioned_execution.hpp"

namespace duckdb {

PartitionedExecution::PartitionedExecution(Optimizer &optimizer_p) : optimizer(optimizer_p) {
}

unique_ptr<LogicalOperator> PartitionedExecution::Optimize(unique_ptr<LogicalOperator> op) {
	return op;
}

} // namespace duckdb
