#include "duckdb/optimizer/partitioned_execution.hpp"

#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/planner/operator/list.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"

namespace duckdb {

PartitionedExecution::PartitionedExecution(Optimizer &optimizer_p, LogicalOperator &root_p)
    : optimizer(optimizer_p), root(root_p) {
}

struct PartitionedExecutionColumn {
	explicit PartitionedExecutionColumn(Expression &expr, OrderType order_type_p = OrderType::INVALID)
	    : column_binding(expr.Cast<BoundColumnRefExpression>().binding), order_type(order_type_p) {
	}

	ColumnBinding column_binding;
	OrderType order_type = OrderType::INVALID;
	StorageIndex storage_index;
};

bool PartitionedExecutionGetColumns(LogicalOperator &op, vector<PartitionedExecutionColumn> &columns) {
	switch (op.type) {
	case LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY: {
		const auto &agg = op.Cast<LogicalAggregate>();
		if (agg.grouping_sets.size() > 1 || !agg.grouping_functions.empty() || agg.groups.empty()) {
			return false; // Only regular grouped aggregations
		}
		for (auto &group : agg.groups) {
			if (group->GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
				columns.emplace_back(*group); // We can partition on any colref
			}
		}
		return true;
	}
	case LogicalOperatorType::LOGICAL_ORDER_BY: {
		auto &order = op.Cast<LogicalOrder>();
		for (auto &order_by_node : order.orders) {
			if (order_by_node.expression->GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF) {
				break; // Have to break on the first non-colref
			}
			columns.emplace_back(*order_by_node.expression, order_by_node.type);
		}
		return true;
	}
	default:
		// TODO: WINDOW/ASOF
		return false;
	}
}

optional_ptr<LogicalGet> PartitionedExecutionTraceColumns(LogicalOperator &op,
                                                          vector<PartitionedExecutionColumn> &columns) {
	reference<LogicalOperator> child_ref(op);
	while (child_ref.get().type != LogicalOperatorType::LOGICAL_GET) {
		switch (child_ref.get().type) {
		case LogicalOperatorType::LOGICAL_PROJECTION: {
			auto &proj = child_ref.get().Cast<LogicalProjection>();
			for (auto it = columns.begin(); it != columns.end(); it++) {
				auto &expr = *proj.expressions[it->column_binding.column_index];
				if (expr.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
					it->column_binding = expr.Cast<BoundColumnRefExpression>().binding;
					continue;
				}
				switch (op.type) {
				case LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY:
					it = columns.erase(it); // We can partition on any colref so we can freely remove any
					break;
				case LogicalOperatorType::LOGICAL_ORDER_BY:
					for (; it != columns.end(); it++) {
						it = columns.erase(it); // Have to remove from the first non-colref onward
					}
					break;
				default:
					throw NotImplementedException("PartitionedExecutionTraceColumns for %s",
					                              EnumUtil::ToString(op.type));
				}
			}
			break;
		}
		case LogicalOperatorType::LOGICAL_FILTER: {
			break; // Don't have to update bindings
		}
		default:
			return nullptr; // Unsupported for partition pass-through
		}
		child_ref = *child_ref.get().children[0];
	}

	D_ASSERT(child_ref.get().type == LogicalOperatorType::LOGICAL_GET);
	if (!child_ref.get().children.empty()) {
		return nullptr; // Table in/out, unsupported
	}
	auto &get = child_ref.get().Cast<LogicalGet>();

	// Get the storage index
	const auto &column_ids = get.GetColumnIds();
	for (auto it = columns.begin(); it != columns.end(); it++) {
		const auto &column_index = column_ids[it->column_binding.column_index];
		if (get.TryGetStorageIndex(column_index, it->storage_index)) {
			continue; // Successfully got a storage index
		}
		switch (op.type) {
		case LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY:
			it = columns.erase(it); // We can partition on any column so we can freely remove any
			break;
		case LogicalOperatorType::LOGICAL_ORDER_BY:
			for (; it != columns.end(); it++) {
				it = columns.erase(it); // Have to remove from the first column without storage index onwards
			}
			break;
		default:
			throw NotImplementedException("PartitionedExecutionTraceColumns for %s", EnumUtil::ToString(op.type));
		}
	}

	return get;
}

struct PartitionedExecutionSplit {
	BaseStatistics &min;
	BaseStatistics &max;
	idx_t count;
};

vector<PartitionedExecutionSplit>
PartitionedExecutionComputeSplits(const vector<PartitionedExecutionColumn> &columns,
                                  const vector<PartitionStatistics> &partition_stats) {
	vector<PartitionedExecutionSplit> result;
	throw NotImplementedException("PartitionedExecutionComputeSplits");
	return result;
}

void PartitionExecutionSplitPipeline(Optimizer &optimizer, LogicalOperator &root, unique_ptr<LogicalOperator> &op,
                                     const vector<PartitionedExecutionSplit> &splits) {
	throw NotImplementedException("PartitionExecutionSplitPipeline");
}

void PartitionedExecution::Optimize(unique_ptr<LogicalOperator> &op) {
	for (auto &child : op->children) {
		Optimize(child); // Depth-first
	}

	vector<PartitionedExecutionColumn> columns;
	if (!PartitionedExecutionGetColumns(*op, columns) || columns.empty()) {
		return; // Unable to get partition columns from this operator
	}

	optional_ptr<LogicalGet> get = PartitionedExecutionTraceColumns(*op, columns);
	if (!get || columns.empty()) {
		return; // Unable to trace any binding down to scan
	}

	GetPartitionStatsInput input(get->function, get->bind_data.get());
	const auto partition_stats = get->function.get_partition_stats(optimizer.context, input);
	if (partition_stats.size() <= 1) {
		return; // Can't split 0 or 1 partitions
	}

	const auto splits = PartitionedExecutionComputeSplits(columns, partition_stats);
	if (splits.empty()) {
		return; // Unable to compute splits
	}

	PartitionExecutionSplitPipeline(optimizer, root, op, splits); // Success!
}

} // namespace duckdb
