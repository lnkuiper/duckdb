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
			break; // Just continue
		}
		default:
			return nullptr; // Unsupported for partition pass-through
		}
		child_ref = *child_ref.get().children[0];
	}
	D_ASSERT(child_ref.get().type == LogicalOperatorType::LOGICAL_GET);
	if (!child_ref.get().children.empty()) {
		return nullptr; // Table in/out, bail
	}
	return child_ref.get().Cast<LogicalGet>();
}

void PartitionedExecution::Optimize(unique_ptr<LogicalOperator> &op) {
	for (auto &child : op->children) {
		Optimize(child); // Depth-first
	}

	if (op->children.size() != 1) {
		return; // Only for straight pipelines to scans
	}

	vector<PartitionedExecutionColumn> columns;
	if (!PartitionedExecutionGetColumns(*op, columns) || columns.empty()) {
		return; // Not able to get partition columns from this operator
	}

	optional_ptr<LogicalGet> get = PartitionedExecutionTraceColumns(*op, columns);
	if (!get || columns.empty()) {
		return; // Unable to trace bindings down to scan
	}

	GetPartitionStatsInput input(get->function, get->bind_data.get());
	auto partition_stats = get->function.get_partition_stats(optimizer.context, input);
}

} // namespace duckdb
