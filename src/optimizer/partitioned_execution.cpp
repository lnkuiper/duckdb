#include "duckdb/optimizer/partitioned_execution.hpp"

#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/planner/operator/list.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

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

static bool PartitionedExecutionGetColumns(LogicalOperator &op, vector<PartitionedExecutionColumn> &columns) {
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

static vector<PartitionedExecutionColumn>::iterator
PartitionedExecutionHandleColumnRemoval(const LogicalOperatorType type, vector<PartitionedExecutionColumn> &columns,
                                        vector<PartitionedExecutionColumn>::iterator &it) {
	switch (type) {
	case LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY:
		return columns.erase(it); // We can partition on any colref so we can freely remove any
	case LogicalOperatorType::LOGICAL_ORDER_BY:
		columns.erase(it, columns.end()); // Have to remove from the first non-colref onward
		return columns.end();
	default:
		throw NotImplementedException("PartitionedExecutionHandleColumnRemoval for %s", EnumUtil::ToString(type));
	}
}

static optional_ptr<LogicalGet> PartitionedExecutionTraceColumns(LogicalOperator &op,
                                                                 vector<PartitionedExecutionColumn> &columns) {
	if (op.children.size() != 1) {
		return nullptr; // Can't handle more than one child (yet)
	}

	reference<LogicalOperator> child_ref(*op.children[0]);
	while (child_ref.get().type != LogicalOperatorType::LOGICAL_GET) {
		switch (child_ref.get().type) {
		case LogicalOperatorType::LOGICAL_PROJECTION: {
			auto &proj = child_ref.get().Cast<LogicalProjection>();
			for (auto it = columns.begin(); it != columns.end();) {
				auto &expr = *proj.expressions[it->column_binding.column_index];
				if (expr.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
					it->column_binding = expr.Cast<BoundColumnRefExpression>().binding;
					it++;
				} else {
					it = PartitionedExecutionHandleColumnRemoval(op.type, columns, it);
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
	for (auto it = columns.begin(); it != columns.end();) {
		const auto &column_index = column_ids[it->column_binding.column_index];
		if (get.TryGetStorageIndex(column_index, it->storage_index)) {
			it++;
		} else {
			it = PartitionedExecutionHandleColumnRemoval(op.type, columns, it); // Did not get a storage index
		}
	}

	return get;
}

struct PartitionedExecutionSplitNode {
	PartitionedExecutionSplitNode(Value &&val_p, int64_t count_delta_p)
	    : val(std::move(val_p)), count_delta(count_delta_p) {
	}
	Value val;
	int64_t count_delta;
};

struct PartitionedExecutionSplit {
	BaseStatistics &min;
	BaseStatistics &max;
	idx_t count;
};

static bool PartitionedExecutionCanUseStats(const unique_ptr<BaseStatistics> &stats) {
	if (!stats || stats->CanHaveNull()) {
		return false; // No stats or contains NULL
	}
	switch (stats->GetStatsType()) {
	case StatisticsType::NUMERIC_STATS:
		return NumericStats::HasMinMax(*stats);
	case StatisticsType::STRING_STATS:
		return true;
	default:
		return false; // Only numeric/string supported for now
	}
}

static void PartitionedExecutionAddSplitNodes(const BaseStatistics &stats, const idx_t count,
                                              vector<PartitionedExecutionSplitNode> &split_nodes) {
	switch (stats.GetStatsType()) {
	case StatisticsType::NUMERIC_STATS:
		split_nodes.emplace_back(NumericStats::Min(stats), NumericCast<int64_t>(count));
		split_nodes.emplace_back(NumericStats::Max(stats), -NumericCast<int64_t>(count));
		break;
	case StatisticsType::STRING_STATS:
		split_nodes.emplace_back(StringStats::Min(stats), NumericCast<int64_t>(count));
		split_nodes.emplace_back(StringStats::Max(stats), -NumericCast<int64_t>(count));
		break;
	default:
		throw NotImplementedException("PartitionedExecutionAddSplitNodes for %s",
		                              EnumUtil::ToString(stats.GetStatsType()));
	}
}

static bool
PartitionedExecutionCompareSplitNodes(const vector<PartitionedExecutionColumn> &columns,
                                      const vector<vector<PartitionedExecutionSplitNode>> &column_split_nodes,
                                      const idx_t &lhs, const idx_t &rhs) {
	idx_t col_idx = 0;
	for (; col_idx < columns.size(); col_idx++) {
		const auto &split_nodes = column_split_nodes[col_idx];
		const auto &lhs_node = split_nodes[lhs];
		const auto &rhs_node = split_nodes[rhs];
		if (lhs_node.val != rhs_node.val) {
			break;
		}
	}

	if (col_idx == columns.size()) {
		return false; // All columns equal
	}

	const auto &split_nodes = column_split_nodes[col_idx];
	const auto &lhs_node = split_nodes[lhs];
	const auto &rhs_node = split_nodes[rhs];

	switch (columns[col_idx].order_type) {
	case OrderType::ASCENDING:
	case OrderType::INVALID:
		return lhs_node.val < rhs_node.val;
	case OrderType::DESCENDING:
		return lhs_node.val > rhs_node.val;
	default:
		throw NotImplementedException("PartitionedExecutionCompareSplitNodes for %s",
		                              EnumUtil::ToString(columns[col_idx].order_type));
	}
}

static vector<PartitionedExecutionSplit>
PartitionedExecutionComputeSplits(LogicalOperator &op, vector<PartitionedExecutionColumn> &columns,
                                  const vector<PartitionStatistics> &partition_stats, const idx_t num_threads) {
	vector<vector<PartitionedExecutionSplitNode>> column_split_nodes;
	for (auto it = columns.begin(); it != columns.end();) {
		vector<PartitionedExecutionSplitNode> split_nodes;
		split_nodes.reserve(partition_stats.size() * 2);

		bool success = true;
		for (auto &ps : partition_stats) {
			if (!ps.partition_row_group) {
				success = false; // No row group to get stats from
				break;
			}
			const auto stats = ps.partition_row_group->GetColumnStatistics(it->storage_index);
			if (!PartitionedExecutionCanUseStats(stats)) {
				success = false; // Unable to use these stats
				break;
			}
			PartitionedExecutionAddSplitNodes(*stats, ps.count, split_nodes);
		}
		if (success) {
			column_split_nodes.emplace_back(std::move(split_nodes));
			it++;
		} else {
			it = PartitionedExecutionHandleColumnRemoval(op.type, columns, it);
		}
		// TODO: Maybe impose a maximum number of columns, don't want to compare many columns in case of "BY ALL"
		//  or, try to split by one column at a time and stop once we're happy with the splits
	}
	D_ASSERT(columns.size() == column_split_nodes.size());

	if (columns.empty()) {
		return {}; // None of the columns were eligible
	}

	// Initialize indices into the split nodes
	vector<idx_t> indices;
	indices.reserve(partition_stats.size() * 2);
	for (idx_t i = 0; i < partition_stats.size() * 2; i++) {
		indices.emplace_back(i);
	}

	// Sort indices based on the values in the split nodes
	std::sort(indices.begin(), indices.end(), [&columns, &column_split_nodes](const idx_t &lhs, const idx_t &rhs) {
		return PartitionedExecutionCompareSplitNodes(columns, column_split_nodes, lhs, rhs);
	});

	// Compute how the count goes up and down at the split points
	vector<int64_t> counts;
	counts.reserve(indices.size());
	int64_t count = 0;
	for (const auto &index : indices) {
		counts.emplace_back(count);
		const auto &split_node = column_split_nodes[0][index];
		count += split_node.count_delta;
	}

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

	const auto num_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(optimizer.context).NumberOfThreads());
	const auto splits = PartitionedExecutionComputeSplits(*op, columns, partition_stats, num_threads);
	if (splits.empty()) {
		return; // Unable to compute splits
	}

	PartitionExecutionSplitPipeline(optimizer, root, op, splits); // Success!
}

} // namespace duckdb
