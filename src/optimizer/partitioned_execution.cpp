#include "duckdb/optimizer/partitioned_execution.hpp"

#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/optimizer/filter_pushdown.hpp"
#include "duckdb/planner/operator/list.hpp"
#include "duckdb/planner/logical_operator_deep_copy.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/expression/bound_window_expression.hpp"
#include "duckdb/optimizer/column_binding_replacer.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

namespace duckdb {

PartitionedExecution::PartitionedExecution(Optimizer &optimizer_p, LogicalOperator &root_p)
    : optimizer(optimizer_p), root(root_p),
      num_threads(NumericCast<idx_t>(TaskScheduler::GetScheduler(optimizer.context).NumberOfThreads())) {
}

struct PartitionedExecutionConfig {
	//! Maximum number of columns to use for splitting ranges
	static constexpr idx_t MAXIMUM_COLUMNS = 3;
	//! Minimum number of row groups per thread per partition to split on
	static constexpr idx_t MIN_ROW_GROUPS_PER_THREAD_PER_PARTITION = 16;
	//! Maximum overlap (as fraction of partition size) that we allow for a split
	static constexpr double MAX_OVERLAP_RATIO = 0.1;
	//! Minimum input cardinality before we even consider splitting
	static constexpr idx_t MINIMUM_INPUT_CARDINALITY = 4194304;
};

struct PartitionedExecutionColumn {
	explicit PartitionedExecutionColumn(idx_t original_idx_p, Expression &expr,
	                                    OrderType order_type_p = OrderType::INVALID)
	    : original_idx(original_idx_p), column_binding(expr.Cast<BoundColumnRefExpression>().binding),
	      order_type(order_type_p) {
	}

	idx_t original_idx;
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
				columns.emplace_back(columns.size(), *group); // We can partition on any colref
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
			columns.emplace_back(columns.size(), *order_by_node.expression, order_by_node.type);
		}
		return true;
	}
	case LogicalOperatorType::LOGICAL_WINDOW: {
		const auto &expr = op.expressions[0]->Cast<BoundWindowExpression>();
		for (idx_t expr_idx = 1; expr_idx < op.expressions.size(); expr_idx++) {
			if (!expr.PartitionsAreEquivalent(op.expressions[expr_idx]->Cast<BoundWindowExpression>())) {
				return false;
			}
		}
		for (auto &partition : expr.partitions) {
			if (partition->GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
				columns.emplace_back(columns.size(), *partition);
			}
		}
		return true;
	}
	default:
		return false;
	}
}

static vector<PartitionedExecutionColumn>::iterator
PartitionedExecutionHandleColumnRemoval(const LogicalOperatorType type, vector<PartitionedExecutionColumn> &columns,
                                        vector<PartitionedExecutionColumn>::iterator &it) {
	switch (type) {
	case LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY:
	case LogicalOperatorType::LOGICAL_WINDOW:
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

	if (!get.function.get_partition_stats) {
		return nullptr; // We need this
	}

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

struct PartitionedExecutionStatsNode {
	PartitionedExecutionStatsNode(Value &&val_p, int64_t count_delta_p)
	    : val(std::move(val_p)), count_delta(count_delta_p) {
	}
	Value val;
	int64_t count_delta;
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

static void PartitionedExecutionAddStatsNodes(const BaseStatistics &stats, const idx_t count,
                                              vector<PartitionedExecutionStatsNode> &stats_nodes) {
	switch (stats.GetStatsType()) {
	case StatisticsType::NUMERIC_STATS:
		stats_nodes.emplace_back(NumericStats::Min(stats), NumericCast<int64_t>(count));
		stats_nodes.emplace_back(NumericStats::Max(stats), -NumericCast<int64_t>(count));
		break;
	case StatisticsType::STRING_STATS:
		stats_nodes.emplace_back(StringStats::Min(stats), NumericCast<int64_t>(count));
		stats_nodes.emplace_back(StringStats::Max(stats), -NumericCast<int64_t>(count));
		break;
	default:
		throw NotImplementedException("PartitionedExecutionAddStatsNodes for %s",
		                              EnumUtil::ToString(stats.GetStatsType()));
	}
}

static vector<vector<PartitionedExecutionStatsNode>>
PartitionedExecutionCollectStatsNodes(LogicalOperator &op, vector<PartitionedExecutionColumn> &columns,
                                      const vector<PartitionStatistics> &partition_stats) {
	// We impose a maximum number of columns otherwise we could be fetching a massive amount of stats from storage
	vector<vector<PartitionedExecutionStatsNode>> column_stats_nodes;
	for (auto it = columns.begin();
	     it != columns.end() && column_stats_nodes.size() < PartitionedExecutionConfig::MAXIMUM_COLUMNS;) {
		vector<PartitionedExecutionStatsNode> stats_nodes;
		stats_nodes.reserve(partition_stats.size() * 2);

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
			PartitionedExecutionAddStatsNodes(*stats, ps.count, stats_nodes);
		}
		if (success) {
			column_stats_nodes.emplace_back(std::move(stats_nodes));
			it++;
		} else {
			it = PartitionedExecutionHandleColumnRemoval(op.type, columns, it);
		}
	}
	while (columns.size() > column_stats_nodes.size()) {
		columns.pop_back();
	}
	D_ASSERT(columns.size() == column_stats_nodes.size());
	return column_stats_nodes;
}

static bool
PartitionedExecutionCompareStatsNodes(const vector<PartitionedExecutionColumn> &columns,
                                      const vector<vector<PartitionedExecutionStatsNode>> &column_stats_nodes,
                                      const idx_t &lhs, const idx_t &rhs) {
	for (idx_t col_idx = 0; col_idx < columns.size(); col_idx++) {
		const auto &stats_nodes = column_stats_nodes[col_idx];
		const auto &lhs_node = stats_nodes[lhs];
		const auto &rhs_node = stats_nodes[rhs];
		if (col_idx < columns.size() - 1 && lhs_node.val == rhs_node.val) {
			continue; // Not the last iteration and values are equal
		}

		// Last iteration or values not equal, need to return
		if (lhs_node.val == rhs_node.val) {
			return lhs_node.count_delta > rhs_node.count_delta; // All columns equal, add start of row groups first
		}

		// Return based on sort order
		switch (columns[col_idx].order_type) {
		case OrderType::ASCENDING:
		case OrderType::INVALID:
			return lhs_node.val < rhs_node.val;
		case OrderType::DESCENDING:
			return lhs_node.val > rhs_node.val;
		default:
			throw NotImplementedException("PartitionedExecutionCompareStatsNodes for %s",
			                              EnumUtil::ToString(columns[col_idx].order_type));
		}
	}
	throw InternalException("PartitionedExecutionCompareStatsNodes failed"); // This should be unreachable
}

static void PartitionedExecutionGrowPartition(
    const vector<PartitionedExecutionStatsNode> &stats_nodes, const vector<idx_t> &indices, idx_t &i,
    int64_t &current_overlap, idx_t &partition_count, double &partition_overlap_ratio,
    const std::function<bool()> &stopping_criterium, const std::function<void()> &update_callback = [] {}) {
	while (i < indices.size() && !stopping_criterium()) {
		const auto &node = stats_nodes[indices[i++]];
		current_overlap += node.count_delta; // Keep track of overall overlap at "i"
		D_ASSERT(current_overlap >= 0);
		if (node.count_delta > 0) {
			partition_count += NumericCast<idx_t>(node.count_delta); // Only add if it's the start of a row group
		}
		partition_overlap_ratio = static_cast<double>(current_overlap) / static_cast<double>(partition_count);
		update_callback();
	}
}

struct PartitionedExecutionRange {
	vector<Value> min;
	vector<Value> max;
	idx_t overlap;
	idx_t count;
};

static vector<PartitionedExecutionRange>
PartitionedExecutionComputeRanges(LogicalOperator &op, vector<PartitionedExecutionColumn> &columns,
                                  const vector<PartitionStatistics> &partition_stats, const idx_t num_threads) {
	auto column_stats_nodes = PartitionedExecutionCollectStatsNodes(op, columns, partition_stats);
	if (columns.empty()) {
		return {}; // None of the columns were eligible
	}

	// Initialize indices into the stats nodes
	vector<idx_t> indices;
	indices.reserve(partition_stats.size() * 2);
	for (idx_t i = 0; i < partition_stats.size() * 2; i++) {
		indices.emplace_back(i);
	}

	// Sort indices based on the values in the stats nodes
	std::sort(indices.begin(), indices.end(), [&columns, &column_stats_nodes](const idx_t &lhs, const idx_t &rhs) {
		return PartitionedExecutionCompareStatsNodes(columns, column_stats_nodes, lhs, rhs);
	});

	// Tuned constants for finding reasonable partitions
	const auto min_row_groups_per_partition =
	    num_threads * PartitionedExecutionConfig::MIN_ROW_GROUPS_PER_THREAD_PER_PARTITION;
	const auto min_partition_count = min_row_groups_per_partition * DEFAULT_ROW_GROUP_SIZE;

	// Compute ranges
	vector<PartitionedExecutionRange> ranges;
	int64_t current_overlap = 0;
	for (idx_t i = 0; i < indices.size();) {
		const auto partition_start_i = i;
		idx_t partition_count = 0;
		double partition_overlap_ratio = NumericLimits<double>::Maximum();

		// Grow partition until "partition_count" is at least "min_partition_count"
		PartitionedExecutionGrowPartition(
		    column_stats_nodes[0], indices, i, current_overlap, partition_count, partition_overlap_ratio,
		    [&partition_count, &min_partition_count] { return partition_count >= min_partition_count; });

		// Grow partition until "partition_overlap" is less than the allowed maximum
		PartitionedExecutionGrowPartition(column_stats_nodes[0], indices, i, current_overlap, partition_count,
		                                  partition_overlap_ratio, [&partition_overlap_ratio] {
			                                  return partition_overlap_ratio <
			                                         PartitionedExecutionConfig::MAX_OVERLAP_RATIO;
		                                  });

		if (partition_overlap_ratio != 0) {
			// Temp variables for looking ahead
			auto temp_current_overlap = current_overlap;
			auto temp_i = i;
			auto temp_partition_count = partition_count;
			auto temp_partition_overlap_ratio = partition_overlap_ratio;

			// To keep track of the lowest overlap
			auto lowest_overlap_ratio = partition_overlap_ratio;
			auto lowest_i = i;

			// Look ahead for at most "min_partition_count" to see if there's a point with less overlap
			const auto end = partition_count + min_partition_count;
			PartitionedExecutionGrowPartition(
			    column_stats_nodes[0], indices, temp_i, temp_current_overlap, temp_partition_count,
			    temp_partition_overlap_ratio, [&temp_partition_count, &end] { return temp_partition_count >= end; },
			    [&temp_i, &temp_partition_overlap_ratio, &lowest_overlap_ratio, &lowest_i] {
				    if (temp_partition_overlap_ratio < lowest_overlap_ratio) {
					    lowest_overlap_ratio = temp_partition_overlap_ratio;
					    lowest_i = temp_i;
				    }
			    });

			// Actually grow if it's a success
			if (lowest_i != i) {
				PartitionedExecutionGrowPartition(column_stats_nodes[0], indices, i, current_overlap, partition_count,
				                                  partition_overlap_ratio, [&i, &lowest_i] { return i == lowest_i; });
			}
		}

		// Finally, grow if we would otherwise leave a remainder that is smaller than min_row_groups_per_partition
		if (indices.size() - 1 < min_row_groups_per_partition) {
			PartitionedExecutionGrowPartition(column_stats_nodes[0], indices, i, current_overlap, partition_count,
			                                  partition_overlap_ratio, [] { return false; });
		}

		// Add the new range to the ranges
		const auto partition_end_i = i;
		PartitionedExecutionRange range;
		if (partition_start_i != 0) {
			for (auto &stats_nodes : column_stats_nodes) {
				range.min.emplace_back(stats_nodes[indices[partition_start_i]].val);
			}
		}
		if (partition_end_i != indices.size()) {
			for (auto &stats_nodes : column_stats_nodes) {
				range.max.emplace_back(stats_nodes[indices[partition_end_i]].val);
			}
		}
		range.overlap = NumericCast<idx_t>(current_overlap);
		range.count = partition_count;
		ranges.emplace_back(std::move(range));
	}

	return ranges;
}

void PartitionExecutionSplitPipeline(Optimizer &optimizer, LogicalOperator &root, unique_ptr<LogicalOperator> &op,
                                     const vector<PartitionedExecutionColumn> &columns,
                                     const vector<PartitionedExecutionRange> &ranges) {
	vector<unique_ptr<LogicalOperator>> children;
	for (const auto &range : ranges) {
		LogicalOperatorDeepCopy deep_copy(optimizer.binder, nullptr);
		unique_ptr<LogicalOperator> copy;
		try {
			copy = deep_copy.DeepCopy(op);
		} catch (NotImplementedException &) {
			return; // Cannot copy this operator
		}

		// Do this again for the copied plan so we can get the new column bindings
		vector<PartitionedExecutionColumn> copy_columns;
		if (!PartitionedExecutionGetColumns(*copy, copy_columns) || copy_columns.empty()) {
			throw InternalException("PartitionExecutionSplitPipeline failed");
		}

		// Create filter operator
		unique_ptr<LogicalOperator> filter = make_uniq<LogicalFilter>();
		D_ASSERT(range.min.empty() || range.min.size() == columns.size());
		for (idx_t col_idx = 0; col_idx < range.min.size(); col_idx++) {
			const auto &val = range.min[col_idx];
			const auto &column_binding = copy_columns[columns[col_idx].original_idx].column_binding;
			filter->expressions.emplace_back(make_uniq<BoundComparisonExpression>(
			    ExpressionType::COMPARE_GREATERTHAN, make_uniq<BoundColumnRefExpression>(val.type(), column_binding),
			    make_uniq<BoundConstantExpression>(val)));
		}
		D_ASSERT(range.max.empty() || range.max.size() == columns.size());
		for (idx_t col_idx = 0; col_idx < range.max.size(); col_idx++) {
			const auto &val = range.max[col_idx];
			const auto &column_binding = copy_columns[columns[col_idx].original_idx].column_binding;
			filter->expressions.emplace_back(
			    make_uniq<BoundComparisonExpression>(ExpressionType::COMPARE_LESSTHANOREQUALTO,
			                                         make_uniq<BoundColumnRefExpression>(val.type(), column_binding),
			                                         make_uniq<BoundConstantExpression>(val)));
		}

		// Add the filter under the operator
		filter->children.emplace_back(std::move(copy->children[0]));
		copy->children[0] = std::move(filter);

		// Now push it down
		FilterPushdown filter_pushdown(optimizer, false);
		copy = filter_pushdown.Rewrite(std::move(copy));

		// Add cardinality estimates for these GETs
		reference<LogicalOperator> get = *copy;
		while (get.get().type != LogicalOperatorType::LOGICAL_GET) {
			get = *get.get().children[0];
		}
		get.get().SetEstimatedCardinality(range.count);

		// Add it the children for the union
		children.emplace_back(std::move(copy));
	}

	// Replace operator with union
	const auto old_bindings = op->GetColumnBindings();
	op = make_uniq<LogicalSetOperation>(optimizer.binder.GenerateTableIndex(), old_bindings.size(), std::move(children),
	                                    LogicalOperatorType::LOGICAL_UNION, true,
	                                    op->type != LogicalOperatorType::LOGICAL_ORDER_BY);
	const auto new_bindings = op->GetColumnBindings();

	// Fix up column bindings
	ColumnBindingReplacer replacer;
	for (idx_t col_idx = 0; col_idx < old_bindings.size(); col_idx++) {
		replacer.replacement_bindings.emplace_back(old_bindings[col_idx], new_bindings[col_idx]);
	}
	replacer.stop_operator = op.get();
	replacer.VisitOperator(root);
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

	if (get->EstimateCardinality(optimizer.context) < PartitionedExecutionConfig::MINIMUM_INPUT_CARDINALITY) {
		return; // Too small
	}

	GetPartitionStatsInput input(get->function, get->bind_data.get());
	const auto partition_stats = get->function.get_partition_stats(optimizer.context, input);
	if (partition_stats.size() <= 1) {
		return; // Can't split 0 or 1 partitions
	}

	const auto ranges = PartitionedExecutionComputeRanges(*op, columns, partition_stats, num_threads);
	if (ranges.size() < 2) {
		return; // Unable to compute useful ranges
	}

	PartitionExecutionSplitPipeline(optimizer, root, op, columns, ranges); // Success!
}

} // namespace duckdb
