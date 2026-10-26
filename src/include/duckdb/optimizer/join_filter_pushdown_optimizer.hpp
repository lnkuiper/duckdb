//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/optimizer/join_filter_pushdown_optimizer.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/planner/logical_operator_visitor.hpp"
#include "duckdb/planner/column_binding_map.hpp"
#include "duckdb/common/unordered_set.hpp"

namespace duckdb {
class Optimizer;
struct JoinFilterPushdownColumn;
struct PushdownFilterTarget;

enum class EqualityPropagation { DISABLED, ENABLED };
enum class CTEFilterTargetMode : uint8_t { IGNORE, INCLUDE };

//! The JoinFilterPushdownOptimizer links comparison joins to data sources to enable dynamic execution-time filter
//! pushdown
class JoinFilterPushdownOptimizer : public LogicalOperatorVisitor {
public:
	explicit JoinFilterPushdownOptimizer(Optimizer &optimizer);

public:
	void VisitOperator(LogicalOperator &op) override;
	void Optimize(LogicalOperator &op);
	//! Bind CTE filter destinations against the final logical plan, before column references are resolved.
	static void FinalizeCTEFilters(LogicalOperator &op);
	//! Find scan targets independently for each column, allowing both NULL-retaining and NULL-rejecting filters.
	//! Equality propagation can return multiple equivalent columns per scan; single-column callers leave it disabled.
	static void GetPushdownFilterTargets(LogicalOperator &op, vector<JoinFilterPushdownColumn> columns,
	                                     vector<PushdownFilterTarget> &targets,
	                                     EqualityPropagation equality_propagation = EqualityPropagation::DISABLED,
	                                     CTEFilterTargetMode cte_targets = CTEFilterTargetMode::IGNORE);

	static bool IsFiltering(const unique_ptr<LogicalOperator> &op);

private:
	void GenerateJoinFilters(LogicalComparisonJoin &join);
	static bool GetJoinFilterColumn(const LogicalComparisonJoin &join, idx_t condition_idx,
	                                JoinFilterPushdownColumn &column);
	static bool ComputeAggregatesForPerfectHash(const LogicalComparisonJoin &join, idx_t filter_count);

private:
	Optimizer &optimizer;
	unordered_set<TableIndex> eligible_ctes;
};
} // namespace duckdb
