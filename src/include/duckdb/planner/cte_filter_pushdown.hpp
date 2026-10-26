//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/planner/cte_filter_pushdown.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/planner/column_binding.hpp"
#include "duckdb/planner/table_filter_set.hpp"

namespace duckdb {

struct CTEFilterPushdownTarget {
	shared_ptr<DynamicTableFilterSet> scan_filters;
	shared_ptr<DynamicTableFilterGroup> group;
};

//! Transient consumer identity proof and filter ownership; deliberately not serialized.
struct CTEFilterPushdownInfo {
	vector<ColumnBinding> producer_bindings;
	vector<LogicalType> producer_types;
	map<TableIndex, vector<LogicalType>> consumers;
	//! Populated against the final logical plan during physical planning.
	vector<CTEFilterPushdownTarget> targets;
};

} // namespace duckdb
