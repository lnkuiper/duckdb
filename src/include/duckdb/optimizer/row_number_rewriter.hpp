//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/optimizer/row_number_rewriter.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once
#include "duckdb/common/unique_ptr.hpp"
#include "duckdb/planner/logical_operator.hpp"
#include "duckdb/optimizer/column_binding_replacer.hpp"

namespace duckdb {

class Optimizer;
class Binder;

class RowNumberRewriter {
public:
	unique_ptr<LogicalOperator> Optimize(unique_ptr<LogicalOperator> op);
	//! Restore ordinary expressions before an already optimized plan re-enters the optimizer.
	static void RestoreScanExpressions(unique_ptr<LogicalOperator> &plan, Binder &binder);

private:
	unique_ptr<LogicalOperator> RewritePlan(unique_ptr<LogicalOperator> op, ColumnBindingReplacer &replacer);
	bool CanOptimize(LogicalOperator &op);
	unique_ptr<LogicalOperator> RewriteGet(unique_ptr<LogicalOperator> op, ColumnBindingReplacer &replacer);
};

} // namespace duckdb
