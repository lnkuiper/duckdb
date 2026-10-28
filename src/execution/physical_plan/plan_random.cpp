#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/expression_iterator.hpp"
#include "duckdb/planner/logical_operator_visitor.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"

namespace duckdb {
namespace {

struct OrdinalColumn {
	ColumnBinding binding;
	ProjectionIndex index;
};

class RandomPlanner {
public:
	void Visit(LogicalOperator &op) {
		// Multiple or correlated input streams do not share a source ordinal.
		if (op.children.size() != 1 || op.type == LogicalOperatorType::LOGICAL_GET) {
			return;
		}
		Visit(*op.children[0]);
		// Aggregate inputs are extracted into a projection during physical planning.
		if ((op.type != LogicalOperatorType::LOGICAL_PROJECTION &&
		     op.type != LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY &&
		     op.type != LogicalOperatorType::LOGICAL_FILTER) ||
		    !SupportsOrdinality(*op.children[0])) {
			return;
		}
		bool needs_ordinality = false;
		LogicalOperatorVisitor::EnumerateExpressions(
		    op, [&](unique_ptr<Expression> *expr) { needs_ordinality |= HasRandom(**expr); });
		if (!needs_ordinality) {
			return;
		}
		auto column_count = op.GetColumnBindings().size();
		auto column = AddOrdinality(*op.children[0]);
		LogicalOperatorVisitor::EnumerateExpressions(
		    op, [&](unique_ptr<Expression> *expr) { BindRandom(**expr, column.binding); });
		if (op.type == LogicalOperatorType::LOGICAL_FILTER && op.GetColumnBindings().size() != column_count) {
			// Keep the helper column out of the filter's output.
			auto &projection_map = op.Cast<LogicalFilter>().projection_map;
			D_ASSERT(projection_map.empty());
			for (idx_t i = 0; i < column_count; i++) {
				projection_map.emplace_back(i);
			}
		}
	}

private:
	static bool SupportsOrdinality(LogicalOperator &op) {
		if (op.type == LogicalOperatorType::LOGICAL_GET) {
			return op.Cast<LogicalGet>().SupportsRowNumber();
		}
		if (op.type == LogicalOperatorType::LOGICAL_WINDOW && op.expressions.empty()) {
			return SupportsOrdinality(*op.children[0]);
		}
		if (op.type != LogicalOperatorType::LOGICAL_PROJECTION && op.type != LogicalOperatorType::LOGICAL_FILTER) {
			return false;
		}
		return op.children.size() == 1 && SupportsOrdinality(*op.children[0]);
	}

	static bool IsRandom(const Expression &expr) {
		if (expr.GetExpressionClass() != ExpressionClass::BOUND_FUNCTION) {
			return false;
		}
		auto &function = expr.Cast<BoundFunctionExpression>().Function();
		return function.GetName() == "random" && function.GetCatalogName() == SYSTEM_CATALOG &&
		       function.GetSchemaName() == Identifier::DefaultSchema() && function.GetArguments().empty();
	}

	static bool HasRandom(const Expression &expr) {
		if (IsRandom(expr)) {
			return true;
		}
		bool result = false;
		ExpressionIterator::EnumerateChildren(expr, [&](const Expression &child) { result |= HasRandom(child); });
		return result;
	}

	void BindRandom(Expression &expr, ColumnBinding binding) {
		ExpressionIterator::EnumerateChildren(expr, [&](Expression &child) { BindRandom(child, binding); });
		if (!IsRandom(expr)) {
			return;
		}
		auto &function = expr.Cast<BoundFunctionExpression>();
		function.FunctionMutable().GetArguments() = {LogicalType::BIGINT, LogicalType::UBIGINT};
		function.GetChildrenMutable().push_back(make_uniq<BoundColumnRefExpression>(LogicalType::BIGINT, binding));
		function.GetChildrenMutable().push_back(make_uniq<BoundConstantExpression>(Value::UBIGINT(next_stream++)));
	}

	OrdinalColumn AddOrdinality(LogicalOperator &op) {
		auto entry = columns.find(op);
		if (entry != columns.end()) {
			return entry->second;
		}
		auto column = AddColumn(op);
		columns.emplace(op, column);
		return column;
	}

	OrdinalColumn AddColumn(LogicalOperator &op) {
		if (op.type == LogicalOperatorType::LOGICAL_GET) {
			auto &get = op.Cast<LogicalGet>();
			auto index = get.TryGetProjectionIndex(COLUMN_IDENTIFIER_ROW_NUMBER);
			if (!index.IsValid()) {
				index = get.AddColumnId(COLUMN_IDENTIFIER_ROW_NUMBER);
			}
			auto output_index = index;
			if (!get.projection_ids.empty()) {
				auto position = std::find(get.projection_ids.begin(), get.projection_ids.end(), index);
				if (position == get.projection_ids.end()) {
					output_index = ProjectionIndex(get.projection_ids.size());
					get.projection_ids.push_back(index);
				} else {
					output_index = ProjectionIndex(NumericCast<idx_t>(position - get.projection_ids.begin()));
				}
			}
			return {{get.table_index, index}, output_index};
		}
		auto column = AddOrdinality(*op.children[0]);
		if (op.type == LogicalOperatorType::LOGICAL_WINDOW) {
			D_ASSERT(op.expressions.empty());
			return column;
		}
		if (op.type == LogicalOperatorType::LOGICAL_PROJECTION) {
			auto &projection = op.Cast<LogicalProjection>();
			auto index = ProjectionIndex(projection.expressions.size());
			projection.expressions.push_back(make_uniq<BoundColumnRefExpression>(LogicalType::BIGINT, column.binding));
			return {{projection.table_index, index}, index};
		}
		auto &filter = op.Cast<LogicalFilter>();
		if (!filter.projection_map.empty()) {
			auto output_index = ProjectionIndex(filter.projection_map.size());
			filter.projection_map.push_back(column.index);
			column.index = output_index;
		}
		return column;
	}

	reference_map_t<LogicalOperator, OrdinalColumn> columns;
	uint64_t next_stream = 0;
};

} // namespace

void PhysicalPlanGenerator::PlanRandom(LogicalOperator &op) {
	RandomPlanner planner;
	planner.Visit(op);
}

} // namespace duckdb
