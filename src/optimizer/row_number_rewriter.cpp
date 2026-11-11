#include "duckdb/optimizer/row_number_rewriter.hpp"

#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/expression/bound_window_expression.hpp"
#include "duckdb/optimizer/column_binding_replacer.hpp"
#include "duckdb/planner/operator/logical_window.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/expression_iterator.hpp"
#include "duckdb/planner/logical_operator_visitor.hpp"
#include "duckdb/planner/binder.hpp"

namespace duckdb {

namespace {

struct RandomColumn {
	ColumnBinding binding;
	ProjectionIndex index;
};

class RandomColumnPushdown {
public:
	void Visit(LogicalOperator &op) {
		if (op.children.size() != 1 || op.type == LogicalOperatorType::LOGICAL_GET) {
			return;
		}
		Visit(*op.children[0]);
		if ((op.type != LogicalOperatorType::LOGICAL_PROJECTION &&
		     op.type != LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY &&
		     op.type != LogicalOperatorType::LOGICAL_FILTER) ||
		    !SupportsPushdown(*op.children[0])) {
			return;
		}
		auto column_count = op.GetColumnBindings().size();
		LogicalOperatorVisitor::EnumerateExpressions(
		    op, [&](unique_ptr<Expression> *expr) { PushRandom(*expr, *op.children[0]); });
		if (op.type == LogicalOperatorType::LOGICAL_FILTER && op.GetColumnBindings().size() != column_count) {
			auto &projection_map = op.Cast<LogicalFilter>().projection_map;
			D_ASSERT(projection_map.empty());
			for (idx_t i = 0; i < column_count; i++) {
				projection_map.emplace_back(i);
			}
		}
	}

private:
	static bool SupportsPushdown(LogicalOperator &op) {
		if (op.type == LogicalOperatorType::LOGICAL_GET) {
			return op.Cast<LogicalGet>().SupportsRowNumber();
		}
		if (op.type != LogicalOperatorType::LOGICAL_PROJECTION && op.type != LogicalOperatorType::LOGICAL_FILTER) {
			return false;
		}
		return op.children.size() == 1 && SupportsPushdown(*op.children[0]);
	}

	static void PushRandom(unique_ptr<Expression> &expr, LogicalOperator &input) {
		ExpressionIterator::EnumerateChildren(*expr, [&](unique_ptr<Expression> &child) { PushRandom(child, input); });
		if (expr->GetExpressionClass() != ExpressionClass::BOUND_FUNCTION) {
			return;
		}
		auto &function = expr->Cast<BoundFunctionExpression>();
		if (function.Function().GetName() != "random" || function.Function().GetCatalogName() != SYSTEM_CATALOG ||
		    function.Function().GetSchemaName() != Identifier::DefaultSchema() || !function.GetChildren().empty() ||
		    expr->GetReturnType() != LogicalType::DOUBLE) {
			return;
		}
		auto alias = expr->GetAlias();
		auto column = AddColumn(input, std::move(expr));
		expr = make_uniq<BoundColumnRefExpression>(alias, LogicalType::DOUBLE, column.binding);
	}

	static RandomColumn AddColumn(LogicalOperator &op, unique_ptr<Expression> expr) {
		if (op.type == LogicalOperatorType::LOGICAL_GET) {
			auto &get = op.Cast<LogicalGet>();
			if (get.GetColumnIds().empty()) {
				get.AddColumnId(get.GetAnyColumn());
			}
			auto output_index = ProjectionIndex(get.GetColumnBindings().size());
			auto binding_index = ProjectionIndex(get.GetColumnIds().size() + get.expressions.size());
			get.expressions.push_back(std::move(expr));
			return {{get.table_index, binding_index}, output_index};
		}
		auto column = AddColumn(*op.children[0], std::move(expr));
		if (op.type == LogicalOperatorType::LOGICAL_PROJECTION) {
			auto &projection = op.Cast<LogicalProjection>();
			auto index = ProjectionIndex(projection.expressions.size());
			projection.expressions.push_back(make_uniq<BoundColumnRefExpression>(LogicalType::DOUBLE, column.binding));
			return {{projection.table_index, index}, index};
		}
		auto &filter = op.Cast<LogicalFilter>();
		if (!filter.projection_map.empty()) {
			auto index = ProjectionIndex(filter.projection_map.size());
			filter.projection_map.push_back(column.index);
			column.index = index;
		}
		return column;
	}
};

void RestoreScanExpressionsInternal(unique_ptr<LogicalOperator> &op, unique_ptr<LogicalOperator> &root,
                                    Binder &binder) {
	for (auto &child : op->children) {
		RestoreScanExpressionsInternal(child, root, binder);
	}
	if (op->type != LogicalOperatorType::LOGICAL_GET || op->expressions.empty()) {
		return;
	}
	auto old_bindings = op->GetColumnBindings();
	auto computed = std::move(op->expressions);
	auto projection = LogicalProjection::CreateIdentity(binder.GenerateTableIndex(), std::move(op));
	for (auto &expression : computed) {
		projection->expressions.push_back(std::move(expression));
	}
	CorrelatedColumnBindingReplacer replacer;
	replacer.AddReplacements(old_bindings, projection->GetColumnBindings());
	replacer.stop_operator = projection.get();
	op = std::move(projection);
	op->ResolveOperatorTypes();
	replacer.VisitOperator(*root);
}

} // namespace

void RowNumberRewriter::RestoreScanExpressions(unique_ptr<LogicalOperator> &plan, Binder &binder) {
	RestoreScanExpressionsInternal(plan, plan, binder);
}

bool RowNumberRewriter::CanOptimize(LogicalOperator &op) {
	if (op.type != LogicalOperatorType::LOGICAL_WINDOW) {
		return false;
	}
	if (op.expressions.size() != 1) {
		return false;
	}

	auto &expression = op.expressions[0];
	if (expression->GetExpressionType() != ExpressionType::WINDOW_ROW_NUMBER) {
		return false;
	}
	auto &window_expr = expression->Cast<BoundWindowExpression>();
	if (!window_expr.Partitions().empty() || !window_expr.OrderBy().empty()) {
		return false;
	}

	// Should be followed by a get
	auto &window_ch = op.children[0];
	if (window_ch->type != LogicalOperatorType::LOGICAL_GET) {
		return false;
	}

	auto &get = window_ch->Cast<LogicalGet>();

	return get.expressions.empty() && get.SupportsRowNumber();
}

unique_ptr<LogicalOperator> RowNumberRewriter::Optimize(unique_ptr<LogicalOperator> op) {
	ColumnBindingReplacer replacer;
	LogicalOperator *root = op.get();
	op = RewritePlan(std::move(op), replacer);

	if (!replacer.replacement_bindings.empty()) {
		replacer.VisitOperator(*root);
	}

	// Computed bindings follow the source columns, so finish adding ROW_NUMBER columns first.
	RandomColumnPushdown random_pushdown;
	random_pushdown.Visit(*op);
	return op;
}

unique_ptr<LogicalOperator> RowNumberRewriter::RewriteGet(unique_ptr<LogicalOperator> op,
                                                          ColumnBindingReplacer &replacer) {
	auto &window = op->Cast<LogicalWindow>();
	auto &child = window.children[0];
	auto &get = child->Cast<LogicalGet>();

	// Add virtual row_number column to the LogicalGet
	auto proj_index = get.AddColumnId(COLUMN_IDENTIFIER_ROW_NUMBER);
	if (!get.projection_ids.empty()) {
		get.projection_ids.push_back(proj_index);
	}
	get.types.push_back(LogicalType::BIGINT);

	const auto child_bindings = get.GetColumnBindings();

	// Get the old window bindings before we remove the window operator
	const auto old_window_bindings = window.GetColumnBindings();

	// Map old window bindings to new get bindings
	for (idx_t i = 0; i < old_window_bindings.size(); i++) {
		ColumnBinding target_binding;
		if (i < child_bindings.size() - 1) {
			// Map existing columns through
			target_binding = child_bindings[i];
		} else {
			// Map the window expression to the virtual ROW_NUMBER column
			target_binding = child_bindings.back();
		}
		replacer.replacement_bindings.emplace_back(old_window_bindings[i], target_binding);
	}

	replacer.stop_operator = child.get();
	return std::move(window.children[0]);
}

unique_ptr<LogicalOperator> RowNumberRewriter::RewritePlan(unique_ptr<LogicalOperator> op,
                                                           ColumnBindingReplacer &replacer) {
	if (CanOptimize(*op)) {
		return RewriteGet(std::move(op), replacer);
	}

	// Recurse into children
	for (auto &child : op->children) {
		child = RewritePlan(std::move(child), replacer);
	}
	return op;
}

} // namespace duckdb
