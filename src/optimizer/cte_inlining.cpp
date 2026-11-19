#include "duckdb/optimizer/cte_inlining.hpp"

#include "duckdb/optimizer/column_binding_replacer.hpp"
#include "duckdb/optimizer/filter_pushdown.hpp"
#include "duckdb/planner/expression/bound_conjunction_expression.hpp"
#include "duckdb/planner/expression/list.hpp"
#include "duckdb/planner/operator/logical_cteref.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"
#include "duckdb/planner/operator/logical_materialized_cte.hpp"

#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/planner/logical_operator_deep_copy.hpp"
#include "duckdb/planner/operator/logical_prepare.hpp"

#include "duckdb/function/scalar/generic_functions.hpp"

namespace duckdb {

CTEInlining::CTEInlining(Optimizer &optimizer_p) : optimizer(optimizer_p) {
}

unique_ptr<LogicalOperator> CTEInlining::Optimize(unique_ptr<LogicalOperator> op) {
	TryInlining(op);
	return op;
}

static void PullUpScanFilters(unique_ptr<LogicalOperator> &op) {
	for (auto &child : op->children) {
		PullUpScanFilters(child);
	}
	if (op->type != LogicalOperatorType::LOGICAL_GET) {
		return;
	}
	auto &get = op->Cast<LogicalGet>();
	if (!get.table_filters.HasFilters() || get.table_filters.HasMultiColumnFilters()) {
		return;
	}
	vector<unique_ptr<Expression>> expressions;
	for (auto &entry : get.table_filters) {
		auto index = entry.GetIndex();
		BoundColumnRefExpression column(get.GetColumnType(get.GetColumnIds()[index]),
		                                ColumnBinding(get.table_index, index));
		expressions.push_back(entry.Filter().ToExpression(column));
	}
	get.table_filters.ClearFilters();
	auto filter = make_uniq<LogicalFilter>();
	filter->expressions = std::move(expressions);
	filter->projection_map = std::move(get.projection_ids);
	filter->children.push_back(std::move(op));
	op = std::move(filter);
}

static void EstimateInlinedPlan(LogicalOperator &op, ClientContext &context) {
	for (auto &child : op.children) {
		EstimateInlinedPlan(*child, context);
	}
	op.SetEstimatedCardinality(op.EstimateCardinality(context));
}

static bool CanStreamCTEReferences(LogicalOperator &op, TableIndex cte_index, bool streaming_sink = false) {
	if (op.type == LogicalOperatorType::LOGICAL_RECURSIVE_CTE || op.type == LogicalOperatorType::LOGICAL_DELIM_JOIN ||
	    op.type == LogicalOperatorType::LOGICAL_EXTENSION_OPERATOR) {
		return false;
	}
	if (op.type == LogicalOperatorType::LOGICAL_CTE_REF) {
		return op.Cast<LogicalCTERef>().cte_index != cte_index || streaming_sink;
	}
	switch (op.type) {
	case LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY: {
		// Aggregate consumers have no build-side dependencies and do not stop their input early.
		streaming_sink = true;
		auto &aggregate = op.Cast<LogicalAggregate>();
		if (aggregate.groups.empty() && aggregate.grouping_sets.size() <= 1) {
			for (auto &expr : op.expressions) {
				if (expr->Cast<BoundAggregateExpression>().Function().GetOrderDependent() ==
				    AggregateOrderDependent::ORDER_DEPENDENT) {
					streaming_sink = false;
					break;
				}
			}
		}
		break;
	}
	case LogicalOperatorType::LOGICAL_PROJECTION:
	case LogicalOperatorType::LOGICAL_FILTER:
		break;
	default:
		// Other pipeline boundaries need their own proof of direct-input support and acyclic dependencies.
		streaming_sink = false;
		break;
	}
	if (streaming_sink) {
		LogicalOperatorVisitor::EnumerateExpressions(op, [&](unique_ptr<Expression> *expr) {
			ExpressionIterator::VisitExpression<BoundFunctionExpression>(
			    **expr, [&](const auto &function) { streaming_sink &= !function.RequiresOrderedExecution(); });
		});
	}
	for (auto &child : op.children) {
		if (!CanStreamCTEReferences(*child, cte_index, streaming_sink)) {
			return false;
		}
	}
	return true;
}

static idx_t CountBaseTableReferences(const LogicalOperator &op) {
	idx_t number_of_references = 0;
	if (op.type == LogicalOperatorType::LOGICAL_GET) {
		number_of_references++;
	}
	for (auto &child : op.children) {
		number_of_references += CountBaseTableReferences(*child);
	}

	return number_of_references;
}

static idx_t CountCTEReferences(const LogicalOperator &op, TableIndex cte_index) {
	if (op.type == LogicalOperatorType::LOGICAL_CTE_REF) {
		auto &cte = op.Cast<LogicalCTERef>();
		if (cte.cte_index == cte_index) {
			return 1;
		}
	}
	idx_t number_of_references = 0;
	for (auto &child : op.children) {
		number_of_references += CountCTEReferences(*child, cte_index);
	}

	return number_of_references;
}

bool CTEInlining::FinalizeStreamingCTEs(unique_ptr<LogicalOperator> &op, bool recursive) {
	if (op->type == LogicalOperatorType::LOGICAL_PREPARE) {
		parameter_data = op->Cast<LogicalPrepare>().prepared->value_map;
	}
	recursive |= op->type == LogicalOperatorType::LOGICAL_RECURSIVE_CTE;
	bool changed = false;
	for (auto &child : op->children) {
		changed |= FinalizeStreamingCTEs(child, recursive);
	}
	if (op->type != LogicalOperatorType::LOGICAL_MATERIALIZED_CTE) {
		return changed;
	}
	auto &cte = op->Cast<LogicalMaterializedCTE>();
	if (cte.materialize != CTEMaterialize::CTE_MATERIALIZE_STREAMING) {
		return changed;
	}
	const auto ref_count = CountCTEReferences(*cte.children[1], cte.table_index);
	if (ref_count > 1 && !recursive && CanStreamCTEReferences(*cte.children[1], cte.table_index)) {
		return changed;
	}
	// Recombine the shared predicates with each consumer's filters after inlining.
	PullUpScanFilters(op->children[0]);
	if (!Inline(op->children[1], *op, ref_count > 1)) {
		throw InternalException("Cannot inline a streaming CTE");
	}
	op = std::move(op->children[1]);
	return true;
}

static bool ContainsLimit(const LogicalOperator &op) {
	if (op.type == LogicalOperatorType::LOGICAL_LIMIT || op.type == LogicalOperatorType::LOGICAL_TOP_N) {
		return true;
	}
	if (op.children.size() != 1) {
		return false;
	}
	for (auto &child : op.children) {
		if (ContainsLimit(*child)) {
			return true;
		}
	}
	return false;
}

bool CTEInlining::EndsInAggregateOrDistinct(const LogicalOperator &op) {
	switch (op.type) {
	case LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY:
	case LogicalOperatorType::LOGICAL_DISTINCT:
	case LogicalOperatorType::LOGICAL_WINDOW:
		return true;
	default:
		break;
	}
	if (op.children.size() != 1) {
		return false;
	}
	for (auto &child : op.children) {
		if (EndsInAggregateOrDistinct(*child)) {
			return true;
		}
	}
	return false;
}

static bool EndsInDummyScan(const LogicalOperator &op) {
	if (op.type == LogicalOperatorType::LOGICAL_DUMMY_SCAN || op.type == LogicalOperatorType::LOGICAL_EMPTY_RESULT ||
	    op.type == LogicalOperatorType::LOGICAL_CTE_REF) {
		return true;
	}
	if (op.children.size() != 1) {
		return false;
	}
	for (auto &child : op.children) {
		if (EndsInDummyScan(*child)) {
			return true;
		}
	}
	return false;
}

static bool ContainsDelimGet(const LogicalOperator &op) {
	if (op.type == LogicalOperatorType::LOGICAL_DELIM_GET) {
		return true;
	}
	for (auto &child : op.children) {
		if (ContainsDelimGet(*child)) {
			return true;
		}
	}
	return false;
}

static bool HasCTEReferenceBelowDelimJoin(const LogicalOperator &op, TableIndex cte_index,
                                          bool below_delim_join = false) {
	if (op.type == LogicalOperatorType::LOGICAL_CTE_REF) {
		auto &cteref = op.Cast<LogicalCTERef>();
		if (cteref.cte_index == cte_index) {
			return below_delim_join;
		}
	}
	auto child_below_delim_join = below_delim_join || op.type == LogicalOperatorType::LOGICAL_DELIM_JOIN;
	for (auto &child : op.children) {
		if (HasCTEReferenceBelowDelimJoin(*child, cte_index, child_below_delim_join)) {
			return true;
		}
	}
	return false;
}

void CTEInlining::TryInlining(unique_ptr<LogicalOperator> &op) {
	if (op->type == LogicalOperatorType::LOGICAL_PREPARE) {
		// we are in a prepare statement, if we have to copy an operator during inlining,
		// we have to be careful to use the correct parameter data
		auto &prepare = op->Cast<LogicalPrepare>();
		parameter_data = prepare.prepared->value_map;
	}

	// traverse children first, so we can inline the deepest CTEs first
	for (auto &child : op->children) {
		TryInlining(child);
	}

	if (op->type == LogicalOperatorType::LOGICAL_MATERIALIZED_CTE) {
		auto &cte = op->Cast<LogicalMaterializedCTE>();
		auto ref_count = CountCTEReferences(*op, cte.table_index);
		if (ref_count == 0) {
			if (cte.children[0]->HasSideEffects()) {
				// Side-effecting CTEs must always execute even when unreferenced
				return;
			}
			// this CTE is not referenced, we can remove it
			op = std::move(op->children[1]);
			return;
		}
		if (cte.children[0]->HasSideEffects()) {
			// Never inline a side-effecting CTE: the LOGICAL_MATERIALIZED_CTE guarantees
			// that it executes exactly once and before the query side.
			return;
		}
		if (ContainsDelimGet(*cte.children[0]) && HasCTEReferenceBelowDelimJoin(*op->children[1], cte.table_index)) {
			// Inlining a CTE that already contains a DELIM_GET stays safe while all matching CTE scans remain outside
			// DELIM_JOIN subtrees, but once a scan is nested below another DELIM_JOIN the inlined DELIM_GETs can attach
			// to the wrong duplicate-elimination source.
			return;
		}
		if (cte.materialize == CTEMaterialize::CTE_MATERIALIZE_ALWAYS) {
			// This CTE is always materialized, we cannot inline it
			return;
		}
		if (ref_count == 1) {
			// this CTE is only referenced once, we can inline it directly without copying
			bool success = Inline(op->children[1], *op, false);
			if (success) {
				op = std::move(op->children[1]);
			}
			return;
		}
		if (ref_count > 1) {
			if (cte.materialize == CTEMaterialize::CTE_MATERIALIZE_NEVER) {
				// this CTE is referenced multiple times, but we are not allowed to materialize it
				// we have to inline it if possible
				bool success = Inline(op->children[1], *op, true);
				if (success) {
					op = std::move(op->children[1]);
				}
				return;
			}
			// check if we can inline this CTE
			PreventInlining prevent_inlining;
			prevent_inlining.VisitOperator(*op->children[0]);

			if (prevent_inlining.prevent_inlining) {
				// we cannot inline this CTE, we have to keep it materialized
				return;
			}

			// Prevent inlining if the CTE ends in an aggregate or distinct operator
			// This mimics the behavior of the CTE materialization in the binder
			if (EndsInAggregateOrDistinct(*op->children[0])) {
				return;
			}

			bool is_cheap_to_inline = op->children[0]->type == LogicalOperatorType::LOGICAL_EMPTY_RESULT ||
			                          op->children[0]->type == LogicalOperatorType::LOGICAL_CTE_REF ||
			                          EndsInDummyScan(*op->children[0]);

			// Check how many base table references the CTE has
			auto base_table_references = CountBaseTableReferences(*op->children[0]);

			if (!is_cheap_to_inline && base_table_references > 2 && base_table_references * ref_count > 10) {
				return;
			}

			// CTEs require full materialization before the CTE scans begin,
			// LIMIT and TOP_N operators cannot abort the materialization,
			// even if only a part of the CTE result is needed.
			// Therefore, we check if the CTE Scans are below the LIMIT or TOP_N operator
			// and if so, we try to inline the CTE definition.
			if (is_cheap_to_inline || ContainsLimit(*op->children[1])) {
				// this CTE is referenced multiple times and has a limit, we want to inline it
				bool success = Inline(op->children[1], *op, true);
				if (success) {
					op = std::move(op->children[1]);
				}
				return;
			}
		}
	}
}

bool CTEInlining::Inline(unique_ptr<LogicalOperator> &op, LogicalOperator &materialized_cte, bool requires_copy) {
	if (op->type == LogicalOperatorType::LOGICAL_CTE_REF) {
		auto &cteref = op->Cast<LogicalCTERef>();
		auto &cte = materialized_cte.Cast<LogicalCTE>();
		if (cteref.cte_index == cte.table_index) {
			unique_ptr<LogicalOperator> &definition = cte.children[0];
			unique_ptr<LogicalOperator> copy;
			if (requires_copy) {
				// there are multiple references to the CTE, we need to copy it
				LogicalOperatorDeepCopy deep_copy(optimizer.binder, parameter_data);
				try {
					copy = deep_copy.DeepCopy(definition);
				} catch (NotImplementedException &ex) {
					// if we have to copy the lhs of a CTE, but we cannot copy the operator, we have to
					// stop inlining and keep the materialized CTE instead
					return false;
				}
			}
			op = LogicalProjection::CreateIdentity(cteref.table_index,
			                                       requires_copy ? std::move(copy) : std::move(definition));
			if (materialized_cte.Cast<LogicalMaterializedCTE>().materialize ==
			    CTEMaterialize::CTE_MATERIALIZE_STREAMING) {
				// Late inlining must restore the estimates that serialization deliberately omits.
				EstimateInlinedPlan(*op, optimizer.context);
			}
			return true;
		}
		return true;
	} else {
		bool success = true;
		for (auto &child : op->children) {
			success &= Inline(child, materialized_cte, requires_copy);
		}
		return success;
	}
}

void PreventInlining::VisitOperator(LogicalOperator &op) {
	VisitOperatorExpressions(op);
	// We can stop checking early if we already know that inlining is not possible
	if (!prevent_inlining) {
		VisitOperatorChildren(op);
	}
}

void PreventInlining::VisitExpression(unique_ptr<Expression> *expression) {
	auto &expr = *expression;

	if (expr->GetExpressionClass() == ExpressionClass::BOUND_FUNCTION) {
		auto &bound_function = expr->Cast<BoundFunctionExpression>();
		// if we encounter the ErrorFun function, we still want to inline
		if (bound_function.Function().GetName() == "error") {
			return;
		}

		if (expr->IsVolatile()) {
			prevent_inlining = true;
			return;
		}
	}
	VisitExpressionChildren(**expression);
}

} // namespace duckdb
