#include "duckdb/execution/operator/projection/physical_projection.hpp"
#include "duckdb/parallel/thread_context.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/expression_iterator.hpp"
#include "duckdb/planner/expression/bound_reference_expression.hpp"
#include "duckdb/execution/operator/scan/physical_table_scan.hpp"
#include "duckdb/parallel/pipeline.hpp"

namespace duckdb {

class ProjectionGlobalState : public GlobalOperatorState {
public:
	unique_ptr<ExpressionExecutorGlobalState> expressions;
};

class ProjectionState : public OperatorState {
public:
	ProjectionState(ExecutionContext &context, const vector<unique_ptr<Expression>> &expressions,
	                optional_ptr<ExpressionExecutorGlobalState> global)
	    : executor(context.client) {
		executor.global_state = global;
		for (auto &expression : expressions) {
			executor.AddExpression(*expression);
		}
	}

	ExpressionExecutor executor;

public:
	void Finalize(const PhysicalOperator &op, ExecutionContext &context) override {
		context.thread.profiler.Flush(op);
	}

	bool SupportsReuse() const override {
		return !executor.global_state;
	}
};

PhysicalProjection::PhysicalProjection(PhysicalPlan &physical_plan, vector<LogicalType> types,
                                       vector<unique_ptr<Expression>> select_list, idx_t estimated_cardinality)
    : PhysicalOperator(physical_plan, PhysicalOperatorType::PROJECTION, std::move(types), estimated_cardinality),
      select_list(std::move(select_list)) {
}

void PhysicalProjection::InitializeOutputChunk(Allocator &allocator, DataChunk &chunk) const {
	// column references only reference their input vector. Only computed columns need a buffer
	vector<bool> initialize;
	bool any_computed = false;
	for (auto &expression : select_list) {
		const bool computed = expression->GetExpressionClass() != ExpressionClass::BOUND_REF;
		initialize.push_back(computed);
		any_computed = any_computed || computed;
	}
	if (!any_computed) {
		chunk.InitializeEmpty(types);
		return;
	}
	chunk.Initialize(allocator, types, initialize);
}

OperatorResultType PhysicalProjection::Execute(ExecutionContext &context, DataChunk &input, DataChunk &chunk,
                                               GlobalOperatorState &gstate, OperatorState &state_p) const {
	auto &state = state_p.Cast<ProjectionState>();
	state.executor.row_offset = context.source_row_offset;
	state.executor.Execute(input, chunk);
	return OperatorResultType::NEED_MORE_INPUT;
}

unique_ptr<OperatorState> PhysicalProjection::GetOperatorState(ExecutionContext &context) const {
	optional_ptr<ExpressionExecutorGlobalState> shared = op_state->Cast<ProjectionGlobalState>().expressions.get();
	if (shared && (!context.pipeline || context.pipeline->IsExternalInput() ||
	               context.pipeline->GetSource().get() != &children[0].get())) {
		shared = nullptr;
	}
	return make_uniq<ProjectionState>(context, select_list, shared);
}

//! Prototype: a source ordinal is an RNG draw ordinal only for unconditional scalar projections.
static bool SupportsDensePositions(const Expression &expr) {
	if (!expr.IsVolatile()) {
		return true;
	}
	if (expr.GetExpressionClass() != ExpressionClass::BOUND_FUNCTION) {
		return false;
	}
	auto &function = expr.Cast<BoundFunctionExpression>().Function();
	if (function.HasBindLambdaCallback() ||
	    (function.GetStability() == FunctionStability::VOLATILE && !function.HasInitGlobalStateCallback())) {
		return false;
	}
	for (auto &child : expr.Cast<BoundFunctionExpression>().GetChildren()) {
		if (!SupportsDensePositions(*child)) {
			return false;
		}
	}
	return true;
}

unique_ptr<GlobalOperatorState> PhysicalProjection::GetGlobalOperatorState(ClientContext &) const {
	auto result = make_uniq<ProjectionGlobalState>();
	// Keep the prototype's row-position contract local: no intervening filter, join, or expansion.
	if (children.size() != 1 || children[0].get().type != PhysicalOperatorType::TABLE_SCAN ||
	    !children[0].get().Cast<PhysicalTableScan>().function.provides_row_offsets) {
		return std::move(result);
	}
	bool has_volatile = false;
	for (auto &expression : select_list) {
		if (!SupportsDensePositions(*expression)) {
			return std::move(result);
		}
		has_volatile = has_volatile || expression->IsVolatile();
	}
	if (has_volatile) {
		result->expressions = make_uniq<ExpressionExecutorGlobalState>();
	}
	return std::move(result);
}

bool PhysicalProjection::ParallelOperator() const {
	for (auto &expression : select_list) {
		bool requires_ordered_execution = false;
		ExpressionIterator::VisitExpression<BoundFunctionExpression>(*expression, [&](const auto &function) {
			if (function.RequiresOrderedExecution()) {
				requires_ordered_execution = true;
			}
		});
		if (requires_ordered_execution) {
			return false;
		}
	}
	return true;
}

InsertionOrderPreservingMap<string> PhysicalProjection::ParamsToString() const {
	InsertionOrderPreservingMap<string> result;
	string projections;
	for (idx_t i = 0; i < select_list.size(); i++) {
		if (i > 0) {
			projections += "\n";
		}
		auto &expr = select_list[i];
		projections += expr->GetName();
	}
	result["__projections__"] = projections;
	SetEstimatedCardinality(result, estimated_cardinality);
	return result;
}

} // namespace duckdb
