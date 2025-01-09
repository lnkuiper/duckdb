#include "duckdb/planner/binder.hpp"
#include "duckdb/planner/tableref/bound_subqueryref.hpp"
#include "duckdb/planner/expression/bound_reference_expression.hpp"

namespace duckdb {

unique_ptr<LogicalOperator> Binder::CreatePlan(BoundSubqueryRef &ref) {
	// generate the logical plan for the subquery
	// this happens separately from the current LogicalPlan generation
	ref.binder->is_outside_flattened = is_outside_flattened;
	auto subquery = ref.binder->CreatePlan(*ref.subquery);
	if (ref.binder->has_unplanned_dependent_joins) {
		has_unplanned_dependent_joins = true;
	} // else if (subquery->type == LogicalOperatorType::LOGICAL_PROJECTION) {
	  // if (subquery->GetColumnBindings().size() == subquery->children[0]->GetColumnBindings().size()) {
	  // 	auto &proj = subquery->Cast<LogicalProjection>();
	  // 	auto &child = subquery->children[0];
	  // 	// check if this projection can be omitted entirely
	  // 	// this happens if a projection simply emits the columns in the same order
	  // 	// e.g. PROJECTION(#0, #1, #2, #3, ...)
	  // 	bool omit_projection = true;
	  // 	auto subquery_bindings = subquery->GetColumnBindings();
	  // 	auto child_bindings = subquery->children[0]->GetColumnBindings();
	  // 	if (subquery_bindings.size() == child_bindings.size()) {
	  // 		for (idx_t i = 0; i < subquery_bindings.size(); i++) {
	  // 			if (proj.expressions[i]->GetExpressionType() != ExpressionType::BOUND_COLUMN_REF) {
	  // 				omit_projection = false;
	  // 				break;
	  // 			}
	  // 			auto &bound_ref = proj.expressions[i]->Cast<BoundColumnRefExpression>();
	  // 			if (bound_ref.binding != child_bindings[i]) {
	  // 				omit_projection = false;
	  // 				break;
	  // 			}
	  // 		}
	  // 	}
	  // 	if (omit_projection) {
	  // 		return std::move(subquery->children[0]);
	  // 	}
	  // }
	//}
	return subquery;
}

} // namespace duckdb
