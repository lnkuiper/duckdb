#include "duckdb/optimizer/join_order/join_order_optimizer.hpp"

#include "duckdb/common/enums/join_type.hpp"
#include "duckdb/common/limits.hpp"
#include "duckdb/common/pair.hpp"
#include "duckdb/optimizer/join_order/cost_model.hpp"
#include "duckdb/optimizer/join_order/plan_enumerator.hpp"
#include "duckdb/planner/expression/list.hpp"
#include "duckdb/planner/operator/list.hpp"
#include "duckdb/optimizer/column_binding_replacer.hpp"

namespace duckdb {

JoinOrderOptimizer::JoinOrderOptimizer(ClientContext &context)
    : context(context), query_graph_manager(context), depth(1) {
}

JoinOrderOptimizer JoinOrderOptimizer::CreateChildOptimizer() {
	JoinOrderOptimizer child_optimizer(context);
	child_optimizer.materialized_cte_stats = materialized_cte_stats;
	child_optimizer.delim_scan_stats = delim_scan_stats;
	child_optimizer.depth = depth + 1;
	child_optimizer.recursive_cte_indexes = recursive_cte_indexes;
	return child_optimizer;
}

unique_ptr<LogicalOperator> RemoveUnnecessaryProjections::RemoveProjectionsChildren(unique_ptr<LogicalOperator> plan) {
	for (idx_t i = 0; i < plan->children.size(); i++) {
		plan->children[i] = RemoveProjections(std::move(plan->children[i]));
	}
	return plan;
}
unique_ptr<LogicalOperator> RemoveUnnecessaryProjections::RemoveProjections(unique_ptr<LogicalOperator> plan) {
	if (plan->type == LogicalOperatorType::LOGICAL_UNION || plan->type == LogicalOperatorType::LOGICAL_EXCEPT ||
	    plan->type == LogicalOperatorType::LOGICAL_INTERSECT ||
	    plan->type == LogicalOperatorType::LOGICAL_RECURSIVE_CTE ||
	    plan->type == LogicalOperatorType::LOGICAL_MATERIALIZED_CTE) {
		// guaranteed to find a projection under this that is meant to keep the column order in the presence of
		// an optimization done by build side probe side.
		for (idx_t i = 0; i < plan->children.size(); i++) {
			first_projection = true;
			plan->children[i] = RemoveProjections(std::move(plan->children[i]));
		}
		return plan;
	}
	if (plan->type != LogicalOperatorType::LOGICAL_PROJECTION) {
		return RemoveProjectionsChildren(std::move(plan));
	}
	// operator is a projection. Remove if possible
	if (first_projection) {
		first_projection = false;
		return RemoveProjectionsChildren(std::move(plan));
	}
	auto &proj = plan->Cast<LogicalProjection>();
	auto child_bindings = plan->children[0]->GetColumnBindings();
	if (proj.GetColumnBindings().size() != child_bindings.size()) {
		return plan;
	}
	idx_t binding_index = 0;
	for (auto &expr : proj.expressions) {
		if (expr->type != ExpressionType::BOUND_COLUMN_REF) {
			return plan;
		}
		auto &bound_ref = expr->Cast<BoundColumnRefExpression>();
		if (bound_ref.binding != child_bindings[binding_index]) {
			return plan;
		}
		binding_index++;
	}
	D_ASSERT(binding_index == plan->GetColumnBindings().size());
	// we have a projection where every expression is a bound column ref, and they are in the same order as the
	// bindings of the child. We can remove this projection
	binding_index = 0;
	for (auto &binding : plan->GetColumnBindings()) {
		replacer.replacement_bindings.push_back(ReplacementBinding(binding, child_bindings[binding_index]));
		binding_index++;
	}
	return RemoveProjectionsChildren(std::move(plan->children[0]));
}

RemoveUnnecessaryProjections::RemoveUnnecessaryProjections() {
	first_projection = true;
}

unique_ptr<LogicalOperator> JoinOrderOptimizer::Optimize(unique_ptr<LogicalOperator> plan,
                                                         optional_ptr<RelationStats> stats, bool remove_projections) {
	if (depth > query_graph_manager.context.config.max_expression_depth) {
		// Very deep plans will eventually consume quite some stack space
		// Returning the current plan is always a valid choice
		return plan;
	}

	// make sure query graph manager has not extracted a relation graph already
	if (remove_projections) {
		RemoveUnnecessaryProjections remover;
		plan = remover.RemoveProjections(std::move(plan));
		remover.replacer.VisitOperator(*plan);

		auto bindings = plan->GetColumnBindings();
	}

	LogicalOperator *op = plan.get();

	// extract the relations that go into the hyper graph.
	// We optimize the children of any non-reorderable operations we come across.
	bool reorderable = query_graph_manager.Build(*this, *op);
	// get relation_stats here since the reconstruction process will move all relations.
	auto relation_stats = query_graph_manager.relation_manager.GetRelationStats();
	unique_ptr<LogicalOperator> new_logical_plan = nullptr;

	if (reorderable) {
		// query graph now has filters and relations
		auto cost_model = CostModel(query_graph_manager);

		// Initialize a plan enumerator.
		auto plan_enumerator =
		    PlanEnumerator(query_graph_manager, cost_model, query_graph_manager.GetQueryGraphEdges());

		// Initialize the leaf/single node plans
		plan_enumerator.InitLeafPlans();
		plan_enumerator.SolveJoinOrder();
		// now reconstruct a logical plan from the query graph plan
		query_graph_manager.plans = &plan_enumerator.GetPlans();

		new_logical_plan = query_graph_manager.Reconstruct(std::move(plan));
	} else {
		new_logical_plan = std::move(plan);
		if (relation_stats.size() == 1) {
			new_logical_plan->estimated_cardinality = relation_stats.at(0).cardinality;
			new_logical_plan->has_estimated_cardinality = true;
		}
	}

	// Propagate up a stats object from the top of the new_logical_plan if stats exist.
	if (stats) {
		auto cardinality = new_logical_plan->EstimateCardinality(context);
		auto bindings = new_logical_plan->GetColumnBindings();
		auto new_stats = RelationStatisticsHelper::CombineStatsOfReorderableOperator(bindings, relation_stats);
		new_stats.cardinality = cardinality;
		RelationStatisticsHelper::CopyRelationStats(*stats, new_stats);
	} else {
		// starts recursively setting cardinality
		new_logical_plan->EstimateCardinality(context);
	}

	if (new_logical_plan->type == LogicalOperatorType::LOGICAL_EXPLAIN) {
		new_logical_plan->SetEstimatedCardinality(3);
	}

	return new_logical_plan;
}

void JoinOrderOptimizer::AddMaterializedCTEStats(idx_t index, RelationStats &&stats) {
	materialized_cte_stats.emplace(index, std::move(stats));
}

RelationStats JoinOrderOptimizer::GetMaterializedCTEStats(idx_t index) {
	auto it = materialized_cte_stats.find(index);
	if (it == materialized_cte_stats.end()) {
		throw InternalException("Unable to find materialized CTE stats with index %llu", index);
	}
	return it->second;
}

void JoinOrderOptimizer::AddDelimScanStats(RelationStats &stats) {
	delim_scan_stats = &stats;
}

RelationStats JoinOrderOptimizer::GetDelimScanStats() {
	if (!delim_scan_stats) {
		throw InternalException("Unable to find delim scan stats!");
	}
	return *delim_scan_stats;
}

} // namespace duckdb
