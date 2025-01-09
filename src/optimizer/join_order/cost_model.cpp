#include "duckdb/optimizer/join_order/join_node.hpp"
#include "duckdb/optimizer/join_order/join_order_optimizer.hpp"
#include "duckdb/optimizer/join_order/cost_model.hpp"

namespace duckdb {

CostModel::CostModel(QueryGraphManager &query_graph_manager)
    : query_graph_manager(query_graph_manager), cardinality_estimator() {
}

double CostModel::ComputeCost(DPJoinNode &left, DPJoinNode &right, NeighborInfo &neighbor_info) {
	bool has_left = false;
	for (auto &filter : neighbor_info.filters) {
		if (filter->join_type == JoinType::LEFT) {
			has_left = true;
		}
		if (filter->join_type == JoinType::INNER) {
			return ComputeInnerJoinCost(left, right);
		}
	}
	if (has_left) {
		return ComputeLeftJoinCost(left, right);
	}
	throw InternalException("umm what");
}

double CostModel::ComputeInnerJoinCost(DPJoinNode &left, DPJoinNode &right) {
	auto &combination = query_graph_manager.set_manager.Union(left.set, right.set);
	auto join_card = cardinality_estimator.EstimateCardinalityWithSet<double>(combination);
	auto join_cost = join_card;
	return join_cost + left.cost + right.cost;
}

double CostModel::ComputeLeftJoinCost(DPJoinNode &left, DPJoinNode &right) {
	auto &combination = query_graph_manager.set_manager.Union(left.set, right.set);
	auto join_card = cardinality_estimator.EstimateCardinalityWithSet<double>(combination);
	auto join_cost = join_card;
	if (join_cost < left.cardinality) {
	}
	return join_cost + left.cost + right.cost;
}

} // namespace duckdb
