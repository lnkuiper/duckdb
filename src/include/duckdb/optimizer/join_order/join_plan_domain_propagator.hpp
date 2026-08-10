//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/optimizer/join_order/join_plan_domain_propagator.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/optional.hpp"
#include "duckdb/common/reference_map.hpp"
#include "duckdb/optimizer/join_order/join_node.hpp"
#include "duckdb/optimizer/join_order/join_predicate.hpp"
#include "duckdb/optimizer/relation_statistics/relation_statistics.hpp"

namespace duckdb {

class JoinPlanDomainPropagator {
public:
	//! Propagate current domains through the selected tree and refine supported SEMI/ANTI node cardinalities.
	static optional<vector<RelationStats>> Propagate(reference_map_t<JoinRelationSet, unique_ptr<DPJoinNode>> &plans,
	                                                 const JoinPredicateModel &predicate_model,
	                                                 const vector<RelationStats> &relation_stats,
	                                                 JoinRelationSet &root);
};

} // namespace duckdb
