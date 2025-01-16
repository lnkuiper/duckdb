//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/optimizer/join_order/query_graph_manager.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/common.hpp"
#include "duckdb/common/enums/join_type.hpp"
#include "duckdb/common/optional_ptr.hpp"
#include "duckdb/common/pair.hpp"
#include "duckdb/common/unordered_map.hpp"
#include "duckdb/common/unordered_set.hpp"
#include "duckdb/common/vector.hpp"
#include "duckdb/optimizer/join_order/join_node.hpp"
#include "duckdb/optimizer/join_order/join_relation.hpp"
#include "duckdb/optimizer/join_order/query_graph.hpp"
#include "duckdb/optimizer/join_order/relation_manager.hpp"
#include "duckdb/planner/column_binding.hpp"
#include "duckdb/planner/logical_operator.hpp"

#include <functional>

namespace duckdb {

class QueryGraphEdges;

struct GenerateJoinRelation {
	GenerateJoinRelation(optional_ptr<JoinRelationSet> set, unique_ptr<LogicalOperator> op_p)
	    : set(set), op(std::move(op_p)) {
	}

	optional_ptr<JoinRelationSet> set;
	unique_ptr<LogicalOperator> op;
};


enum class FilterInfoApplicationRule : uint8_t {
	AS_JOIN = 0,                 // the join filter can join two relation sets
	AS_STRICT_FILTER = 1         // the join filter is not allowed to join two relation sets (even though it may seem like it)
	                             // needed for filters above left joins.
};


//! FilterInfo models stores filter information so that edges between relations can be made
//! with the original ColumnBinding information available so that the cardinality estimator can
//! view the statistics of the underlying base tables.
class FilterInfo {
public:
	FilterInfo(unique_ptr<Expression> filter, JoinRelationSet &set, idx_t filter_index,
	           JoinType join_type = JoinType::INNER)
	    : filter(std::move(filter)), set(set), filter_index(filter_index), join_type(join_type), application_rule(FilterInfoApplicationRule::AS_JOIN) {
	}

public:
	unique_ptr<Expression> filter;
	reference<JoinRelationSet> set;
	idx_t filter_index;
	JoinType join_type;
	optional_ptr<JoinRelationSet> left_relation_set;
	optional_ptr<JoinRelationSet> right_relation_set;
	//! required relations present in a join tree before applying
	//! this FilterInfo. Needed to make sure filters above left joins
	//! are applied above the left join
	FilterInfoApplicationRule application_rule;
	// TODO: change this to be a binding set
	ColumnBinding left_binding;
	ColumnBinding right_binding;

	void SetLeftSet(optional_ptr<JoinRelationSet> left_set_new);
	void SetRightSet(optional_ptr<JoinRelationSet> right_set_new);
};

//! The QueryGraphManager manages the process of extracting the reorderable and nonreorderable operations
//! from the logical plan and creating the intermediate structures needed by the plan enumerator.
//! When the plan enumerator finishes, the Query Graph Manger can then recreate the logical plan.
class QueryGraphManager {
public:
	explicit QueryGraphManager(ClientContext &context) : relation_manager(context), context(context) {
	}

	//! manage relations and the logical operators they represent
	RelationManager relation_manager;

	//! A structure holding all the created JoinRelationSet objects
	JoinRelationSetManager set_manager;

	ClientContext &context;

	//! Extract the join relations, optimizing non-reoderable relations when encountered
	bool Build(JoinOrderOptimizer &optimizer, LogicalOperator &op);

	//! Reconstruct the logical plan using the plan found by the plan enumerator
	unique_ptr<LogicalOperator> Reconstruct(unique_ptr<LogicalOperator> plan);

	//! Get a reference to the QueryGraphEdges structure that stores edges between
	//! nodes and hypernodes.
	const QueryGraphEdges &GetQueryGraphEdges() const;

	//! Get a list of the join filters in the join plan than eventually are
	//! transformed into the query graph edges
	const vector<unique_ptr<FilterInfo>> &GetFilterBindings() const;

	//! Plan enumerator may not find a full plan and therefore will need to create cross
	//! products to create edges.
	void CreateQueryGraphCrossProduct(JoinRelationSet &left, JoinRelationSet &right);

	//! A map to store the optimal join plan found for a specific JoinRelationSet*
	optional_ptr<const reference_map_t<JoinRelationSet, unique_ptr<DPJoinNode>>> plans;

private:
	vector<reference<LogicalOperator>> filter_operators;

	//! Filter information including the column_bindings that join filters
	//! used by the cardinality estimator to estimate distinct counts
	vector<unique_ptr<FilterInfo>> filters_and_bindings;

	QueryGraphEdges query_graph;

	// void GetColumnBinding(Expression &expression, ColumnBinding &binding);

	void CreateHyperGraphEdges();

	GenerateJoinRelation GenerateJoins(vector<unique_ptr<LogicalOperator>> &extracted_relations, JoinRelationSet &set);
};

} // namespace duckdb
