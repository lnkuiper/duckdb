#include "duckdb/optimizer/join_order/join_plan_domain_propagator.hpp"

#include "duckdb/common/limits.hpp"
#include "duckdb/common/unordered_map.hpp"
#include "duckdb/optimizer/join_order/join_order_operator.hpp"
#include "duckdb/optimizer/relation_statistics/relation_statistics_helper.hpp"

#include <math.h>

namespace duckdb {

struct JoinPlanDomainState {
	idx_t cardinality = 0;
	unordered_map<idx_t, RelationStats> relations;
};

struct EqualityDomainEvidence {
	idx_t class_index;
	idx_t intersection;
};

struct ColumnIntersectionEvidence {
	ColumnBinding binding;
	idx_t intersection;
};

static optional_ptr<RelationColumnStats> GetColumn(JoinPlanDomainState &state, const ColumnBinding &binding) {
	if (!binding.table_index.IsValid() || !binding.column_index.IsValid()) {
		return nullptr;
	}
	auto entry = state.relations.find(binding.table_index.index);
	if (entry == state.relations.end() || binding.column_index.GetIndex() >= entry->second.columns.size()) {
		return nullptr;
	}
	return entry->second.columns[binding.column_index.GetIndex()];
}

static optional_ptr<const RelationColumnStats> GetColumn(const JoinPlanDomainState &state,
                                                         const ColumnBinding &binding) {
	if (!binding.table_index.IsValid() || !binding.column_index.IsValid()) {
		return nullptr;
	}
	auto entry = state.relations.find(binding.table_index.index);
	if (entry == state.relations.end() || binding.column_index.GetIndex() >= entry->second.columns.size()) {
		return nullptr;
	}
	return entry->second.columns[binding.column_index.GetIndex()];
}

static optional_ptr<const RelationColumnStats> GetColumn(const vector<RelationStats> &relation_stats,
                                                         const ColumnBinding &binding) {
	if (!binding.table_index.IsValid() || !binding.column_index.IsValid() ||
	    binding.table_index.index >= relation_stats.size() ||
	    binding.column_index.GetIndex() >= relation_stats[binding.table_index.index].columns.size()) {
		return nullptr;
	}
	return relation_stats[binding.table_index.index].columns[binding.column_index.GetIndex()];
}

static void InvalidateCurrentDomainEvidence(JoinPlanDomainState &state) {
	for (auto &relation : state.relations) {
		for (auto &column : relation.second.columns) {
			column.current_domain_info.InvalidateStructuralEvidence();
		}
	}
}

static bool ContainsRelation(const JoinRelationSet &set, TableIndex table_index) {
	for (idx_t set_idx = 0; set_idx < set.count; set_idx++) {
		if (set.relations[set_idx].index == table_index.index) {
			return true;
		}
	}
	return false;
}

static optional_idx GetClassTotalDomain(const JoinEqualityClass &equality_class,
                                        const vector<RelationStats> &relation_stats) {
	vector<DistinctCount> domains;
	for (auto &binding : equality_class.columns) {
		auto column = GetColumn(relation_stats, binding);
		if (column) {
			domains.push_back(column->total_domain);
		}
	}
	auto domain = RelationStatisticsHelper::SelectTotalDomain(domains);
	return domain ? optional_idx(domain->distinct_count) : optional_idx();
}

struct EqualitySideEvidence {
	double fraction = 1;
	idx_t minimum_current = NumericLimits<idx_t>::Maximum();
};

static optional<EqualitySideEvidence> GetSideEqualityEvidence(const vector<RelationStats> &relation_stats,
                                                              const JoinEqualityClass &equality_class,
                                                              const JoinRelationSet &side, idx_t total_domain) {
	EqualitySideEvidence evidence;
	bool found = false;
	for (auto &binding : equality_class.columns) {
		if (!ContainsRelation(side, binding.table_index)) {
			continue;
		}
		auto column = GetColumn(relation_stats, binding);
		if (!column) {
			continue;
		}
		found = true;
		evidence.minimum_current = MinValue(evidence.minimum_current, column->current_domain.distinct_count);
		evidence.fraction *= MinValue(1.0, double(column->current_domain.distinct_count) / double(total_domain));
	}
	return found ? optional<EqualitySideEvidence>(evidence) : optional<EqualitySideEvidence>();
}

static idx_t EstimateIntersection(idx_t total_domain, double fraction, idx_t minimum_current,
                                  idx_t output_cardinality) {
	if (total_domain == 0 || fraction <= 0 || minimum_current == 0 || output_cardinality == 0) {
		return 0;
	}
	auto estimate = round(double(total_domain) * MinValue(fraction, 1.0));
	auto result = LossyNumericCast<idx_t>(MinValue(estimate, double(output_cardinality)));
	result = MinValue(result, MinValue(total_domain, minimum_current));
	return MaxValue<idx_t>(result, 1);
}

static idx_t EstimateIntersection(idx_t total_domain, idx_t left_current, idx_t right_current,
                                  idx_t output_cardinality) {
	if (total_domain == 0) {
		return 0;
	}
	auto fraction = MinValue(1.0, double(left_current) / double(total_domain)) *
	                MinValue(1.0, double(right_current) / double(total_domain));
	return EstimateIntersection(total_domain, fraction, MinValue(left_current, right_current), output_cardinality);
}

static idx_t EstimateRetainedRows(idx_t input_cardinality, idx_t output_cardinality, double coverage) {
	if (input_cardinality == 0 || output_cardinality == 0 || coverage <= 0) {
		return 0;
	}
	auto estimate = round(double(input_cardinality) * MinValue(coverage, 1.0));
	auto retained = LossyNumericCast<idx_t>(MinValue(estimate, double(input_cardinality)));
	retained = MinValue(retained, output_cardinality);
	return MaxValue<idx_t>(retained, 1);
}

static void ApplyRowRetention(JoinPlanDomainState &state, const JoinRelationSet &set, idx_t input_cardinality,
                              idx_t retained_cardinality) {
	for (idx_t set_idx = 0; set_idx < set.count; set_idx++) {
		auto relation_index = set.relations[set_idx].index;
		auto entry = state.relations.find(relation_index);
		if (entry == state.relations.end()) {
			continue;
		}
		for (auto &column : entry->second.columns) {
			column.current_domain = RelationStatisticsHelper::EstimateCurrentDomain(
			    column.current_domain, input_cardinality, retained_cardinality);
			column.current_domain.distinct_count =
			    MinValue(column.current_domain.distinct_count, column.total_domain.distinct_count);
			column.current_domain_info.InvalidateStructuralEvidence();
		}
	}
}

static void ApplyClassIntersection(JoinPlanDomainState &state, const JoinEqualityClass &equality_class,
                                   const JoinRelationSet &left, const JoinRelationSet &right, idx_t intersection,
                                   bool update_left, bool update_right) {
	for (auto &binding : equality_class.columns) {
		auto in_left = update_left && ContainsRelation(left, binding.table_index);
		auto in_right = update_right && ContainsRelation(right, binding.table_index);
		if (!in_left && !in_right) {
			continue;
		}
		auto column = GetColumn(state, binding);
		if (!column) {
			continue;
		}
		column->current_domain = DistinctCount(MinValue(column->current_domain.distinct_count, intersection),
		                                       DistinctCountSource::CARDINALITY);
		column->current_domain_info.InvalidateStructuralEvidence();
	}
}

static vector<reference<JoinPredicate>> GetNodePredicates(const DPJoinNode &node,
                                                          const JoinPredicateModel &predicate_model) {
	vector<reference<JoinPredicate>> result;
	unordered_set<idx_t> seen;
	for (auto predicate_ref : node.predicates) {
		auto &predicate = predicate_ref.get();
		if (seen.insert(predicate.GetIndex()).second) {
			result.push_back(predicate);
		}
	}
	if (!node.join_operator) {
		return result;
	}
	for (auto predicate_ref : predicate_model.GetPredicates()) {
		auto &predicate = predicate_ref.get();
		auto source_index = predicate.GetFilter().source_operator_index;
		if (source_index.IsValid() && source_index.GetIndex() == node.join_operator->index &&
		    seen.insert(predicate.GetIndex()).second) {
			result.push_back(predicate);
		}
	}
	return result;
}

static JoinOrderOperatorType GetNodeType(const DPJoinNode &node) {
	if (node.join_operator) {
		return node.join_operator->type;
	}
	if (node.generated_cross_product || node.predicates.empty()) {
		return JoinOrderOperatorType::CROSS_PRODUCT;
	}
	return JoinOrderOperatorType::INNER;
}

static optional<SemiAntiJoinCardinalityEstimate> EstimateSelectedSemiAntiJoin(const DPJoinNode &node,
                                                                              const JoinPredicateModel &predicate_model,
                                                                              const JoinPlanDomainState &left,
                                                                              const JoinPlanDomainState &right) {
	if (!node.join_operator || (node.join_operator->type != JoinOrderOperatorType::SEMI &&
	                            node.join_operator->type != JoinOrderOperatorType::ANTI)) {
		return {};
	}
	auto &join_operator = *node.join_operator;
	if (!JoinRelationSet::IsSubset(node.left_set, join_operator.left_total_set) ||
	    !JoinRelationSet::IsSubset(node.right_set, join_operator.right_total_set)) {
		return {};
	}

	vector<SemiAntiJoinDomain> domains;
	unordered_set<string> seen_domains;
	bool has_residual = false;
	auto join_type = join_operator.type == JoinOrderOperatorType::SEMI ? JoinType::SEMI : JoinType::ANTI;
	for (auto predicate_ref : predicate_model.GetPredicates()) {
		auto &predicate = predicate_ref.get();
		auto source_index = predicate.GetFilter().source_operator_index;
		if (!source_index.IsValid() || source_index.GetIndex() != join_operator.index) {
			continue;
		}
		if (predicate.GetJoinType() != join_type ||
		    (predicate.GetComparisonType() != ExpressionType::COMPARE_EQUAL &&
		     predicate.GetComparisonType() != ExpressionType::COMPARE_NOT_DISTINCT_FROM)) {
			has_residual = true;
			continue;
		}
		auto preserved_binding = predicate.GetStatsBinding(true);
		auto rhs_binding = predicate.GetStatsBinding(false);
		auto lhs_is_rhs = ContainsRelation(node.right_set, preserved_binding.table_index);
		auto rhs_is_rhs = ContainsRelation(node.right_set, rhs_binding.table_index);
		if (lhs_is_rhs == rhs_is_rhs) {
			has_residual = true;
			continue;
		}
		if (lhs_is_rhs) {
			std::swap(preserved_binding, rhs_binding);
		}
		auto preserved_column = GetColumn(left, preserved_binding);
		auto rhs_column = GetColumn(right, rhs_binding);
		if (!preserved_column || !rhs_column) {
			has_residual = true;
			continue;
		}
		auto evidence_key = preserved_binding.ToString() + "=" + rhs_binding.ToString();
		if (!seen_domains.insert(evidence_key).second) {
			continue;
		}
		auto rhs_current_domain = rhs_column->GetSemiAntiCurrentDomain();
		auto total_domain =
		    RelationStatisticsHelper::SelectTotalDomain({preserved_column->total_domain, rhs_column->total_domain});
		if (!rhs_current_domain.IsValid() || !total_domain || total_domain->distinct_count == 0) {
			has_residual = true;
			continue;
		}
		domains.push_back(
		    {total_domain->distinct_count, rhs_current_domain.GetIndex(), preserved_binding, rhs_binding});
	}

	double rhs_row_retention = 1;
	if (right.relations.size() == 1) {
		rhs_row_retention = right.relations.begin()->second.filter_strength;
	}
	return RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(
	    left.cardinality, right.cardinality, rhs_row_retention, join_type, domains, has_residual);
}

static optional<JoinPlanDomainState> PropagateNode(reference_map_t<JoinRelationSet, unique_ptr<DPJoinNode>> &plans,
                                                   const JoinPredicateModel &predicate_model,
                                                   const vector<RelationStats> &relation_stats, JoinRelationSet &set) {
	auto entry = plans.find(set);
	if (entry == plans.end() || entry->second->cardinality == DConstants::INVALID_INDEX) {
		return {};
	}
	auto &node = *entry->second;
	if (node.is_leaf) {
		if (node.set.count != 1 || node.set.relations[0].index >= relation_stats.size()) {
			return {};
		}
		JoinPlanDomainState result;
		result.cardinality = node.cardinality;
		auto relation_index = node.set.relations[0].index;
		result.relations.emplace(relation_index, relation_stats[relation_index]);
		return result;
	}

	auto left = PropagateNode(plans, predicate_model, relation_stats, node.left_set);
	auto right = PropagateNode(plans, predicate_model, relation_stats, node.right_set);
	if (!left || !right) {
		return {};
	}
	auto semi_anti_estimate = EstimateSelectedSemiAntiJoin(node, predicate_model, *left, *right);
	if (semi_anti_estimate && !semi_anti_estimate->IsFallback()) {
		node.cardinality = LossyNumericCast<idx_t>(semi_anti_estimate->cardinality);
	}
	JoinPlanDomainState result = std::move(*left);
	for (auto &entry : right->relations) {
		result.relations.emplace(entry.first, std::move(entry.second));
	}
	result.cardinality = node.cardinality;
	InvalidateCurrentDomainEvidence(result);

	auto node_type = GetNodeType(node);
	if (node_type == JoinOrderOperatorType::CROSS_PRODUCT) {
		return result;
	}
	auto predicates = GetNodePredicates(node, predicate_model);
	if (node_type == JoinOrderOperatorType::INNER) {
		vector<EqualityDomainEvidence> equality_evidence;
		unordered_set<idx_t> seen_classes;
		double left_coverage = 1;
		double right_coverage = 1;
		for (auto predicate_ref : predicates) {
			auto &predicate = predicate_ref.get();
			auto class_index = predicate.GetEqualityClassIndex();
			if (predicate.GetJoinType() != JoinType::INNER || !class_index.IsValid() ||
			    !seen_classes.insert(class_index.GetIndex()).second ||
			    class_index.GetIndex() >= predicate_model.GetEqualityClasses().size()) {
				continue;
			}
			auto &equality_class = predicate_model.GetEqualityClasses()[class_index.GetIndex()];
			auto total_domain = GetClassTotalDomain(equality_class, relation_stats);
			if (!total_domain.IsValid() || total_domain.GetIndex() == 0) {
				continue;
			}
			auto left_evidence =
			    GetSideEqualityEvidence(relation_stats, equality_class, node.left_set, total_domain.GetIndex());
			auto right_evidence =
			    GetSideEqualityEvidence(relation_stats, equality_class, node.right_set, total_domain.GetIndex());
			if (!left_evidence || !right_evidence) {
				continue;
			}
			left_coverage *= right_evidence->fraction;
			right_coverage *= left_evidence->fraction;
			equality_evidence.push_back(
			    {class_index.GetIndex(),
			     EstimateIntersection(total_domain.GetIndex(), left_evidence->fraction * right_evidence->fraction,
			                          MinValue(left_evidence->minimum_current, right_evidence->minimum_current),
			                          node.cardinality)});
		}
		if (equality_evidence.empty()) {
			ApplyRowRetention(result, node.left_set, left->cardinality, MinValue(left->cardinality, node.cardinality));
			ApplyRowRetention(result, node.right_set, right->cardinality,
			                  MinValue(right->cardinality, node.cardinality));
		} else {
			ApplyRowRetention(result, node.left_set, left->cardinality,
			                  EstimateRetainedRows(left->cardinality, node.cardinality, left_coverage));
			ApplyRowRetention(result, node.right_set, right->cardinality,
			                  EstimateRetainedRows(right->cardinality, node.cardinality, right_coverage));
			for (auto &evidence : equality_evidence) {
				auto &equality_class = predicate_model.GetEqualityClasses()[evidence.class_index];
				ApplyClassIntersection(result, equality_class, node.left_set, node.right_set, evidence.intersection,
				                       true, true);
			}
		}
		return result;
	}

	optional_ptr<const JoinRelationSet> preserved_set;
	optional_ptr<const JoinRelationSet> filtered_set;
	for (auto predicate_ref : predicates) {
		auto &predicate = predicate_ref.get();
		if (predicate.GetJoinType() != JoinType::LEFT && predicate.GetJoinType() != JoinType::SEMI &&
		    predicate.GetJoinType() != JoinType::ANTI) {
			continue;
		}
		if (JoinRelationSet::IsSubset(node.left_set, predicate.GetLeftSet())) {
			preserved_set = node.left_set;
			filtered_set = node.right_set;
		} else if (JoinRelationSet::IsSubset(node.right_set, predicate.GetLeftSet())) {
			preserved_set = node.right_set;
			filtered_set = node.left_set;
		}
		if (preserved_set) {
			break;
		}
	}
	if (!preserved_set || !filtered_set) {
		return result;
	}
	auto preserved_cardinality = preserved_set == &node.left_set ? left->cardinality : right->cardinality;
	auto filtered_cardinality = filtered_set == &node.left_set ? left->cardinality : right->cardinality;
	if (node_type == JoinOrderOperatorType::SEMI || node_type == JoinOrderOperatorType::ANTI) {
		ApplyRowRetention(result, *preserved_set, preserved_cardinality,
		                  MinValue(preserved_cardinality, node.cardinality));
		return result;
	}

	double filtered_coverage = 1;
	bool has_equality = false;
	vector<ColumnIntersectionEvidence> intersection_evidence;
	unordered_set<string> seen_evidence;
	for (auto predicate_ref : predicates) {
		auto &predicate = predicate_ref.get();
		if (predicate.GetJoinType() != JoinType::LEFT || !predicate.HasValidEqualityBindings()) {
			continue;
		}
		auto preserved_binding = predicate.GetEqualityBinding(true);
		auto filtered_binding = predicate.GetEqualityBinding(false);
		if (!ContainsRelation(*preserved_set, preserved_binding.table_index)) {
			std::swap(preserved_binding, filtered_binding);
		}
		auto preserved_column = GetColumn(result, preserved_binding);
		auto filtered_column = GetColumn(result, filtered_binding);
		if (!preserved_column || !filtered_column) {
			continue;
		}
		auto evidence_key = to_string(preserved_binding.table_index.index) + ":" +
		                    to_string(preserved_binding.column_index.GetIndex()) + "=" +
		                    to_string(filtered_binding.table_index.index) + ":" +
		                    to_string(filtered_binding.column_index.GetIndex());
		if (!seen_evidence.insert(evidence_key).second) {
			continue;
		}
		auto total_domain = RelationStatisticsHelper::SelectTotalDomain(
		    {preserved_column->total_domain, filtered_column->total_domain});
		if (!total_domain || total_domain->distinct_count == 0) {
			continue;
		}
		has_equality = true;
		filtered_coverage *= MinValue(1.0, double(preserved_column->current_domain.distinct_count) /
		                                       double(total_domain->distinct_count));
		auto intersection =
		    EstimateIntersection(total_domain->distinct_count, preserved_column->current_domain.distinct_count,
		                         filtered_column->current_domain.distinct_count, node.cardinality);
		intersection_evidence.push_back({filtered_binding, intersection});
	}
	auto retained_filtered = has_equality
	                             ? EstimateRetainedRows(filtered_cardinality, node.cardinality, filtered_coverage)
	                             : MinValue(filtered_cardinality, node.cardinality);
	ApplyRowRetention(result, *filtered_set, filtered_cardinality, retained_filtered);
	for (auto &evidence : intersection_evidence) {
		auto column = GetColumn(result, evidence.binding);
		if (column) {
			column->current_domain =
			    DistinctCount(MinValue(column->current_domain.distinct_count, evidence.intersection),
			                  DistinctCountSource::CARDINALITY);
			column->current_domain_info.InvalidateStructuralEvidence();
		}
	}
	return result;
}

optional<vector<RelationStats>>
JoinPlanDomainPropagator::Propagate(reference_map_t<JoinRelationSet, unique_ptr<DPJoinNode>> &plans,
                                    const JoinPredicateModel &predicate_model,
                                    const vector<RelationStats> &relation_stats, JoinRelationSet &root) {
	auto propagated = PropagateNode(plans, predicate_model, relation_stats, root);
	if (!propagated || propagated->relations.size() != relation_stats.size()) {
		return {};
	}
	vector<RelationStats> result;
	result.reserve(relation_stats.size());
	for (idx_t relation_idx = 0; relation_idx < relation_stats.size(); relation_idx++) {
		auto entry = propagated->relations.find(relation_idx);
		if (entry == propagated->relations.end()) {
			return {};
		}
		auto stats = std::move(entry->second);
		stats.cardinality = propagated->cardinality;
		for (auto &column : stats.columns) {
			column.current_domain.distinct_count =
			    MinValue(column.current_domain.distinct_count, propagated->cardinality);
		}
		vector<ColumnBinding> bindings;
		for (auto &column : stats.columns) {
			bindings.push_back(column.binding);
		}
		stats.Verify(bindings);
		result.push_back(std::move(stats));
	}
	return result;
}

} // namespace duckdb
