#include "catch.hpp"
#include "test_helpers.hpp"

#include "duckdb.hpp"
#include "duckdb/optimizer/join_order/cardinality_estimator.hpp"
#include "duckdb/optimizer/join_order/join_order_optimizer.hpp"
#include "duckdb/optimizer/join_order/relation_manager.hpp"
#include "duckdb/optimizer/relation_statistics/relation_statistics_extractor.hpp"
#include "duckdb/optimizer/relation_statistics/relation_statistics_helper.hpp"
#include "duckdb/parser/parser.hpp"
#include "duckdb/planner/planner.hpp"
#include "duckdb/planner/expression/bound_between_expression.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/expression/bound_conjunction_expression.hpp"
#include "duckdb/planner/expression/bound_operator_expression.hpp"
#include "duckdb/planner/operator/logical_aggregate.hpp"
#include "duckdb/planner/operator/logical_column_data_get.hpp"
#include "duckdb/planner/operator/logical_cteref.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_cross_product.hpp"
#include "duckdb/planner/operator/logical_dummy_scan.hpp"
#include "duckdb/planner/operator/logical_distinct.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_limit.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"
#include "duckdb/planner/operator/logical_set_operation.hpp"
#include "duckdb/planner/operator/logical_window.hpp"

using namespace duckdb;

static RelationStats CreateStats(const vector<ColumnBinding> &bindings, const vector<idx_t> &distinct_counts,
                                 idx_t cardinality) {
	REQUIRE(bindings.size() == distinct_counts.size());
	RelationStats result;
	result.cardinality = cardinality;
	result.stats_initialized = true;
	for (idx_t column_idx = 0; column_idx < bindings.size(); column_idx++) {
		auto domain = DistinctCount(distinct_counts[column_idx], DistinctCountSource::EXACT);
		result.columns.emplace_back(bindings[column_idx], domain, domain, Identifier("column_" + to_string(column_idx)),
		                            CurrentDomainInfo(CurrentDomainProvenance::BASE));
	}
	result.Verify(bindings);
	return result;
}

static optional_ptr<LogicalOperator> FindOperator(LogicalOperator &op, LogicalOperatorType type) {
	if (op.type == type) {
		return op;
	}
	for (auto &child : op.children) {
		auto result = FindOperator(*child, type);
		if (result) {
			return result;
		}
	}
	return nullptr;
}

class UnsupportedOutputOperator : public LogicalOperator {
public:
	explicit UnsupportedOutputOperator(TableIndex output_index)
	    : LogicalOperator(LogicalOperatorType::LOGICAL_PIVOT), output_index(output_index) {
	}

public:
	vector<ColumnBinding> GetColumnBindings() override {
		return {ColumnBinding(output_index, ProjectionIndex(0))};
	}

protected:
	void ResolveTypes() override {
		types = {LogicalType::INTEGER};
	}

private:
	TableIndex output_index;
};

static unique_ptr<Expression> CreateComparison(ColumnBinding binding, ExpressionType type, int32_t value) {
	return BoundComparisonExpression::Create(type, make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, binding),
	                                         make_uniq<BoundConstantExpression>(Value::INTEGER(value)));
}

TEST_CASE("Current domains use multiplicity-aware row survival", "[optimizer][relation_statistics]") {
	auto unique_domain =
	    RelationStatisticsHelper::EstimateCurrentDomain(DistinctCount(1000, DistinctCountSource::EXACT), 1000, 100);
	REQUIRE(unique_domain.distinct_count == 100);
	REQUIRE(unique_domain.source == DistinctCountSource::CARDINALITY);

	auto repeated_domain =
	    RelationStatisticsHelper::EstimateCurrentDomain(DistinctCount(100, DistinctCountSource::EXACT), 1000, 100);
	REQUIRE(repeated_domain.distinct_count == 65);
	REQUIRE(RelationStatisticsHelper::EstimateCurrentDomain(DistinctCount(100, DistinctCountSource::EXACT), 1000, 1000)
	            .distinct_count == 100);
	REQUIRE(RelationStatisticsHelper::EstimateCurrentDomain(DistinctCount(100, DistinctCountSource::EXACT), 1000, 0)
	            .distinct_count == 0);
}

TEST_CASE("Direct filter domains recognize equality IN and conjunctions", "[optimizer][relation_statistics]") {
	auto binding = ColumnBinding(TableIndex(10), ProjectionIndex(0));
	auto equality = CreateComparison(binding, ExpressionType::COMPARE_EQUAL, 42);
	auto equality_bound = RelationStatisticsHelper::EstimateDirectFilterDomain(*equality, binding, 1000);
	REQUIRE(equality_bound.IsValid());
	REQUIRE(equality_bound.GetIndex() == 1);

	auto in = make_uniq<BoundOperatorExpression>(ExpressionType::COMPARE_IN, LogicalType::BOOLEAN);
	in->GetChildrenMutable().push_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, binding));
	in->GetChildrenMutable().push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(2)));
	in->GetChildrenMutable().push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(4)));
	in->GetChildrenMutable().push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(4)));
	auto in_bound = RelationStatisticsHelper::EstimateDirectFilterDomain(*in, binding, 1000);
	REQUIRE(in_bound.IsValid());
	REQUIRE(in_bound.GetIndex() == 2);

	auto range = make_uniq<BoundConjunctionExpression>(
	    ExpressionType::CONJUNCTION_AND, CreateComparison(binding, ExpressionType::COMPARE_GREATERTHANOREQUALTO, 10),
	    CreateComparison(binding, ExpressionType::COMPARE_LESSTHAN, 20));
	auto range_bound = RelationStatisticsHelper::EstimateDirectFilterDomain(*range, binding, 1000);
	REQUIRE_FALSE(range_bound.IsValid());

	auto conjunction = make_uniq<BoundConjunctionExpression>(
	    ExpressionType::CONJUNCTION_AND, CreateComparison(binding, ExpressionType::COMPARE_EQUAL, 10),
	    CreateComparison(binding, ExpressionType::COMPARE_LESSTHAN, 20));
	auto conjunction_bound = RelationStatisticsHelper::EstimateDirectFilterDomain(*conjunction, binding, 1000);
	REQUIRE(conjunction_bound.IsValid());
	REQUIRE(conjunction_bound.GetIndex() == 1);

	auto disjunction = make_uniq<BoundConjunctionExpression>(
	    ExpressionType::CONJUNCTION_OR, CreateComparison(binding, ExpressionType::COMPARE_EQUAL, 10),
	    CreateComparison(binding, ExpressionType::COMPARE_EQUAL, 20));
	auto or_bound = RelationStatisticsHelper::EstimateDirectFilterDomain(*disjunction, binding, 1000);
	REQUIRE(or_bound.IsValid());
	REQUIRE(or_bound.GetIndex() == 2);
}

TEST_CASE("Single-column groups use the current domain directly", "[optimizer][relation_statistics]") {
	auto child_binding = ColumnBinding(TableIndex(10), ProjectionIndex(0));
	auto child_stats = CreateStats({child_binding}, {1000}, 1000);
	child_stats.columns[0].current_domain = DistinctCount(100, DistinctCountSource::CARDINALITY);

	LogicalAggregate aggregate(TableIndex(20), TableIndex(21), {});
	aggregate.groups.push_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, child_binding));
	aggregate.grouping_sets.emplace_back();
	aggregate.grouping_sets.back().insert(ProjectionIndex(0));
	auto stats = RelationStatisticsHelper::ExtractAggregationStats(aggregate, child_stats);
	REQUIRE(stats);
	REQUIRE(stats->cardinality == 100);
	REQUIRE(stats->columns[0].total_domain.distinct_count == 1000);
	REQUIRE(stats->columns[0].current_domain.distinct_count == 100);
	REQUIRE(stats->columns[0].current_domain_info.is_unique);
	REQUIRE(stats->columns[0].current_domain_info.IsEligibleForSemiAnti());
}

TEST_CASE("Current-domain eligibility distinguishes direct and cross-column filters",
          "[optimizer][relation_statistics]") {
	auto bindings = LogicalOperator::GenerateColumnBindings(TableIndex(10), 2);
	auto child_stats = CreateStats(bindings, {1000, 1000}, 1000);

	vector<unique_ptr<Expression>> child_expressions;
	child_expressions.push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(1)));
	child_expressions.push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(2)));
	auto child = make_uniq<LogicalProjection>(TableIndex(10), std::move(child_expressions));
	child->children.push_back(make_uniq<LogicalDummyScan>(TableIndex(11)));
	LogicalFilter filter(CreateComparison(bindings[0], ExpressionType::COMPARE_EQUAL, 42));
	filter.children.push_back(std::move(child));

	auto stats = RelationStatisticsHelper::ExtractFilterStats(filter, child_stats);
	REQUIRE(stats);
	REQUIRE(stats->columns[0].current_domain_info.provenance == CurrentDomainProvenance::MODELED);
	REQUIRE(stats->columns[0].current_domain_info.direct_bound.IsValid());
	REQUIRE(stats->columns[0].current_domain_info.direct_bound.GetIndex() == 1);
	REQUIRE(stats->columns[0].GetSemiAntiCurrentDomain().GetIndex() == 1);
	REQUIRE(stats->columns[0].current_domain_info.IsEligibleForSemiAnti());
	REQUIRE(stats->columns[1].current_domain_info.provenance == CurrentDomainProvenance::MODELED);
	REQUIRE(!stats->columns[1].current_domain_info.direct_bound.IsValid());
	REQUIRE(!stats->columns[1].current_domain_info.IsEligibleForSemiAnti());

	auto range = BoundBetweenExpression::Create(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, bindings[0]),
	                                            make_uniq<BoundConstantExpression>(Value::INTEGER(1)),
	                                            make_uniq<BoundConstantExpression>(Value::INTEGER(500)), true, true);
	LogicalFilter range_filter(std::move(range));
	vector<unique_ptr<Expression>> range_child_expressions;
	range_child_expressions.push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(1)));
	range_child_expressions.push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(2)));
	auto range_child = make_uniq<LogicalProjection>(TableIndex(10), std::move(range_child_expressions));
	range_child->children.push_back(make_uniq<LogicalDummyScan>(TableIndex(12)));
	range_filter.children.push_back(std::move(range_child));
	auto range_stats = RelationStatisticsHelper::ExtractFilterStats(range_filter, child_stats);
	REQUIRE(range_stats);
	REQUIRE(range_stats->columns[0].current_domain.distinct_count == 200);
	REQUIRE_FALSE(range_stats->columns[0].current_domain_info.direct_bound.IsValid());
	REQUIRE_FALSE(range_stats->columns[0].GetSemiAntiCurrentDomain().IsValid());
}

TEST_CASE("Single-key aggregate uniqueness survives HAVING row estimation", "[optimizer][relation_statistics]") {
	auto child_binding = ColumnBinding(TableIndex(10), ProjectionIndex(0));
	auto child_stats = CreateStats({child_binding}, {1000}, 1000);
	vector<unique_ptr<Expression>> aggregates;
	aggregates.push_back(make_uniq<BoundConstantExpression>(Value::BIGINT(1)));
	auto aggregate = make_uniq<LogicalAggregate>(TableIndex(20), TableIndex(21), std::move(aggregates));
	aggregate->groups.push_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, child_binding));
	aggregate->grouping_sets.emplace_back();
	aggregate->grouping_sets.back().insert(ProjectionIndex(0));

	auto aggregate_stats = RelationStatisticsHelper::ExtractAggregationStats(*aggregate, child_stats);
	REQUIRE(aggregate_stats);
	auto aggregate_bindings = aggregate->GetColumnBindings();
	REQUIRE(aggregate_stats->columns[0].current_domain_info.is_unique);

	LogicalFilter having(CreateComparison(aggregate_bindings[1], ExpressionType::COMPARE_GREATERTHAN, 0));
	having.children.push_back(std::move(aggregate));
	auto having_stats = RelationStatisticsHelper::ExtractFilterStats(having, *aggregate_stats);
	REQUIRE(having_stats);
	REQUIRE(having_stats->columns[0].current_domain_info.provenance == CurrentDomainProvenance::MODELED);
	REQUIRE(having_stats->columns[0].current_domain_info.is_unique);
	REQUIRE(having_stats->columns[0].current_domain_info.IsEligibleForSemiAnti());
}

TEST_CASE("Limit statistics include offsets and percentages", "[optimizer][relation_statistics]") {
	auto binding = ColumnBinding(TableIndex(10), ProjectionIndex(0));
	auto child_stats = CreateStats({binding}, {1000}, 1000);

	LogicalLimit offset_limit(BoundLimitNode::ConstantValue(10), BoundLimitNode::ConstantValue(999));
	offset_limit.children.push_back(make_uniq<LogicalDummyScan>(TableIndex(10)));
	auto offset_stats = RelationStatisticsHelper::ExtractLimitStats(offset_limit, child_stats);
	REQUIRE(offset_stats);
	REQUIRE(offset_stats->cardinality == 1);
	REQUIRE(offset_stats->columns[0].current_domain.distinct_count == 1);

	LogicalLimit percent_limit(BoundLimitNode::ConstantPercentage(10), BoundLimitNode());
	percent_limit.children.push_back(make_uniq<LogicalDummyScan>(TableIndex(10)));
	auto percent_stats = RelationStatisticsHelper::ExtractLimitStats(percent_limit, child_stats);
	REQUIRE(percent_stats);
	REQUIRE(percent_stats->cardinality == 100);
	REQUIRE(percent_stats->columns[0].current_domain.distinct_count == 100);
}

TEST_CASE("Semi and anti joins use RHS domain coverage", "[optimizer][relation_statistics]") {
	vector<SemiAntiJoinDomain> single_domain {{1000, 100}};
	auto supported_semi =
	    RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 100, 1, JoinType::SEMI, single_domain, false);
	REQUIRE(supported_semi.cardinality == 100);
	REQUIRE(supported_semi.source == SemiAntiJoinCardinalitySource::SUPPORTED_DOMAIN);
	auto supported_anti =
	    RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 100, 1, JoinType::ANTI, single_domain, false);
	REQUIRE(supported_anti.cardinality == 900);
	REQUIRE(supported_anti.source == SemiAntiJoinCardinalitySource::SUPPORTED_DOMAIN);
	for (auto join_type : {JoinType::SEMI, JoinType::ANTI}) {
		auto residual =
		    RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 100, 1, join_type, single_domain, true);
		REQUIRE(residual.cardinality == 200);
		REQUIRE(residual.source == SemiAntiJoinCardinalitySource::FALLBACK);
		auto missing = RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 100, 1, join_type, {}, false);
		REQUIRE(missing.cardinality == 200);
		REQUIRE(missing.source == SemiAntiJoinCardinalitySource::FALLBACK);
	}
	REQUIRE(RelationStatisticsHelper::EstimateSemiAntiJoinFallback(1, 1, JoinType::SEMI).cardinality == 1);
	REQUIRE(RelationStatisticsHelper::EstimateSemiAntiJoinFallback(2, 1, JoinType::ANTI).cardinality == 1);

	auto left_binding = ColumnBinding(TableIndex(10), ProjectionIndex(0));
	auto second_left_binding = ColumnBinding(TableIndex(10), ProjectionIndex(1));
	auto right_binding = ColumnBinding(TableIndex(20), ProjectionIndex(0));
	auto second_right_binding = ColumnBinding(TableIndex(20), ProjectionIndex(1));
	vector<SemiAntiJoinDomain> composite_domains {{1000, 100, left_binding, right_binding},
	                                              {1000, 100, second_left_binding, second_right_binding}};
	REQUIRE(RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 1000, 0.5, JoinType::SEMI,
	                                                                  composite_domains, false)
	            .cardinality == 100);
	REQUIRE(RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 1000, 0.5, JoinType::SEMI,
	                                                                  composite_domains, true)
	            .cardinality == 200);
	vector<SemiAntiJoinDomain> shared_endpoint_domains {{1000, 100, left_binding, right_binding},
	                                                    {1000, 100, left_binding, second_right_binding}};
	REQUIRE(RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 1000, 0.5, JoinType::SEMI,
	                                                                  shared_endpoint_domains, false)
	            .cardinality == 10);
	auto empty_semi = RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 0, 0, JoinType::SEMI, {}, false);
	REQUIRE(empty_semi.cardinality == 0);
	REQUIRE(empty_semi.source == SemiAntiJoinCardinalitySource::EMPTY_RHS);
	auto empty_anti = RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 0, 0, JoinType::ANTI, {}, false);
	REQUIRE(empty_anti.cardinality == 1000);
	REQUIRE(empty_anti.source == SemiAntiJoinCardinalitySource::EMPTY_RHS);

	vector<SemiAntiJoinDomain> full_count_coverage {{1000, 1000}};
	REQUIRE(RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(1000, 1000, 1, JoinType::ANTI,
	                                                                  full_count_coverage, false)
	            .cardinality == 1);
}

TEST_CASE("Logical semi and anti estimates require supported equality evidence", "[optimizer][relation_statistics]") {
	auto left_binding = ColumnBinding(TableIndex(10), ProjectionIndex(0));
	auto right_binding = ColumnBinding(TableIndex(20), ProjectionIndex(0));
	auto second_right_binding = ColumnBinding(TableIndex(21), ProjectionIndex(0));
	auto left_stats = CreateStats({left_binding}, {1000}, 1000);
	auto right_stats = CreateStats({right_binding}, {1000}, 100);
	right_stats.columns[0].current_domain = DistinctCount(100, DistinctCountSource::EXACT);
	right_stats.columns[0].current_domain_info.is_unique = true;

	LogicalComparisonJoin null_safe_join(JoinType::SEMI);
	null_safe_join.conditions.emplace_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, left_binding),
	                                       make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, right_binding),
	                                       ExpressionType::COMPARE_NOT_DISTINCT_FROM);
	auto null_safe_estimate =
	    RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(null_safe_join, left_stats, right_stats);
	REQUIRE(null_safe_estimate.cardinality == 100);
	REQUIRE(null_safe_estimate.source == SemiAntiJoinCardinalitySource::SUPPORTED_DOMAIN);

	LogicalComparisonJoin residual_join(JoinType::ANTI);
	residual_join.conditions.push_back(null_safe_join.conditions[0].Copy());
	residual_join.conditions.emplace_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, left_binding),
	                                      make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, right_binding),
	                                      ExpressionType::COMPARE_LESSTHAN);
	auto residual_estimate =
	    RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(residual_join, left_stats, right_stats);
	REQUIRE(residual_estimate.cardinality == 200);
	REQUIRE(residual_estimate.source == SemiAntiJoinCardinalitySource::FALLBACK);

	auto empty_rhs = CreateStats({right_binding, second_right_binding}, {0, 0}, 0);
	for (auto join_type : {JoinType::SEMI, JoinType::ANTI}) {
		LogicalComparisonJoin empty_join(join_type);
		empty_join.conditions.push_back(null_safe_join.conditions[0].Copy());
		auto estimate = RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(empty_join, left_stats, empty_rhs);
		REQUIRE(estimate.cardinality == (join_type == JoinType::SEMI ? 0 : 1000));
		REQUIRE(estimate.source == SemiAntiJoinCardinalitySource::EMPTY_RHS);
	}
}

TEST_CASE("Nested empty SEMI remains exact in the enumeration cache", "[optimizer][relation_statistics]") {
	for (auto outer_join_type : {JoinType::SEMI, JoinType::ANTI}) {
		JoinRelationSetManager set_manager;
		auto &left_set = set_manager.GetJoinRelation(RelationIndex(0));
		auto &right_set = set_manager.GetJoinRelation(RelationIndex(1));
		auto &empty_set = set_manager.GetJoinRelation(RelationIndex(2));
		auto &inner_set = set_manager.Union(right_set, empty_set);
		auto &root_set = set_manager.Union(left_set, inner_set);
		auto left_binding = ColumnBinding(TableIndex(0), ProjectionIndex(0));
		auto right_binding = ColumnBinding(TableIndex(1), ProjectionIndex(0));
		auto empty_binding = ColumnBinding(TableIndex(2), ProjectionIndex(0));

		auto inner_comparison = BoundComparisonExpression::Create(
		    ExpressionType::COMPARE_EQUAL, make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, right_binding),
		    make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, empty_binding));
		FilterInfo inner_filter(std::move(inner_comparison), inner_set, 0, JoinType::SEMI);
		inner_filter.SetLeftSet(right_set);
		inner_filter.SetRightSet(empty_set);
		inner_filter.left_binding = right_binding;
		inner_filter.right_binding = empty_binding;
		inner_filter.source_operator_index = optional_idx(0);
		inner_filter.semantic_left_set = right_set;
		inner_filter.semantic_right_set = empty_set;

		auto outer_comparison = BoundComparisonExpression::Create(
		    ExpressionType::COMPARE_EQUAL, make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, left_binding),
		    make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, right_binding));
		FilterInfo outer_filter(std::move(outer_comparison), root_set, 1, outer_join_type);
		outer_filter.SetLeftSet(left_set);
		outer_filter.SetRightSet(inner_set);
		outer_filter.left_binding = left_binding;
		outer_filter.right_binding = right_binding;
		outer_filter.source_operator_index = optional_idx(1);
		outer_filter.semantic_left_set = left_set;
		outer_filter.semantic_right_set = inner_set;

		JoinPredicateModel predicate_model;
		predicate_model.RegisterPredicate(inner_filter, JoinPredicateClass::SEMI_ANTI_JOIN, right_binding,
		                                  empty_binding);
		predicate_model.RegisterPredicate(outer_filter, JoinPredicateClass::SEMI_ANTI_JOIN, left_binding,
		                                  right_binding);
		CardinalityEstimator estimator(set_manager, predicate_model);
		estimator.InitEquivalentRelations();
		auto left_stats = CreateStats({left_binding}, {1000}, 1000);
		auto right_stats = CreateStats({right_binding}, {100}, 100);
		auto empty_stats = CreateStats({empty_binding}, {0}, 0);
		estimator.InitCardinalityEstimatorProps(left_set, left_stats);
		estimator.InitCardinalityEstimatorProps(right_set, right_stats);
		estimator.InitCardinalityEstimatorProps(empty_set, empty_stats);

		REQUIRE(estimator.EstimateCardinalityWithSet<idx_t>(root_set) ==
		        (outer_join_type == JoinType::SEMI ? 0 : 1000));
		REQUIRE(estimator.EstimateCardinalityWithSet<idx_t>(inner_set) == 0);
	}
}

TEST_CASE("Empty multi-relation preserved SEMI remains exact", "[optimizer][relation_statistics]") {
	JoinRelationSetManager set_manager;
	auto &nonempty_set = set_manager.GetJoinRelation(RelationIndex(0));
	auto &empty_set = set_manager.GetJoinRelation(RelationIndex(1));
	auto &rhs_set = set_manager.GetJoinRelation(RelationIndex(2));
	auto &preserved_set = set_manager.Union(nonempty_set, empty_set);
	auto &root_set = set_manager.Union(preserved_set, rhs_set);
	auto nonempty_binding = ColumnBinding(TableIndex(0), ProjectionIndex(0));
	auto empty_binding = ColumnBinding(TableIndex(1), ProjectionIndex(0));
	auto rhs_binding = ColumnBinding(TableIndex(2), ProjectionIndex(0));

	auto comparison = BoundComparisonExpression::Create(
	    ExpressionType::COMPARE_EQUAL, make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, nonempty_binding),
	    make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, rhs_binding));
	FilterInfo filter(std::move(comparison), root_set, 0, JoinType::SEMI);
	filter.SetLeftSet(preserved_set);
	filter.SetRightSet(rhs_set);
	filter.left_binding = nonempty_binding;
	filter.right_binding = rhs_binding;
	filter.source_operator_index = optional_idx(0);
	filter.semantic_left_set = preserved_set;
	filter.semantic_right_set = rhs_set;

	JoinPredicateModel predicate_model;
	predicate_model.RegisterPredicate(filter, JoinPredicateClass::SEMI_ANTI_JOIN, nonempty_binding, rhs_binding);
	CardinalityEstimator estimator(set_manager, predicate_model);
	estimator.InitEquivalentRelations();
	auto nonempty_stats = CreateStats({nonempty_binding}, {100}, 100);
	auto empty_stats = CreateStats({empty_binding}, {0}, 0);
	auto rhs_stats = CreateStats({rhs_binding}, {100}, 100);
	estimator.InitCardinalityEstimatorProps(nonempty_set, nonempty_stats);
	estimator.InitCardinalityEstimatorProps(empty_set, empty_stats);
	estimator.InitCardinalityEstimatorProps(rhs_set, rhs_stats);

	REQUIRE(estimator.EstimateCardinalityWithSet<idx_t>(root_set) == 0);
	REQUIRE(estimator.EstimateCardinalityWithSet<idx_t>(preserved_set) == 100);
}

TEST_CASE("Relation statistics follow projection output bindings", "[optimizer][relation_statistics]") {
	auto child_table = TableIndex(10);
	auto child_bindings = LogicalOperator::GenerateColumnBindings(child_table, 2);
	auto child_stats = CreateStats(child_bindings, {7, 13}, 20);
	child_stats.columns[1].current_domain_info.UpdateDirectBound(13);
	child_stats.columns[1].current_domain_info.is_unique = true;

	vector<unique_ptr<Expression>> expressions;
	expressions.push_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, child_bindings[1]));
	expressions.push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(42)));
	expressions.push_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, child_bindings[0]));
	LogicalProjection projection(TableIndex(20), std::move(expressions));

	auto stats = RelationStatisticsHelper::ExtractProjectionStats(projection, child_stats);
	REQUIRE(stats);
	auto output_bindings = projection.GetColumnBindings();
	REQUIRE(stats->MatchesBindings(output_bindings));
	REQUIRE(stats->columns[0].total_domain.distinct_count == 13);
	REQUIRE(stats->columns[0].current_domain.distinct_count == 13);
	REQUIRE(stats->columns[0].current_domain_info.direct_bound.IsValid());
	REQUIRE(stats->columns[0].current_domain_info.direct_bound.GetIndex() == 13);
	REQUIRE(stats->columns[0].current_domain_info.is_unique);
	REQUIRE(stats->columns[1].total_domain.distinct_count == 1);
	REQUIRE(stats->columns[1].current_domain.source == DistinctCountSource::EXACT);
	REQUIRE(!stats->columns[1].current_domain_info.IsEligibleForSemiAnti());
	REQUIRE(stats->columns[2].current_domain.distinct_count == 7);
}

TEST_CASE("Aggregate group and result statistics use distinct bindings", "[optimizer][relation_statistics]") {
	auto child_table = TableIndex(10);
	auto child_bindings = LogicalOperator::GenerateColumnBindings(child_table, 2);
	auto child_stats = CreateStats(child_bindings, {5, 17}, 100);

	vector<unique_ptr<Expression>> aggregates;
	aggregates.push_back(make_uniq<BoundConstantExpression>(Value::BIGINT(1)));
	LogicalAggregate aggregate(TableIndex(20), TableIndex(21), std::move(aggregates));
	aggregate.groups.push_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, child_bindings[1]));
	aggregate.grouping_sets.emplace_back();
	aggregate.grouping_sets.back().insert(ProjectionIndex(0));

	auto stats = RelationStatisticsHelper::ExtractAggregationStats(aggregate, child_stats);
	REQUIRE(stats);
	auto output_bindings = aggregate.GetColumnBindings();
	REQUIRE(output_bindings.size() == 2);
	REQUIRE(stats->MatchesBindings(output_bindings));
	REQUIRE(stats->columns.size() == 2);
	REQUIRE(stats->columns[0].binding == ColumnBinding(TableIndex(20), ProjectionIndex(0)));
	REQUIRE(stats->columns[0].total_domain.distinct_count == 17);
	REQUIRE(stats->columns[0].current_domain.distinct_count == 17);
	REQUIRE(stats->columns[1].binding == ColumnBinding(TableIndex(21), ProjectionIndex(0)));

	DuckDB db;
	Connection connection(db);
	RelationManager relation_manager(*connection.context);
	REQUIRE(relation_manager.AddRelation(aggregate, nullptr, *stats));
	ColumnBinding normalized_group;
	ColumnBinding normalized_result;
	REQUIRE(relation_manager.TryNormalizeBinding(output_bindings[0], normalized_group));
	REQUIRE(relation_manager.TryNormalizeBinding(output_bindings[1], normalized_result));
	REQUIRE(normalized_group == ColumnBinding(TableIndex(0), ProjectionIndex(0)));
	REQUIRE(normalized_result == ColumnBinding(TableIndex(0), ProjectionIndex(1)));
}

static unique_ptr<LogicalOperator> CreateProjectedDummy(TableIndex input_index, TableIndex projection_index) {
	auto dummy = make_uniq<LogicalDummyScan>(input_index);
	vector<unique_ptr<Expression>> expressions;
	expressions.push_back(
	    make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, ColumnBinding(input_index, ProjectionIndex(0))));
	auto projection = make_uniq<LogicalProjection>(projection_index, std::move(expressions));
	projection->children.push_back(std::move(dummy));
	return std::move(projection);
}

static unique_ptr<LogicalOperator> CreateConstantProjection(TableIndex table_index, idx_t column_count) {
	vector<unique_ptr<Expression>> expressions;
	for (idx_t column_idx = 0; column_idx < column_count; column_idx++) {
		expressions.push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(NumericCast<int32_t>(column_idx))));
	}
	auto projection = make_uniq<LogicalProjection>(table_index, std::move(expressions));
	projection->children.push_back(make_uniq<LogicalDummyScan>(TableIndex(table_index.index + 100)));
	return std::move(projection);
}

TEST_CASE("Set operations invalidate SEMI and ANTI coverage evidence", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto left = CreateConstantProjection(TableIndex(10), 1);
	auto right = CreateConstantProjection(TableIndex(20), 1);
	auto left_stats = CreateStats(left->GetColumnBindings(), {10}, 10);
	auto right_stats = CreateStats(right->GetColumnBindings(), {10}, 10);
	left_stats.columns[0].current_domain_info.UpdateDirectBound(10);
	left_stats.columns[0].current_domain_info.is_unique = true;
	LogicalSetOperation setop(TableIndex(30), 1, std::move(left), std::move(right), LogicalOperatorType::LOGICAL_UNION,
	                          true);
	vector<reference<const RelationStats>> child_stats {left_stats, right_stats};
	auto stats = RelationStatisticsHelper::ExtractOperatorStats(setop, *connection.context, child_stats);
	REQUIRE(stats);
	REQUIRE(stats->columns[0].current_domain_info.provenance == CurrentDomainProvenance::MODELED);
	REQUIRE(!stats->columns[0].current_domain_info.IsEligibleForSemiAnti());
}

TEST_CASE("Relation statistics extraction is cached by operator", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto plan = CreateProjectedDummy(TableIndex(10), TableIndex(20));
	RelationStatsExtractor extractor(*connection.context);

	auto first = extractor.Extract(*plan);
	REQUIRE(first);
	REQUIRE(first->MatchesBindings(plan->GetColumnBindings()));
	REQUIRE(extractor.ExtractedOperatorCount() == 2);
	auto second = extractor.Extract(*plan);
	REQUIRE(second == first);
	REQUIRE(extractor.ExtractedOperatorCount() == 2);
}

TEST_CASE("Relation statistics follow SQL scan, filter and projection bindings", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	REQUIRE_NO_FAIL(connection.Query("CREATE TABLE stats_scan AS SELECT i % 10 AS k, i AS v FROM range(100) t(i)"));
	REQUIRE_NO_FAIL(connection.Query("ANALYZE stats_scan"));
	connection.BeginTransaction();
	auto plan = connection.ExtractPlan("SELECT v, k FROM stats_scan WHERE k = 3");
	RelationStatsExtractor extractor(*connection.context);

	auto stats = extractor.Extract(*plan);
	REQUIRE(stats);
	REQUIRE(stats->MatchesBindings(plan->GetColumnBindings()));
	REQUIRE(stats->cardinality > 0);
	REQUIRE(stats->cardinality < 100);
	REQUIRE(stats->columns.size() == 2);
	REQUIRE(stats->columns[0].name == Identifier("v"));
	REQUIRE(stats->columns[1].name == Identifier("k"));
	auto filtered_scan = FindOperator(*plan, LogicalOperatorType::LOGICAL_GET);
	REQUIRE(filtered_scan);
	auto filtered_scan_stats = extractor.Extract(*filtered_scan);
	REQUIRE(filtered_scan_stats);
	idx_t direct_columns = 0;
	for (auto &column : filtered_scan_stats->columns) {
		REQUIRE(column.current_domain_info.provenance == CurrentDomainProvenance::MODELED);
		if (column.current_domain_info.direct_bound.IsValid()) {
			direct_columns++;
			REQUIRE(column.current_domain_info.direct_bound.GetIndex() == 1);
			REQUIRE(column.GetSemiAntiCurrentDomain().GetIndex() == 1);
		} else {
			REQUIRE(!column.current_domain_info.IsEligibleForSemiAnti());
		}
	}
	REQUIRE(direct_columns == 1);
	auto scan_plan = connection.ExtractPlan("SELECT k, v FROM stats_scan");
	RelationStatsExtractor scan_extractor(*connection.context);
	auto scan = FindOperator(*scan_plan, LogicalOperatorType::LOGICAL_GET);
	REQUIRE(scan);
	auto scan_stats = scan_extractor.Extract(*scan);
	REQUIRE(scan_stats);
	REQUIRE(scan_stats->cardinality == 100);
	REQUIRE(scan_stats->columns.size() == 2);
	REQUIRE(scan_stats->columns[0].total_domain.distinct_count > 0);
	REQUIRE(scan_stats->columns[0].total_domain.distinct_count <= 100);
	REQUIRE(scan_stats->columns[0].current_domain.distinct_count <= scan_stats->columns[0].total_domain.distinct_count);
	REQUIRE(scan_stats->columns[0].current_domain_info.provenance == CurrentDomainProvenance::BASE);
	REQUIRE(!scan_stats->columns[0].current_domain_info.IsEligibleForSemiAnti());
	connection.Rollback();
}

TEST_CASE("Operator statistics remain aligned across unary and join operators", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto left_bindings = LogicalOperator::GenerateColumnBindings(TableIndex(10), 2);
	auto right_bindings = LogicalOperator::GenerateColumnBindings(TableIndex(20), 1);
	auto left_stats = CreateStats(left_bindings, {7, 11}, 100);
	auto right_stats = CreateStats(right_bindings, {3}, 4);
	left_stats.columns[0].current_domain_info.UpdateDirectBound(7);
	left_stats.columns[0].current_domain_info.is_unique = true;

	vector<unique_ptr<Expression>> projection_expressions;
	projection_expressions.push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(1)));
	projection_expressions.push_back(make_uniq<BoundConstantExpression>(Value::INTEGER(2)));
	auto filter_child = make_uniq<LogicalProjection>(TableIndex(10), std::move(projection_expressions));
	LogicalFilter filter;
	filter.projection_map.push_back(ProjectionIndex(1));
	filter.children.push_back(std::move(filter_child));
	vector<reference<const RelationStats>> unary_children {left_stats};
	auto filter_stats = RelationStatisticsHelper::ExtractOperatorStats(filter, *connection.context, unary_children);
	REQUIRE(filter_stats);
	REQUIRE(filter_stats->cardinality == 20);
	REQUIRE(filter_stats->MatchesBindings(filter.GetColumnBindings()));
	REQUIRE(filter_stats->columns[0].total_domain.distinct_count == 11);
	REQUIRE(filter_stats->columns[0].current_domain.distinct_count <= 11);

	auto window_child = CreateProjectedDummy(TableIndex(30), TableIndex(10));
	LogicalWindow window(TableIndex(40));
	window.expressions.push_back(make_uniq<BoundConstantExpression>(Value::BIGINT(1)));
	window.children.push_back(std::move(window_child));
	auto single_child_stats = CreateStats({ColumnBinding(TableIndex(10), ProjectionIndex(0))}, {7}, 100);
	vector<reference<const RelationStats>> window_children {single_child_stats};
	auto window_stats = RelationStatisticsHelper::ExtractOperatorStats(window, *connection.context, window_children);
	REQUIRE(window_stats);
	REQUIRE(window_stats->MatchesBindings(window.GetColumnBindings()));
	REQUIRE(window_stats->columns.size() == 2);

	auto cross_product = make_uniq<LogicalCrossProduct>(CreateConstantProjection(TableIndex(10), 2),
	                                                    CreateConstantProjection(TableIndex(20), 1));
	vector<reference<const RelationStats>> join_children {left_stats, right_stats};
	auto cross_stats =
	    RelationStatisticsHelper::ExtractOperatorStats(*cross_product, *connection.context, join_children);
	REQUIRE(cross_stats);
	REQUIRE(cross_stats->cardinality == 400);
	REQUIRE(cross_stats->MatchesBindings(cross_product->GetColumnBindings()));
	REQUIRE(!cross_stats->columns[0].current_domain_info.IsEligibleForSemiAnti());

	LogicalComparisonJoin semi_join(JoinType::SEMI);
	semi_join.children.push_back(CreateConstantProjection(TableIndex(10), 2));
	semi_join.children.push_back(CreateConstantProjection(TableIndex(20), 1));
	auto semi_stats = RelationStatisticsHelper::ExtractOperatorStats(semi_join, *connection.context, join_children);
	REQUIRE(semi_stats);
	REQUIRE(semi_stats->cardinality == 20);
	REQUIRE(semi_stats->MatchesBindings(semi_join.GetColumnBindings()));
	REQUIRE(!semi_stats->columns[0].current_domain_info.IsEligibleForSemiAnti());

	for (auto join_type :
	     {JoinType::INNER, JoinType::LEFT, JoinType::RIGHT, JoinType::OUTER, JoinType::ANTI, JoinType::RIGHT_ANTI}) {
		LogicalComparisonJoin join(join_type);
		join.children.push_back(CreateConstantProjection(TableIndex(10), 2));
		join.children.push_back(CreateConstantProjection(TableIndex(20), 1));
		auto join_stats = RelationStatisticsHelper::ExtractOperatorStats(join, *connection.context, join_children);
		REQUIRE(join_stats);
		REQUIRE(join_stats->MatchesBindings(join.GetColumnBindings()));
		idx_t expected_cardinality = left_stats.cardinality;
		if (join_type == JoinType::ANTI) {
			expected_cardinality = 20;
		} else if (join_type == JoinType::RIGHT_ANTI) {
			expected_cardinality = 1;
		}
		REQUIRE(join_stats->cardinality == expected_cardinality);
	}
}

TEST_CASE("Semi and anti operator statistics preserve the correct side", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto left_binding = ColumnBinding(TableIndex(10), ProjectionIndex(0));
	auto right_binding = ColumnBinding(TableIndex(20), ProjectionIndex(0));
	auto full_left_stats = CreateStats({left_binding}, {1000}, 1000);
	auto full_right_stats = CreateStats({right_binding}, {1000}, 1000);
	auto sparse_left_stats = full_left_stats;
	sparse_left_stats.cardinality = 100;
	sparse_left_stats.filter_strength = 0.1;
	sparse_left_stats.columns[0].current_domain = DistinctCount(100, DistinctCountSource::EXACT);
	auto sparse_right_stats = full_right_stats;
	sparse_right_stats.cardinality = 100;
	sparse_right_stats.filter_strength = 0.1;
	sparse_right_stats.columns[0].current_domain = DistinctCount(100, DistinctCountSource::EXACT);

	auto extract = [&](JoinType join_type, const RelationStats &left_stats, const RelationStats &right_stats,
	                   bool duplicate_evidence = false, bool has_residual = false) {
		LogicalComparisonJoin join(join_type);
		join.children.push_back(CreateConstantProjection(TableIndex(10), 1));
		join.children.push_back(CreateConstantProjection(TableIndex(20), 1));
		join.conditions.emplace_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, left_binding),
		                             make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, right_binding),
		                             ExpressionType::COMPARE_EQUAL);
		if (duplicate_evidence) {
			join.conditions.push_back(join.conditions[0].Copy());
		}
		if (has_residual) {
			join.conditions.emplace_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, left_binding),
			                             make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, right_binding),
			                             ExpressionType::COMPARE_LESSTHAN);
		}
		vector<reference<const RelationStats>> child_stats {left_stats, right_stats};
		return RelationStatisticsHelper::ExtractOperatorStats(join, *connection.context, child_stats);
	};

	auto ineligible_semi_stats = extract(JoinType::SEMI, full_left_stats, sparse_right_stats);
	REQUIRE(ineligible_semi_stats);
	REQUIRE(ineligible_semi_stats->cardinality == 200);
	auto ineligible_anti_stats = extract(JoinType::ANTI, full_left_stats, sparse_right_stats);
	REQUIRE(ineligible_anti_stats);
	REQUIRE(ineligible_anti_stats->cardinality == 200);
	sparse_left_stats.columns[0].current_domain_info.UpdateDirectBound(100);
	sparse_right_stats.columns[0].current_domain_info.UpdateDirectBound(100);
	auto semi_stats = extract(JoinType::SEMI, full_left_stats, sparse_right_stats);
	REQUIRE(semi_stats);
	REQUIRE(semi_stats->cardinality == 100);
	REQUIRE(semi_stats->filter_strength == Approx(0.1));
	auto anti_stats = extract(JoinType::ANTI, full_left_stats, sparse_right_stats);
	REQUIRE(anti_stats);
	REQUIRE(anti_stats->cardinality == 900);
	auto duplicate_stats = extract(JoinType::SEMI, full_left_stats, sparse_right_stats, true);
	REQUIRE(duplicate_stats);
	REQUIRE(duplicate_stats->cardinality == 100);
	auto residual_stats = extract(JoinType::SEMI, full_left_stats, sparse_right_stats, false, true);
	REQUIRE(residual_stats);
	REQUIRE(residual_stats->cardinality == 200);
	auto residual_anti_stats = extract(JoinType::ANTI, full_left_stats, sparse_right_stats, false, true);
	REQUIRE(residual_anti_stats);
	REQUIRE(residual_anti_stats->cardinality == 200);

	auto right_semi_stats = extract(JoinType::RIGHT_SEMI, sparse_left_stats, full_right_stats);
	REQUIRE(right_semi_stats);
	REQUIRE(right_semi_stats->cardinality == 100);
	REQUIRE(right_semi_stats->MatchesBindings({right_binding}));
	auto right_anti_stats = extract(JoinType::RIGHT_ANTI, sparse_left_stats, full_right_stats);
	REQUIRE(right_anti_stats);
	REQUIRE(right_anti_stats->cardinality == 900);
}

TEST_CASE("Delim semi and anti joins use domain coverage", "[optimizer][relation_statistics]") {
	auto left_binding = ColumnBinding(TableIndex(10), ProjectionIndex(0));
	auto right_binding = ColumnBinding(TableIndex(20), ProjectionIndex(0));
	auto left_stats = CreateStats({left_binding}, {1000}, 1000);
	auto right_stats = CreateStats({right_binding}, {1000}, 100);
	right_stats.filter_strength = 0.1;
	right_stats.columns[0].current_domain = DistinctCount(100, DistinctCountSource::EXACT);
	right_stats.columns[0].current_domain_info.UpdateDirectBound(100);
	LogicalComparisonJoin delim_join(JoinType::SEMI, LogicalOperatorType::LOGICAL_DELIM_JOIN);
	delim_join.conditions.emplace_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, left_binding),
	                                   make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, right_binding),
	                                   ExpressionType::COMPARE_EQUAL);
	REQUIRE(
	    RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(delim_join, left_stats, right_stats).cardinality ==
	    100);
	delim_join.join_type = JoinType::ANTI;
	REQUIRE(
	    RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(delim_join, left_stats, right_stats).cardinality ==
	    900);
}

TEST_CASE("Operator statistics reject incomplete child layouts", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto projection = CreateProjectedDummy(TableIndex(10), TableIndex(20));
	auto incomplete_stats = CreateStats({ColumnBinding(TableIndex(11), ProjectionIndex(0))}, {1}, 1);
	vector<reference<const RelationStats>> child_stats {incomplete_stats};
	REQUIRE_FALSE(RelationStatisticsHelper::ExtractOperatorStats(*projection, *connection.context, child_stats));
}

TEST_CASE("Unsupported output operators do not receive guessed statistics", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto child = make_uniq<LogicalDummyScan>(TableIndex(10));
	auto child_stats = CreateStats(child->GetColumnBindings(), {1}, 1);
	UnsupportedOutputOperator unsupported(TableIndex(20));
	unsupported.children.push_back(std::move(child));

	RelationStatsExtractor extractor(*connection.context);
	REQUIRE_FALSE(extractor.Extract(unsupported));

	RelationManager relation_manager(*connection.context);
	REQUIRE_FALSE(relation_manager.AddRelation(unsupported, nullptr, child_stats));
	REQUIRE_FALSE(relation_manager.HasCompleteStats());

	auto same_binding_child = make_uniq<LogicalDummyScan>(TableIndex(30));
	auto same_binding_stats = CreateStats(same_binding_child->GetColumnBindings(), {1}, 1);
	UnsupportedOutputOperator same_binding_unsupported(TableIndex(30));
	same_binding_unsupported.children.push_back(std::move(same_binding_child));
	RelationManager same_binding_manager(*connection.context);
	REQUIRE_FALSE(same_binding_manager.AddRelation(same_binding_unsupported, nullptr, same_binding_stats));
	REQUIRE_FALSE(same_binding_manager.HasCompleteStats());
}

TEST_CASE("Distinct statistics derive cardinality before relation registration", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto child = CreateConstantProjection(TableIndex(30), 2);
	auto child_stats = CreateStats(child->GetColumnBindings(), {5, 7}, 20);
	LogicalDistinct distinct(DistinctType::DISTINCT);
	distinct.children.push_back(std::move(child));

	vector<reference<const RelationStats>> children {child_stats};
	auto stats = RelationStatisticsHelper::ExtractOperatorStats(distinct, *connection.context, children);
	REQUIRE(stats);
	REQUIRE(stats->MatchesBindings(distinct.GetColumnBindings()));
	REQUIRE(stats->cardinality < child_stats.cardinality);
	REQUIRE(RelationStatisticsHelper::EstimateDistinctCardinality({DistinctCount(1, DistinctCountSource::EXACT)}, 1) ==
	        1);
	vector<unique_ptr<Expression>> matching_targets;
	matching_targets.push_back(
	    make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, distinct.GetColumnBindings()[0]));
	LogicalDistinct matching_distinct(std::move(matching_targets), DistinctType::DISTINCT_ON);
	matching_distinct.children.push_back(CreateConstantProjection(TableIndex(30), 2));
	auto matching_stats = RelationStatisticsHelper::ExtractDistinctStats(matching_distinct, child_stats);
	REQUIRE(matching_stats);
	REQUIRE(matching_stats->columns[0].current_domain_info.is_unique);
	REQUIRE(matching_stats->columns[0].current_domain_info.IsEligibleForSemiAnti());
	REQUIRE(!matching_stats->columns[1].current_domain_info.is_unique);

	auto single_child = CreateConstantProjection(TableIndex(40), 1);
	auto single_child_stats = CreateStats(single_child->GetColumnBindings(), {5}, 20);
	LogicalDistinct single_distinct(DistinctType::DISTINCT);
	single_distinct.children.push_back(std::move(single_child));
	vector<reference<const RelationStats>> single_children {single_child_stats};
	auto single_stats =
	    RelationStatisticsHelper::ExtractOperatorStats(single_distinct, *connection.context, single_children);
	REQUIRE(single_stats);
	REQUIRE(single_stats->columns[0].current_domain_info.is_unique);
	REQUIRE(single_stats->columns[0].current_domain_info.IsEligibleForSemiAnti());

	auto distinct_on_child = CreateConstantProjection(TableIndex(50), 1);
	auto output_binding = distinct_on_child->GetColumnBindings()[0];
	auto target_binding = ColumnBinding(TableIndex(51), ProjectionIndex(0));
	auto distinct_on_stats = CreateStats({output_binding, target_binding}, {5, 3}, 20);
	vector<unique_ptr<Expression>> distinct_on_targets;
	distinct_on_targets.push_back(make_uniq<BoundColumnRefExpression>(LogicalType::INTEGER, target_binding));
	LogicalDistinct unrelated_distinct(std::move(distinct_on_targets), DistinctType::DISTINCT_ON);
	unrelated_distinct.children.push_back(std::move(distinct_on_child));
	auto unrelated_stats = RelationStatisticsHelper::ExtractDistinctStats(unrelated_distinct, distinct_on_stats);
	REQUIRE(unrelated_stats);
	REQUIRE(!unrelated_stats->columns[0].current_domain_info.is_unique);
	REQUIRE(!unrelated_stats->columns[0].current_domain_info.IsEligibleForSemiAnti());

	RelationManager relation_manager(*connection.context);
	REQUIRE(relation_manager.AddRelation(distinct, nullptr, child_stats));
	REQUIRE(relation_manager.HasCompleteStats());
}

TEST_CASE("Generated scans and explain outputs have explicit statistics", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto in_plan = connection.ExtractPlan("SELECT * FROM range(10) t(i) WHERE i IN (1, 2, 3, 4, 5, 6)");
	auto chunk_get = FindOperator(*in_plan, LogicalOperatorType::LOGICAL_CHUNK_GET);
	REQUIRE(chunk_get);
	RelationStatsExtractor chunk_extractor(*connection.context);
	auto chunk_stats = chunk_extractor.Extract(*chunk_get);
	REQUIRE(chunk_stats);
	REQUIRE(chunk_stats->MatchesBindings(chunk_get->GetColumnBindings()));
	REQUIRE(chunk_stats->cardinality == chunk_get->Cast<LogicalColumnDataGet>().collection->Count());
	REQUIRE(chunk_stats->cardinality > 0);

	auto explain_plan = connection.ExtractPlan("EXPLAIN SELECT 42");
	REQUIRE(explain_plan->type == LogicalOperatorType::LOGICAL_EXPLAIN);
	auto unsupported_child = make_uniq<UnsupportedOutputOperator>(TableIndex(50));
	unsupported_child->children.push_back(make_uniq<LogicalDummyScan>(TableIndex(40)));
	explain_plan->children[0] = std::move(unsupported_child);
	RelationStatsExtractor explain_extractor(*connection.context);
	auto explain_stats = explain_extractor.Extract(*explain_plan);
	REQUIRE(explain_stats);
	REQUIRE(explain_stats->MatchesBindings(explain_plan->GetColumnBindings()));
	REQUIRE(explain_stats->cardinality == 3);
}

TEST_CASE("Relation statistics extraction rebinds CTE outputs", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto definition = CreateProjectedDummy(TableIndex(10), TableIndex(20));
	auto cte_index = TableIndex(30);
	LogicalCTERef cte_ref(TableIndex(40), cte_index, {LogicalType::INTEGER}, {Identifier("i")});
	RelationStatsExtractor extractor(*connection.context, [&](TableIndex index) -> optional_ptr<LogicalOperator> {
		return index == cte_index ? optional_ptr<LogicalOperator>(*definition) : nullptr;
	});

	auto stats = extractor.Extract(cte_ref);
	REQUIRE(stats);
	REQUIRE(stats->MatchesBindings(cte_ref.GetColumnBindings()));
	REQUIRE(stats->columns[0].binding == ColumnBinding(TableIndex(40), ProjectionIndex(0)));
	REQUIRE(extractor.ExtractedOperatorCount() == 3);
	REQUIRE(extractor.Extract(cte_ref) == stats);
	REQUIRE(extractor.ExtractedOperatorCount() == 3);
	LogicalCTERef second_ref(TableIndex(41), cte_index, {LogicalType::INTEGER}, {Identifier("i")});
	auto second_stats = extractor.Extract(second_ref);
	REQUIRE(second_stats);
	REQUIRE(second_stats->MatchesBindings(second_ref.GetColumnBindings()));
	REQUIRE(extractor.ExtractedOperatorCount() == 4);
}

TEST_CASE("Linear recursive CTE statistics preserve their output layout", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	connection.BeginTransaction();
	Parser parser(connection.context->GetParserOptions());
	parser.ParseQuery(R"(
		WITH RECURSIVE values(i) AS (
			SELECT 1
			UNION ALL
			SELECT i + 1 FROM values WHERE i < 10
		)
		SELECT i FROM values
	)");
	REQUIRE(parser.statements.size() == 1);
	Planner planner(*connection.context);
	planner.CreatePlan(std::move(parser.statements[0]));
	auto plan = std::move(planner.plan);
	JoinOrderOptimizer optimizer(*connection.context);
	RelationStats stats;
	plan = optimizer.Optimize(std::move(plan), stats);

	REQUIRE(stats.stats_initialized);
	REQUIRE(stats.cardinality == 1001);
	REQUIRE(stats.MatchesBindings(plan->GetColumnBindings()));
	REQUIRE(stats.columns.size() == 1);
	REQUIRE(stats.columns[0].total_domain.distinct_count == stats.cardinality);
	REQUIRE(stats.columns[0].current_domain.distinct_count == stats.cardinality);
	REQUIRE(stats.columns[0].total_domain.source == DistinctCountSource::CARDINALITY);
	REQUIRE(!stats.columns[0].current_domain_info.IsEligibleForSemiAnti());
	connection.Rollback();
}

TEST_CASE("Recursive CTE join terms use fixpoint cardinality fallbacks", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	connection.BeginTransaction();
	Parser parser(connection.context->GetParserOptions());
	parser.ParseQuery(R"(
		WITH RECURSIVE values(i) AS (
			SELECT 1
			UNION ALL
			SELECT i + 1 FROM values, (VALUES (1)) extra(j) WHERE i < 10
		)
		SELECT i FROM values
	)");
	REQUIRE(parser.statements.size() == 1);
	Planner planner(*connection.context);
	planner.CreatePlan(std::move(parser.statements[0]));
	auto plan = std::move(planner.plan);
	JoinOrderOptimizer optimizer(*connection.context);
	RelationStats stats;
	plan = optimizer.Optimize(std::move(plan), stats);

	REQUIRE(stats.stats_initialized);
	REQUIRE(stats.cardinality == 1001);
	REQUIRE(stats.MatchesBindings(plan->GetColumnBindings()));
	REQUIRE(stats.columns.size() == 1);
	REQUIRE(stats.columns[0].total_domain.distinct_count == stats.cardinality);
	REQUIRE(stats.columns[0].current_domain.distinct_count == stats.cardinality);
	REQUIRE(stats.columns[0].total_domain.source == DistinctCountSource::CARDINALITY);
	REQUIRE(!stats.columns[0].current_domain_info.IsEligibleForSemiAnti());
	connection.Rollback();
}

TEST_CASE("Relation statistics extraction rejects recurring and cyclic CTEs", "[optimizer][relation_statistics]") {
	DuckDB db;
	Connection connection(db);
	auto cte_index = TableIndex(30);
	LogicalCTERef recurring(TableIndex(40), cte_index, {LogicalType::INTEGER}, {Identifier("i")}, true);
	RelationStatsExtractor recurring_extractor(*connection.context,
	                                           [&](TableIndex) -> optional_ptr<LogicalOperator> { return recurring; });
	REQUIRE_FALSE(recurring_extractor.Extract(recurring));

	LogicalCTERef cyclic(TableIndex(50), cte_index, {LogicalType::INTEGER}, {Identifier("i")});
	RelationStatsExtractor cyclic_extractor(*connection.context,
	                                        [&](TableIndex) -> optional_ptr<LogicalOperator> { return cyclic; });
	REQUIRE_FALSE(cyclic_extractor.Extract(cyclic));
	REQUIRE(cyclic_extractor.ExtractedOperatorCount() == 1);
}
