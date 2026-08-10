#include "duckdb/optimizer/relation_statistics/relation_statistics_helper.hpp"

#include "duckdb/common/operator/add.hpp"
#include "duckdb/common/operator/subtract.hpp"
#include "duckdb/common/types/hugeint.hpp"
#include "duckdb/common/types/value_map.hpp"
#include "duckdb/common/unordered_set.hpp"
#include "duckdb/planner/expression/bound_between_expression.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/expression/list.hpp"
#include "duckdb/planner/expression_iterator.hpp"
#include "duckdb/planner/operator/list.hpp"

#include <math.h>

namespace duckdb {

optional<DistinctCount> RelationStatisticsHelper::SelectTotalDomain(const vector<DistinctCount> &domains) {
	optional<DistinctCount> reliable;
	optional<DistinctCount> min_max;
	optional<DistinctCount> fallback;
	for (auto &domain : domains) {
		switch (domain.source) {
		case DistinctCountSource::HLL:
		case DistinctCountSource::EXACT:
			if (!reliable || domain.distinct_count > reliable->distinct_count) {
				reliable = domain;
			}
			break;
		case DistinctCountSource::MIN_MAX:
			if (!min_max || domain.distinct_count > min_max->distinct_count) {
				min_max = domain;
			}
			break;
		case DistinctCountSource::CARDINALITY:
			if (!fallback || domain.distinct_count < fallback->distinct_count) {
				fallback = domain;
			}
			break;
		}
	}
	if (reliable) {
		return reliable;
	}
	if (min_max) {
		return min_max;
	}
	return fallback;
}

struct DiscreteInterval {
	optional<int64_t> lower;
	optional<int64_t> upper;
};

static bool MatchesTargetColumn(const Expression &expression, optional<ColumnBinding> target_binding) {
	if (expression.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
		return !target_binding || expression.Cast<BoundColumnRefExpression>().Binding() == *target_binding;
	}
	return !target_binding && expression.GetExpressionClass() == ExpressionClass::BOUND_REF;
}

static bool TryGetDiscreteValue(const Expression &expression, int64_t &result) {
	if (expression.GetExpressionClass() != ExpressionClass::BOUND_CONSTANT) {
		return false;
	}
	auto &value = expression.Cast<BoundConstantExpression>().GetValue();
	if (value.IsNull()) {
		return false;
	}
	switch (value.type().id()) {
	case LogicalTypeId::BOOLEAN:
		result = value.GetValueUnsafe<bool>() ? 1 : 0;
		return true;
	case LogicalTypeId::TINYINT:
		result = value.GetValueUnsafe<int8_t>();
		return true;
	case LogicalTypeId::SMALLINT:
		result = value.GetValueUnsafe<int16_t>();
		return true;
	case LogicalTypeId::INTEGER:
		result = value.GetValueUnsafe<int32_t>();
		return true;
	case LogicalTypeId::BIGINT:
		result = value.GetValueUnsafe<int64_t>();
		return true;
	case LogicalTypeId::UTINYINT:
		result = value.GetValueUnsafe<uint8_t>();
		return true;
	case LogicalTypeId::USMALLINT:
		result = value.GetValueUnsafe<uint16_t>();
		return true;
	case LogicalTypeId::UINTEGER:
		result = value.GetValueUnsafe<uint32_t>();
		return true;
	case LogicalTypeId::UBIGINT: {
		auto input = value.GetValueUnsafe<uint64_t>();
		if (input > static_cast<uint64_t>(NumericLimits<int64_t>::Maximum())) {
			return false;
		}
		result = static_cast<int64_t>(input);
		return true;
	}
	case LogicalTypeId::HUGEINT:
		return Hugeint::TryCast(value.GetValueUnsafe<hugeint_t>(), result);
	case LogicalTypeId::DATE:
		result = value.GetValueUnsafe<date_t>().days;
		return true;
	default:
		return false;
	}
}

static ExpressionType FlipComparison(ExpressionType type) {
	switch (type) {
	case ExpressionType::COMPARE_LESSTHAN:
		return ExpressionType::COMPARE_GREATERTHAN;
	case ExpressionType::COMPARE_LESSTHANOREQUALTO:
		return ExpressionType::COMPARE_GREATERTHANOREQUALTO;
	case ExpressionType::COMPARE_GREATERTHAN:
		return ExpressionType::COMPARE_LESSTHAN;
	case ExpressionType::COMPARE_GREATERTHANOREQUALTO:
		return ExpressionType::COMPARE_LESSTHANOREQUALTO;
	default:
		return type;
	}
}

static bool TryGetComparisonInterval(const Expression &expression, optional<ColumnBinding> target_binding,
                                     DiscreteInterval &result) {
	if (!BoundComparisonExpression::IsComparison(expression)) {
		return false;
	}
	auto &comparison = expression.Cast<BoundFunctionExpression>();
	auto &left = BoundComparisonExpression::Left(comparison);
	auto &right = BoundComparisonExpression::Right(comparison);
	ExpressionType comparison_type = comparison.GetExpressionType();
	optional_ptr<const Expression> constant;
	if (MatchesTargetColumn(left, target_binding)) {
		constant = right;
	} else if (MatchesTargetColumn(right, target_binding)) {
		constant = left;
		comparison_type = FlipComparison(comparison_type);
	} else {
		return false;
	}
	int64_t value;
	if (!TryGetDiscreteValue(*constant, value)) {
		return false;
	}
	switch (comparison_type) {
	case ExpressionType::COMPARE_EQUAL:
		result.lower = value;
		result.upper = value;
		return true;
	case ExpressionType::COMPARE_GREATERTHAN:
		if (!TryAddOperator::Operation(value, int64_t(1), value)) {
			return false;
		}
		result.lower = value;
		return true;
	case ExpressionType::COMPARE_GREATERTHANOREQUALTO:
		result.lower = value;
		return true;
	case ExpressionType::COMPARE_LESSTHAN:
		if (!TrySubtractOperator::Operation(value, int64_t(1), value)) {
			return false;
		}
		result.upper = value;
		return true;
	case ExpressionType::COMPARE_LESSTHANOREQUALTO:
		result.upper = value;
		return true;
	default:
		return false;
	}
}

static void NormalizeIntervals(vector<DiscreteInterval> &intervals) {
	std::sort(intervals.begin(), intervals.end(), [](const DiscreteInterval &left, const DiscreteInterval &right) {
		if (!left.lower || !right.lower) {
			return !left.lower && right.lower;
		}
		if (*left.lower != *right.lower) {
			return *left.lower < *right.lower;
		}
		if (!left.upper || !right.upper) {
			return left.upper && !right.upper;
		}
		return *left.upper < *right.upper;
	});
	vector<DiscreteInterval> merged;
	for (auto &interval : intervals) {
		if (merged.empty()) {
			merged.push_back(interval);
			continue;
		}
		auto &last = merged.back();
		auto overlaps = !last.upper || !interval.lower || *interval.lower <= *last.upper;
		auto adjacent = last.upper && interval.lower && *last.upper < NumericLimits<int64_t>::Maximum() &&
		                *interval.lower == *last.upper + 1;
		if (!overlaps && !adjacent) {
			merged.push_back(interval);
			continue;
		}
		if (!last.upper || !interval.upper) {
			last.upper.reset();
		} else {
			last.upper = MaxValue(*last.upper, *interval.upper);
		}
	}
	intervals = std::move(merged);
}

static void IntersectIntervals(vector<DiscreteInterval> left, vector<DiscreteInterval> right,
                               vector<DiscreteInterval> &result) {
	NormalizeIntervals(left);
	NormalizeIntervals(right);
	idx_t left_idx = 0;
	idx_t right_idx = 0;
	while (left_idx < left.size() && right_idx < right.size()) {
		auto &left_interval = left[left_idx];
		auto &right_interval = right[right_idx];
		DiscreteInterval intersection;
		if (left_interval.lower && right_interval.lower) {
			intersection.lower = MaxValue(*left_interval.lower, *right_interval.lower);
		} else {
			intersection.lower = left_interval.lower ? left_interval.lower : right_interval.lower;
		}
		if (left_interval.upper && right_interval.upper) {
			intersection.upper = MinValue(*left_interval.upper, *right_interval.upper);
		} else {
			intersection.upper = left_interval.upper ? left_interval.upper : right_interval.upper;
		}
		if (!intersection.lower || !intersection.upper || *intersection.lower <= *intersection.upper) {
			result.push_back(intersection);
		}
		if (!left_interval.upper) {
			right_idx++;
		} else if (!right_interval.upper) {
			left_idx++;
		} else if (*left_interval.upper < *right_interval.upper) {
			left_idx++;
		} else {
			right_idx++;
		}
	}
}

static bool TryGetDiscreteIntervals(const Expression &expression, optional<ColumnBinding> target_binding,
                                    vector<DiscreteInterval> &result) {
	DiscreteInterval comparison_interval;
	if (TryGetComparisonInterval(expression, target_binding, comparison_interval)) {
		result.push_back(comparison_interval);
		return true;
	}
	if (expression.GetExpressionType() == ExpressionType::COMPARE_BETWEEN &&
	    expression.GetExpressionClass() == ExpressionClass::BOUND_FUNCTION) {
		auto &between = expression.Cast<BoundFunctionExpression>();
		if (!MatchesTargetColumn(BoundBetweenExpression::Input(between), target_binding)) {
			return false;
		}
		int64_t lower;
		int64_t upper;
		if (!TryGetDiscreteValue(BoundBetweenExpression::LowerBound(between), lower) ||
		    !TryGetDiscreteValue(BoundBetweenExpression::UpperBound(between), upper)) {
			return false;
		}
		if (!BoundBetweenExpression::LowerInclusive(between) && !TryAddOperator::Operation(lower, int64_t(1), lower)) {
			return false;
		}
		if (!BoundBetweenExpression::UpperInclusive(between) &&
		    !TrySubtractOperator::Operation(upper, int64_t(1), upper)) {
			return false;
		}
		if (lower <= upper) {
			result.push_back({lower, upper});
		}
		return true;
	}
	if (expression.GetExpressionType() == ExpressionType::COMPARE_IN &&
	    expression.GetExpressionClass() == ExpressionClass::BOUND_OPERATOR) {
		auto &children = expression.Cast<BoundOperatorExpression>().GetChildren();
		if (children.empty() || !MatchesTargetColumn(*children[0], target_binding)) {
			return false;
		}
		for (idx_t child_idx = 1; child_idx < children.size(); child_idx++) {
			int64_t value;
			if (!TryGetDiscreteValue(*children[child_idx], value)) {
				return false;
			}
			result.push_back({value, value});
		}
		return true;
	}
	if (expression.GetExpressionClass() != ExpressionClass::BOUND_CONJUNCTION) {
		return false;
	}
	auto &conjunction = expression.Cast<BoundConjunctionExpression>();
	if (conjunction.GetChildren().empty()) {
		return false;
	}
	if (expression.GetExpressionType() == ExpressionType::CONJUNCTION_OR) {
		for (auto &child : conjunction.GetChildren()) {
			vector<DiscreteInterval> child_intervals;
			if (!TryGetDiscreteIntervals(*child, target_binding, child_intervals)) {
				return false;
			}
			result.insert(result.end(), child_intervals.begin(), child_intervals.end());
		}
		NormalizeIntervals(result);
		return true;
	}
	if (expression.GetExpressionType() != ExpressionType::CONJUNCTION_AND) {
		return false;
	}
	result.push_back({});
	bool found_supported = false;
	for (auto &child : conjunction.GetChildren()) {
		vector<DiscreteInterval> child_intervals;
		if (!TryGetDiscreteIntervals(*child, target_binding, child_intervals)) {
			continue;
		}
		found_supported = true;
		vector<DiscreteInterval> intersection;
		IntersectIntervals(result, child_intervals, intersection);
		result = std::move(intersection);
	}
	return found_supported;
}

static optional_idx CountDiscreteIntervals(vector<DiscreteInterval> intervals, idx_t maximum_domain,
                                           optional<int64_t> minimum_value, optional<int64_t> maximum_value) {
	if (intervals.empty()) {
		return optional_idx(0);
	}
	for (auto &interval : intervals) {
		if (!interval.lower) {
			interval.lower = minimum_value;
		}
		if (!interval.upper) {
			interval.upper = maximum_value;
		}
		if (!interval.lower || !interval.upper) {
			return {};
		}
		interval.lower = MaxValue(*interval.lower, minimum_value ? *minimum_value : *interval.lower);
		interval.upper = MinValue(*interval.upper, maximum_value ? *maximum_value : *interval.upper);
	}
	intervals.erase(std::remove_if(intervals.begin(), intervals.end(),
	                               [](const DiscreteInterval &interval) { return *interval.lower > *interval.upper; }),
	                intervals.end());
	if (intervals.empty()) {
		return optional_idx(0);
	}
	std::sort(intervals.begin(), intervals.end(), [](const DiscreteInterval &left, const DiscreteInterval &right) {
		return *left.lower < *right.lower || (*left.lower == *right.lower && *left.upper < *right.upper);
	});
	vector<DiscreteInterval> merged;
	for (auto &interval : intervals) {
		if (merged.empty()) {
			merged.push_back(interval);
			continue;
		}
		auto &last = merged.back();
		auto adjacent = *last.upper < NumericLimits<int64_t>::Maximum() && *interval.lower == *last.upper + 1;
		if (*interval.lower <= *last.upper || adjacent) {
			last.upper = MaxValue(*last.upper, *interval.upper);
		} else {
			merged.push_back(interval);
		}
	}
	idx_t count = 0;
	for (auto &interval : merged) {
		hugeint_t span;
		if (!TrySubtractOperator::Operation(hugeint_t(*interval.upper), hugeint_t(*interval.lower), span) ||
		    !TryAddOperator::Operation(span, hugeint_t(1), span)) {
			return {};
		}
		uint64_t interval_count;
		if (!Hugeint::TryCast(span, interval_count) || interval_count >= maximum_domain - count) {
			return optional_idx(maximum_domain);
		}
		count += static_cast<idx_t>(interval_count);
	}
	return optional_idx(count);
}

optional_idx RelationStatisticsHelper::EstimateDirectFilterDomain(const Expression &expression,
                                                                  optional<ColumnBinding> target_binding,
                                                                  idx_t maximum_domain, optional<int64_t> minimum_value,
                                                                  optional<int64_t> maximum_value) {
	vector<DiscreteInterval> intervals;
	if (TryGetDiscreteIntervals(expression, target_binding, intervals)) {
		auto count = CountDiscreteIntervals(std::move(intervals), maximum_domain, minimum_value, maximum_value);
		if (count.IsValid()) {
			return count;
		}
	}

	if (BoundComparisonExpression::IsComparison(expression)) {
		auto &comparison = expression.Cast<BoundFunctionExpression>();
		if (comparison.GetExpressionType() != ExpressionType::COMPARE_EQUAL) {
			return {};
		}
		auto &left = BoundComparisonExpression::Left(comparison);
		auto &right = BoundComparisonExpression::Right(comparison);
		optional_ptr<const Expression> constant;
		if (MatchesTargetColumn(left, target_binding)) {
			constant = right;
		} else if (MatchesTargetColumn(right, target_binding)) {
			constant = left;
		}
		if (constant && constant->GetExpressionClass() == ExpressionClass::BOUND_CONSTANT &&
		    !constant->Cast<BoundConstantExpression>().GetValue().IsNull()) {
			return optional_idx(MinValue<idx_t>(maximum_domain, 1));
		}
		return {};
	}

	if (expression.GetExpressionType() == ExpressionType::COMPARE_IN &&
	    expression.GetExpressionClass() == ExpressionClass::BOUND_OPERATOR) {
		auto &children = expression.Cast<BoundOperatorExpression>().GetChildren();
		if (children.empty() || !MatchesTargetColumn(*children[0], target_binding)) {
			return {};
		}
		value_set_t constants;
		for (idx_t child_idx = 1; child_idx < children.size(); child_idx++) {
			if (children[child_idx]->GetExpressionClass() != ExpressionClass::BOUND_CONSTANT) {
				return {};
			}
			auto &value = children[child_idx]->Cast<BoundConstantExpression>().GetValue();
			if (!value.IsNull()) {
				constants.insert(value);
			}
		}
		return optional_idx(MinValue<idx_t>(maximum_domain, constants.size()));
	}

	if (expression.GetExpressionClass() != ExpressionClass::BOUND_CONJUNCTION) {
		return {};
	}
	auto &conjunction = expression.Cast<BoundConjunctionExpression>();
	if (expression.GetExpressionType() == ExpressionType::CONJUNCTION_AND) {
		optional_idx result;
		for (auto &child : conjunction.GetChildren()) {
			auto child_bound =
			    EstimateDirectFilterDomain(*child, target_binding, maximum_domain, minimum_value, maximum_value);
			if (!child_bound.IsValid()) {
				continue;
			}
			result = result.IsValid() ? MinValue(result.GetIndex(), child_bound.GetIndex()) : child_bound;
		}
		return result;
	}
	if (expression.GetExpressionType() != ExpressionType::CONJUNCTION_OR) {
		return {};
	}
	idx_t result = 0;
	for (auto &child : conjunction.GetChildren()) {
		auto child_bound =
		    EstimateDirectFilterDomain(*child, target_binding, maximum_domain, minimum_value, maximum_value);
		if (!child_bound.IsValid()) {
			return {};
		}
		if (child_bound.GetIndex() >= maximum_domain - result) {
			return optional_idx(maximum_domain);
		}
		result += child_bound.GetIndex();
	}
	return optional_idx(result);
}

struct ExpressionBinding {
	bool FoundExpression() const {
		return expression;
	}

	bool FoundColumnRef() const {
		return FoundExpression() && expression->GetExpressionType() == ExpressionType::BOUND_COLUMN_REF;
	}

	optional_ptr<Expression> expression;
	ColumnBinding child_binding;
	bool expression_is_constant = false;
};

static ExpressionBinding GetChildColumnBinding(Expression &expr) {
	ExpressionBinding result;
	if (expr.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
		result.expression = expr;
		result.child_binding = expr.Cast<BoundColumnRefExpression>().Binding();
		return result;
	}
	if (expr.IsFoldable()) {
		result.expression = expr;
		result.expression_is_constant = true;
		return result;
	}
	if (expr.IsVolatile()) {
		return result;
	}
	ExpressionIterator::EnumerateChildren(expr, [&](unique_ptr<Expression> &child) {
		if (result.FoundColumnRef()) {
			return;
		}
		auto child_result = GetChildColumnBinding(*child);
		if (child_result.FoundColumnRef()) {
			result = child_result;
		}
	});
	return result;
}

optional<RelationStats> RelationStatisticsHelper::ExtractProjectionStats(LogicalProjection &projection,
                                                                         const RelationStats &child_stats) {
	RelationStats result;
	result.cardinality = child_stats.cardinality;
	result.table_name = Identifier(projection.GetName());
	result.stats_initialized = true;
	auto bindings = projection.GetColumnBindings();
	D_ASSERT(bindings.size() == projection.expressions.size());
	for (idx_t expression_idx = 0; expression_idx < projection.expressions.size(); expression_idx++) {
		auto &expression = *projection.expressions[expression_idx];
		auto expression_binding = GetChildColumnBinding(expression);
		DistinctCount distinct_count(result.cardinality, DistinctCountSource::CARDINALITY);
		CurrentDomainInfo current_domain_info;
		if (expression_binding.expression_is_constant) {
			distinct_count = DistinctCount(MinValue<idx_t>(result.cardinality, 1), DistinctCountSource::EXACT);
		} else if (expression_binding.FoundColumnRef()) {
			auto child_column = child_stats.GetColumnStats(expression_binding.child_binding);
			if (!child_column) {
				return {};
			}
			distinct_count = child_column->total_domain;
		}
		auto current_domain = distinct_count;
		if (expression_binding.FoundColumnRef()) {
			auto child_column = child_stats.GetColumnStats(expression_binding.child_binding);
			current_domain = child_column->current_domain;
			if (expression.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
				current_domain_info = child_column->current_domain_info;
			}
		}
		result.columns.emplace_back(bindings[expression_idx], distinct_count, current_domain,
		                            Identifier(expression.GetName()), current_domain_info);
	}
	result.Verify(bindings);
	return result;
}

idx_t RelationStatisticsHelper::EstimateDistinctCardinality(const vector<DistinctCount> &distinct_counts,
                                                            idx_t input_cardinality) {
	if (distinct_counts.empty()) {
		return input_cardinality / 2;
	}
	if (distinct_counts.size() == 1) {
		return MinValue(distinct_counts[0].distinct_count, input_cardinality);
	}
	double product = 1;
	for (auto &distinct_count : distinct_counts) {
		product *= static_cast<double>(MaxValue<idx_t>(distinct_count.distinct_count, 1));
	}
	product *= pow(0.95, static_cast<double>(distinct_counts.size() - 1));
	const auto multiplier = 1.0 - exp(-static_cast<double>(input_cardinality) / product);
	const auto estimate = multiplier == 0 ? static_cast<double>(input_cardinality) : product * multiplier;
	auto result = LossyNumericCast<idx_t>(MinValue(estimate, static_cast<double>(input_cardinality)));
	return input_cardinality > 0 ? MaxValue<idx_t>(result, 1) : 0;
}

DistinctCount RelationStatisticsHelper::EstimateCurrentDomain(const DistinctCount &current_domain,
                                                              idx_t input_cardinality, idx_t output_cardinality) {
	auto count = MinValue(current_domain.distinct_count, input_cardinality);
	if (count == 0 || output_cardinality == 0 || input_cardinality == 0) {
		return DistinctCount(0, DistinctCountSource::CARDINALITY);
	}
	if (output_cardinality >= input_cardinality) {
		return DistinctCount(count, current_domain.source);
	}

	const auto row_survival = static_cast<double>(output_cardinality) / static_cast<double>(input_cardinality);
	const auto multiplicity = static_cast<double>(input_cardinality) / static_cast<double>(count);
	const auto retained_fraction = -expm1(multiplicity * log1p(-row_survival));
	const auto estimate = round(static_cast<double>(count) * retained_fraction);
	auto result = LossyNumericCast<idx_t>(MinValue(estimate, static_cast<double>(output_cardinality)));
	result = MinValue(result, count);
	return DistinctCount(MaxValue<idx_t>(result, 1), DistinctCountSource::CARDINALITY);
}

void RelationStatisticsHelper::ApplyCardinality(RelationStats &stats, idx_t output_cardinality) {
	const auto input_cardinality = stats.cardinality;
	for (auto &column : stats.columns) {
		column.current_domain = EstimateCurrentDomain(column.current_domain, input_cardinality, output_cardinality);
		column.current_domain.distinct_count =
		    MinValue(column.current_domain.distinct_count, column.total_domain.distinct_count);
		if (output_cardinality < input_cardinality) {
			column.current_domain_info.MarkModeled();
		}
	}
	if (input_cardinality > 0 && output_cardinality < input_cardinality) {
		stats.filter_strength *= double(output_cardinality) / double(input_cardinality);
	}
	stats.cardinality = output_cardinality;
}

static double RoundSemiAntiCardinality(double preserved_cardinality, double selectivity) {
	if (preserved_cardinality <= 0) {
		return 0;
	}
	auto estimate = round(preserved_cardinality * MaxValue(0.0, MinValue(selectivity, 1.0)));
	if (estimate <= 0) {
		estimate = 1;
	}
	return MinValue(estimate, preserved_cardinality);
}

SemiAntiJoinCardinalityEstimate RelationStatisticsHelper::EstimateSemiAntiJoinFallback(double preserved_cardinality,
                                                                                       idx_t rhs_cardinality,
                                                                                       JoinType join_type) {
	auto is_anti = join_type == JoinType::ANTI || join_type == JoinType::RIGHT_ANTI;
	if (rhs_cardinality == 0) {
		return {is_anti ? MaxValue(0.0, preserved_cardinality) : 0, SemiAntiJoinCardinalitySource::EMPTY_RHS};
	}
	return {RoundSemiAntiCardinality(preserved_cardinality, DEFAULT_SELECTIVITY),
	        SemiAntiJoinCardinalitySource::FALLBACK};
}

SemiAntiJoinCardinalityEstimate RelationStatisticsHelper::EstimateSemiAntiJoinCardinality(
    double preserved_cardinality, idx_t rhs_cardinality, double rhs_row_retention, JoinType join_type,
    const vector<SemiAntiJoinDomain> &domains, bool has_residual) {
	auto fallback = EstimateSemiAntiJoinFallback(preserved_cardinality, rhs_cardinality, join_type);
	if (preserved_cardinality <= 0 || rhs_cardinality == 0 || has_residual || domains.empty()) {
		return fallback;
	}

	double match_fraction = 1;
	double strongest_fraction = 1;
	idx_t strongest_total_domain = 0;
	idx_t maximum_total_domain = 0;
	unordered_set<string> preserved_bindings;
	unordered_set<string> rhs_bindings;
	bool has_complete_binding_evidence = !domains.empty();
	for (auto &domain : domains) {
		if (domain.total_domain == 0) {
			return fallback;
		}
		double fraction = MinValue(1.0, double(domain.rhs_current_domain) / double(domain.total_domain));
		match_fraction *= fraction;
		if (fraction < strongest_fraction) {
			strongest_fraction = fraction;
			strongest_total_domain = domain.total_domain;
		}
		maximum_total_domain = MaxValue(maximum_total_domain, domain.total_domain);
		if (!domain.preserved_binding || !domain.rhs_binding) {
			has_complete_binding_evidence = false;
			continue;
		}
		preserved_bindings.insert(domain.preserved_binding->ToString());
		rhs_bindings.insert(domain.rhs_binding->ToString());
	}
	if (domains.size() >= 2 && has_complete_binding_evidence && preserved_bindings.size() >= 2 &&
	    rhs_bindings.size() >= 2 && maximum_total_domain > 0) {
		auto retention = MaxValue(0.0, MinValue(rhs_row_retention, 1.0));
		auto rhs_base_cardinality = retention > 0 ? double(rhs_cardinality) / retention : double(rhs_cardinality);
		if (rhs_base_cardinality <= double(maximum_total_domain) * 8 &&
		    (strongest_total_domain == 0 || double(strongest_total_domain) <= rhs_base_cardinality)) {
			match_fraction = MaxValue(match_fraction, MinValue(strongest_fraction, retention));
		}
	}
	match_fraction = MaxValue(0.0, MinValue(match_fraction, 1.0));
	auto is_anti = join_type == JoinType::ANTI || join_type == JoinType::RIGHT_ANTI;
	auto selectivity = is_anti ? 1 - match_fraction : match_fraction;
	return {RoundSemiAntiCardinality(preserved_cardinality, selectivity),
	        SemiAntiJoinCardinalitySource::SUPPORTED_DOMAIN};
}

optional<RelationStats> RelationStatisticsHelper::ExtractAggregationStats(LogicalAggregate &aggregate,
                                                                          const RelationStats &child_stats) {
	vector<DistinctCount> cardinality_counts;
	for (auto &grouping_set : aggregate.grouping_sets) {
		vector<DistinctCount> set_counts;
		for (auto group_idx : grouping_set) {
			auto &group = aggregate.GetGroupExpression(group_idx);
			if (group.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF) {
				continue;
			}
			auto column = child_stats.GetColumnStats(group.Cast<BoundColumnRefExpression>().Binding());
			if (!column) {
				return {};
			}
			auto count = column->current_domain;
			count.distinct_count = MaxValue<idx_t>(count.distinct_count, 1);
			set_counts.push_back(count);
		}
		if (set_counts.size() > cardinality_counts.size()) {
			cardinality_counts = std::move(set_counts);
		}
	}

	RelationStats result;
	result.cardinality =
	    aggregate.groups.empty() ? 1 : EstimateDistinctCardinality(cardinality_counts, child_stats.cardinality);
	result.table_name = Identifier(aggregate.GetName());
	result.stats_initialized = true;
	auto bindings = aggregate.GetColumnBindings();
	for (idx_t group_idx = 0; group_idx < aggregate.groups.size(); group_idx++) {
		auto &group = *aggregate.groups[group_idx];
		DistinctCount distinct_count(result.cardinality, DistinctCountSource::CARDINALITY);
		CurrentDomainInfo current_domain_info;
		if (group.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
			auto child_column = child_stats.GetColumnStats(group.Cast<BoundColumnRefExpression>().Binding());
			if (!child_column) {
				return {};
			}
			distinct_count = child_column->total_domain;
			current_domain_info = child_column->current_domain_info;
		}
		auto current_domain = distinct_count;
		if (group.GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF) {
			current_domain =
			    child_stats.GetColumnStats(group.Cast<BoundColumnRefExpression>().Binding())->current_domain;
			current_domain.distinct_count = MinValue(current_domain.distinct_count, result.cardinality);
		}
		if (aggregate.groups.size() == 1 && aggregate.grouping_sets.size() == 1 &&
		    aggregate.grouping_sets[0].size() == 1 && aggregate.grouping_sets[0].count(ProjectionIndex(0)) > 0) {
			current_domain_info.is_unique = true;
		}
		result.columns.emplace_back(bindings[result.columns.size()], distinct_count, current_domain,
		                            Identifier(group.GetName()), current_domain_info);
	}
	for (auto &expression : aggregate.expressions) {
		result.columns.emplace_back(bindings[result.columns.size()],
		                            DistinctCount(result.cardinality, DistinctCountSource::CARDINALITY),
		                            Identifier(expression->GetName()));
	}
	for (idx_t grouping_idx = 0; grouping_idx < aggregate.grouping_functions.size(); grouping_idx++) {
		auto grouping_count = MinValue<idx_t>(result.cardinality, MaxValue<idx_t>(aggregate.grouping_sets.size(), 1));
		result.columns.emplace_back(bindings[result.columns.size()],
		                            DistinctCount(grouping_count, DistinctCountSource::CARDINALITY),
		                            Identifier("grouping"));
	}
	result.Verify(bindings);
	return result;
}

optional<RelationStats> RelationStatisticsHelper::ExtractWindowStats(LogicalWindow &window,
                                                                     const RelationStats &child_stats) {
	RelationStats result;
	result.cardinality = child_stats.cardinality;
	result.table_name = Identifier(window.GetName());
	result.stats_initialized = true;
	for (auto &binding : window.GetColumnBindings()) {
		auto child_column = child_stats.GetColumnStats(binding);
		if (child_column) {
			result.columns.emplace_back(binding, child_column->total_domain, child_column->current_domain,
			                            child_column->name, child_column->current_domain_info);
		} else if (binding.table_index == window.window_index) {
			result.columns.emplace_back(binding, DistinctCount(result.cardinality, DistinctCountSource::CARDINALITY),
			                            Identifier("window"));
		} else {
			return {};
		}
	}
	result.Verify(window.GetColumnBindings());
	return result;
}

static optional<DistinctCount> GetDistinctTargetCount(Expression &target, const RelationStats &child_stats) {
	switch (target.GetExpressionClass()) {
	case ExpressionClass::BOUND_COLUMN_REF: {
		auto column = child_stats.GetColumnStats(target.Cast<BoundColumnRefExpression>().Binding());
		return column ? optional<DistinctCount>(column->current_domain) : optional<DistinctCount>();
	}
	case ExpressionClass::BOUND_REF: {
		auto index = target.Cast<BoundReferenceExpression>().Index();
		return index < child_stats.columns.size() ? optional<DistinctCount>(child_stats.columns[index].current_domain)
		                                          : optional<DistinctCount>();
	}
	case ExpressionClass::BOUND_CONSTANT:
		return DistinctCount(MinValue<idx_t>(child_stats.cardinality, 1), DistinctCountSource::EXACT);
	default:
		return {};
	}
}

static optional<ColumnBinding> GetDistinctTargetBinding(Expression &target, const RelationStats &child_stats) {
	switch (target.GetExpressionClass()) {
	case ExpressionClass::BOUND_COLUMN_REF:
		return target.Cast<BoundColumnRefExpression>().Binding();
	case ExpressionClass::BOUND_REF: {
		auto index = target.Cast<BoundReferenceExpression>().Index();
		return index < child_stats.columns.size() ? optional<ColumnBinding>(child_stats.columns[index].binding)
		                                          : optional<ColumnBinding>();
	}
	default:
		return {};
	}
}

optional<RelationStats> RelationStatisticsHelper::ExtractDistinctStats(LogicalDistinct &distinct,
                                                                       const RelationStats &child_stats) {
	auto result = ProjectOutputStats(child_stats, distinct);
	if (!result) {
		return {};
	}
	vector<DistinctCount> distinct_counts;
	if (distinct.distinct_targets.empty()) {
		for (auto &column : child_stats.columns) {
			distinct_counts.push_back(column.current_domain);
		}
	} else {
		for (auto &target : distinct.distinct_targets) {
			auto count = GetDistinctTargetCount(*target, child_stats);
			if (!count) {
				return result;
			}
			distinct_counts.push_back(*count);
		}
	}
	result->cardinality = EstimateDistinctCardinality(distinct_counts, child_stats.cardinality);
	for (auto &column : result->columns) {
		column.current_domain.distinct_count = MinValue(column.current_domain.distinct_count, result->cardinality);
		column.current_domain_info.MarkModeled();
	}
	if (result->columns.size() == 1 && distinct.distinct_targets.empty()) {
		result->columns[0].current_domain_info.is_unique = true;
	} else if (distinct.distinct_targets.size() == 1) {
		auto target_binding = GetDistinctTargetBinding(*distinct.distinct_targets[0], child_stats);
		if (target_binding) {
			for (auto &column : result->columns) {
				if (column.binding == *target_binding) {
					column.current_domain_info.is_unique = true;
					break;
				}
			}
		}
	}
	return result;
}

optional<RelationStats> RelationStatisticsHelper::ExtractFilterStats(LogicalFilter &filter,
                                                                     const RelationStats &child_stats) {
	auto result = ProjectOutputStats(child_stats, filter);
	if (!result) {
		return {};
	}
	auto cardinality = child_stats.cardinality;
	if (cardinality > 0) {
		cardinality = MaxValue<idx_t>(LossyNumericCast<idx_t>(double(cardinality) * DEFAULT_SELECTIVITY), 1);
	}
	ApplyCardinality(*result, cardinality);
	for (auto &column : result->columns) {
		for (auto &expression : filter.expressions) {
			auto direct_bound =
			    EstimateDirectFilterDomain(*expression, column.binding, column.total_domain.distinct_count);
			if (!direct_bound.IsValid()) {
				continue;
			}
			column.current_domain =
			    DistinctCount(MinValue(column.current_domain.distinct_count, direct_bound.GetIndex()),
			                  DistinctCountSource::CARDINALITY);
			column.current_domain_info.UpdateDirectBound(direct_bound.GetIndex());
		}
	}
	return result;
}

optional<RelationStats> RelationStatisticsHelper::ExtractLimitStats(LogicalLimit &limit,
                                                                    const RelationStats &child_stats) {
	auto result = ProjectOutputStats(child_stats, limit);
	if (!result) {
		return {};
	}
	auto available_cardinality = child_stats.cardinality;
	if (limit.offset_val.Type() == LimitNodeType::CONSTANT_VALUE) {
		auto offset = limit.offset_val.GetConstantValue();
		available_cardinality = offset >= available_cardinality ? 0 : available_cardinality - offset;
	}
	auto cardinality = available_cardinality;
	switch (limit.limit_val.Type()) {
	case LimitNodeType::CONSTANT_VALUE:
		cardinality = MinValue(cardinality, limit.limit_val.GetConstantValue());
		break;
	case LimitNodeType::CONSTANT_PERCENTAGE: {
		auto percentage = limit.limit_val.GetConstantPercentage();
		auto percentage_cardinality = child_stats.cardinality;
		if (percentage < 100) {
			percentage_cardinality = LossyNumericCast<idx_t>(double(child_stats.cardinality) * percentage / 100.0);
		}
		cardinality = MinValue(cardinality, percentage_cardinality);
		break;
	}
	default:
		break;
	}
	ApplyCardinality(*result, cardinality);
	return result;
}

RelationStats RelationStatisticsHelper::ExtractEmptyResultStats(LogicalEmptyResult &empty) {
	RelationStats result;
	result.cardinality = 0;
	result.table_name = Identifier(empty.GetName());
	result.stats_initialized = true;
	for (auto &binding : empty.GetColumnBindings()) {
		result.columns.emplace_back(binding, DistinctCount(0, DistinctCountSource::CARDINALITY),
		                            Identifier("empty_result_column"));
	}
	result.Verify(empty.GetColumnBindings());
	return result;
}

optional<RelationStats> RelationStatisticsHelper::ProjectOutputStats(const RelationStats &stats, LogicalOperator &op) {
	if (!stats.stats_initialized) {
		return {};
	}
	RelationStats result;
	result.cardinality = stats.cardinality;
	result.filter_strength = stats.filter_strength;
	result.stats_initialized = true;
	result.table_name = stats.table_name;
	for (auto &binding : op.GetColumnBindings()) {
		auto column = stats.GetColumnStats(binding);
		if (!column) {
			return {};
		}
		result.columns.push_back(*column);
	}
	result.Verify(op.GetColumnBindings());
	return result;
}

optional<RelationStats> RelationStatisticsHelper::RebindOutputStats(const RelationStats &stats, LogicalOperator &op) {
	auto bindings = op.GetColumnBindings();
	if (!stats.stats_initialized || bindings.size() != stats.columns.size()) {
		return {};
	}
	auto result = stats;
	for (idx_t column_idx = 0; column_idx < bindings.size(); column_idx++) {
		result.columns[column_idx].binding = bindings[column_idx];
	}
	result.Verify(bindings);
	return result;
}

} // namespace duckdb
