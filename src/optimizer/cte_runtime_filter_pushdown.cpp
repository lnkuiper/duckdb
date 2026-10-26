#include "duckdb/optimizer/join_filter_pushdown_optimizer.hpp"

#include "duckdb/execution/operator/join/join_filter_pushdown.hpp"
#include "duckdb/planner/expression/bound_aggregate_expression.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_reference_expression.hpp"
#include "duckdb/planner/expression/expression_barrier.hpp"
#include "duckdb/planner/expression_iterator.hpp"
#include "duckdb/planner/filter/expression_filter.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_cteref.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_materialized_cte.hpp"

namespace duckdb {

using cte_filter_references_t = unordered_map<TableIndex, vector<reference<LogicalCTERef>>>;

static void CollectRuntimeCTEReferences(LogicalOperator &op, vector<reference<LogicalMaterializedCTE>> &ctes,
                                        cte_filter_references_t &references) {
	if (op.type == LogicalOperatorType::LOGICAL_MATERIALIZED_CTE) {
		ctes.push_back(op.Cast<LogicalMaterializedCTE>());
	} else if (op.type == LogicalOperatorType::LOGICAL_CTE_REF) {
		auto &ref = op.Cast<LogicalCTERef>();
		references[ref.cte_index].push_back(ref);
	}
	for (auto &child : op.children) {
		CollectRuntimeCTEReferences(*child, ctes, references);
	}
}

static bool ContainsRecursiveCTE(const LogicalOperator &op) {
	if (op.type == LogicalOperatorType::LOGICAL_RECURSIVE_CTE) {
		return true;
	}
	for (auto &child : op.children) {
		if (ContainsRecursiveCTE(*child)) {
			return true;
		}
	}
	return false;
}

void JoinFilterPushdownOptimizer::Optimize(LogicalOperator &op) {
	eligible_ctes.clear();
	if (!ContainsRecursiveCTE(op)) {
		vector<reference<LogicalMaterializedCTE>> ctes;
		cte_filter_references_t references;
		CollectRuntimeCTEReferences(op, ctes, references);
		if (!ctes.empty()) {
			op.ResolveOperatorTypes();
		}
		for (auto &cte_ref : ctes) {
			auto &cte = cte_ref.get();
			cte.filter_pushdown.reset();
			auto entry = references.find(cte.table_index);
			if (entry == references.end() || entry->second.size() < 2) {
				continue;
			}
			auto info = make_uniq<CTEFilterPushdownInfo>();
			info->producer_bindings = cte.children[0]->GetColumnBindings();
			info->producer_types = cte.children[0]->types;
			bool valid = true;
			for (auto &consumer_ref : entry->second) {
				auto &consumer = consumer_ref.get();
				valid &= !consumer.is_recurring && consumer.correlated_columns == 0 &&
				         consumer.chunk_types == info->producer_types;
				valid &= info->consumers.emplace(consumer.table_index, consumer.chunk_types).second;
			}
			if (valid) {
				eligible_ctes.insert(cte.table_index);
				cte.filter_pushdown = std::move(info);
			}
		}
	}
	VisitOperator(op);
}

static bool MatchesRuntimeCTEReferences(LogicalMaterializedCTE &cte,
                                        const vector<reference<LogicalCTERef>> &references) {
	auto &info = *cte.filter_pushdown;
	if (info.producer_bindings != cte.children[0]->GetColumnBindings() ||
	    info.producer_types != cte.children[0]->types || info.consumers.size() != references.size()) {
		return false;
	}
	unordered_set<TableIndex> seen;
	for (auto &ref : references) {
		auto &consumer = ref.get();
		auto entry = info.consumers.find(consumer.table_index);
		if (entry == info.consumers.end() || entry->second != consumer.chunk_types || consumer.is_recurring ||
		    consumer.correlated_columns != 0 || !seen.insert(consumer.table_index).second) {
			return false;
		}
	}
	return true;
}

struct CTERuntimeFilterContribution {
	CTERuntimeFilterContribution(LogicalComparisonJoin &join, idx_t consumer_idx,
	                             vector<JoinFilterPushdownColumn> columns)
	    : join(join), consumer_idx(consumer_idx), columns(std::move(columns)) {
	}
	LogicalComparisonJoin &join;
	idx_t consumer_idx;
	vector<JoinFilterPushdownColumn> columns;
};

struct CTENativeFilterTarget {
	CTENativeFilterTarget(LogicalGet &get, idx_t consumer_count) : get(get), static_filters(consumer_count) {
	}
	LogicalGet &get;
	vector<unique_ptr<TableFilterSet>> static_filters;
	vector<CTERuntimeFilterContribution> contributions;
};

struct CTERuntimeFilterContext {
	explicit CTERuntimeFilterContext(LogicalMaterializedCTE &cte) : cte(cte) {
		for (auto &consumer : cte.filter_pushdown->consumers) {
			consumer_indexes.emplace(consumer.first, consumer_indexes.size());
		}
	}

	CTENativeFilterTarget &GetTarget(LogicalGet &get) {
		auto &target = targets[get];
		if (!target) {
			target = make_uniq<CTENativeFilterTarget>(get, consumer_indexes.size());
		}
		return *target;
	}

	LogicalMaterializedCTE &cte;
	unordered_map<TableIndex, idx_t> consumer_indexes;
	reference_map_t<LogicalGet, unique_ptr<CTENativeFilterTarget>> targets;
};

static void VisitRuntimeFilterPlan(LogicalOperator &op, const std::function<void(LogicalOperator &)> &callback) {
	callback(op);
	for (auto &child : op.children) {
		VisitRuntimeFilterPlan(*child, callback);
	}
}

static vector<PushdownFilterTarget> FindCTEProducerTargets(CTERuntimeFilterContext &context,
                                                           vector<JoinFilterPushdownColumn> columns) {
	auto &bindings = context.cte.filter_pushdown->producer_bindings;
	for (auto &column : columns) {
		auto index = column.probe_column_index.column_index.GetIndex();
		if (index >= bindings.size()) {
			return {};
		}
		column.probe_column_index = bindings[index];
	}
	vector<PushdownFilterTarget> result;
	JoinFilterPushdownOptimizer::GetPushdownFilterTargets(*context.cte.children[0], std::move(columns), result);
	return result;
}

static bool GetStaticFilterColumns(const Expression &expression, vector<ColumnBinding> &bindings,
                                   vector<LogicalType> &types, vector<JoinFilterPushdownColumn> &columns) {
	if (expression.CanThrow() || expression.IsVolatile() || ExpressionBarrier::Contains(expression)) {
		return false;
	}
	bool valid = true;
	ExpressionIterator::VisitExpression<BoundColumnRefExpression>(expression, [&](const BoundColumnRefExpression &ref) {
		if (ref.Depth() != 0) {
			valid = false;
			return;
		}
		for (auto binding : bindings) {
			if (binding == ref.Binding()) {
				return;
			}
		}
		JoinFilterPushdownColumn column;
		if (!JoinFilterPushdownUtil::PushdownJoinFilterExpression(ref, column)) {
			valid = false;
			return;
		}
		column.join_filter_idx = bindings.size();
		bindings.push_back(ref.Binding());
		types.push_back(ref.GetReturnType());
		columns.push_back(std::move(column));
	});
	return valid && !columns.empty();
}

static void AddStaticCTEFilter(CTENativeFilterTarget &target, idx_t consumer_idx, const Expression &expression,
                               const vector<ColumnBinding> &bindings, const vector<LogicalType> &types,
                               const vector<JoinFilterPushdownColumn> &columns) {
	vector<optional_idx> indexes(bindings.size());
	vector<ProjectionIndex> inputs;
	for (auto &column : columns) {
		auto origin = column.join_filter_idx;
		if (origin >= bindings.size() || indexes[origin].IsValid() || column.storage_type != types[origin] ||
		    !column.runtime_filter_casts.empty() || column.mode != JoinFilterPushdownMode::RECONSTRUCT_EXPRESSION) {
			return;
		}
		auto scan_column = column.probe_column_index.column_index;
		idx_t index = 0;
		for (; index < inputs.size(); index++) {
			if (inputs[index] == scan_column) {
				break;
			}
		}
		if (index == inputs.size()) {
			inputs.push_back(scan_column);
		}
		indexes[origin] = index;
	}
	for (auto &index : indexes) {
		if (!index.IsValid()) {
			return;
		}
	}
	auto copy = expression.Copy();
	ExpressionIterator::VisitExpressionMutable<BoundColumnRefExpression>(
	    copy, [&](BoundColumnRefExpression &ref, unique_ptr<Expression> &owned) {
		    for (idx_t i = 0; i < bindings.size(); i++) {
			    if (bindings[i] == ref.Binding()) {
				    owned = make_uniq<BoundReferenceExpression>(ref.GetReturnType(), indexes[i].GetIndex());
				    return;
			    }
		    }
		    throw InternalException("CTE static filter contains an unmapped binding");
	    });
	auto &filters = target.static_filters[consumer_idx];
	if (!filters) {
		filters = make_uniq<TableFilterSet>();
	}
	if (inputs.size() == 1) {
		filters->PushFilter(inputs[0], make_uniq<ExpressionFilter>(std::move(copy)));
	} else {
		filters->PushMultiColumnFilter(make_uniq<ExpressionFilter>(std::move(copy), std::move(inputs)));
	}
}

void JoinFilterPushdownOptimizer::FinalizeCTEFilters(LogicalOperator &op) {
	if (ContainsRecursiveCTE(op)) {
		return;
	}
	vector<reference<LogicalMaterializedCTE>> ctes;
	cte_filter_references_t references;
	CollectRuntimeCTEReferences(op, ctes, references);
	if (ctes.empty()) {
		return;
	}
	unordered_map<TableIndex, unique_ptr<CTERuntimeFilterContext>> contexts;
	for (auto &cte_ref : ctes) {
		auto &cte = cte_ref.get();
		// Table-parameter subqueries can enter physical planning again.
		if (!cte.filter_pushdown || !cte.filter_pushdown->targets.empty()) {
			continue;
		}
		auto entry = references.find(cte.table_index);
		if (entry == references.end() || !MatchesRuntimeCTEReferences(cte, entry->second)) {
			cte.filter_pushdown.reset();
			continue;
		}
		contexts.emplace(cte.table_index, make_uniq<CTERuntimeFilterContext>(cte));
	}
	VisitRuntimeFilterPlan(op, [&](LogicalOperator &current) {
		if (current.type != LogicalOperatorType::LOGICAL_COMPARISON_JOIN) {
			return;
		}
		auto &join = current.Cast<LogicalComparisonJoin>();
		if (!join.filter_pushdown || !JoinFilterPushdownUtil::JoinTypeIsSupported(join.join_type)) {
			return;
		}
		auto &info = *join.filter_pushdown;
		if (info.min_max_aggregates.size() != 2 * info.join_condition.size()) {
			return;
		}
		vector<JoinFilterPushdownColumn> columns;
		for (idx_t i = 0; i < info.join_condition.size(); i++) {
			auto condition_idx = info.join_condition[i];
			JoinFilterPushdownColumn column;
			if (!GetJoinFilterColumn(join, condition_idx, column)) {
				continue;
			}
			auto &condition = join.conditions[condition_idx];
			if (condition.GetComparisonType() != ExpressionType::COMPARE_EQUAL) {
				continue;
			}
			bool valid = true;
			for (idx_t offset = 0; offset < 2; offset++) {
				auto &aggregate = info.min_max_aggregates[2 * i + offset]->Cast<BoundAggregateExpression>();
				valid &= aggregate.GetChildren().size() == 1 && aggregate.GetChildren()[0]->Equals(condition.GetRHS());
			}
			if (valid) {
				column.join_filter_idx = i;
				columns.push_back(std::move(column));
			}
		}
		vector<PushdownFilterTarget> targets;
		GetPushdownFilterTargets(*join.children[0], std::move(columns), targets, EqualityPropagation::DISABLED,
		                         CTEFilterTargetMode::INCLUDE);
		for (auto &target : targets) {
			if (!target.IsCTEReference()) {
				continue;
			}
			auto &ref = target.GetCTEReference();
			auto entry = contexts.find(ref.cte_index);
			if (entry == contexts.end()) {
				continue;
			}
			auto &context = *entry->second;
			auto consumer_idx = context.consumer_indexes.at(ref.table_index);
			for (auto &native : FindCTEProducerTargets(context, std::move(target.columns))) {
				context.GetTarget(native.GetScan())
				    .contributions.emplace_back(join, consumer_idx, std::move(native.columns));
			}
		}
	});
	VisitRuntimeFilterPlan(op, [&](LogicalOperator &current) {
		if (current.type != LogicalOperatorType::LOGICAL_FILTER) {
			return;
		}
		for (auto &expression : current.expressions) {
			vector<ColumnBinding> bindings;
			vector<LogicalType> types;
			vector<JoinFilterPushdownColumn> columns;
			if (!GetStaticFilterColumns(*expression, bindings, types, columns)) {
				continue;
			}
			vector<PushdownFilterTarget> targets;
			GetPushdownFilterTargets(current, std::move(columns), targets, EqualityPropagation::DISABLED,
			                         CTEFilterTargetMode::INCLUDE);
			for (auto &target : targets) {
				if (!target.IsCTEReference() || target.columns.size() != bindings.size()) {
					continue;
				}
				auto &ref = target.GetCTEReference();
				auto entry = contexts.find(ref.cte_index);
				if (entry == contexts.end()) {
					continue;
				}
				auto &context = *entry->second;
				auto consumer_idx = context.consumer_indexes.at(ref.table_index);
				for (auto &native : FindCTEProducerTargets(context, std::move(target.columns))) {
					AddStaticCTEFilter(context.GetTarget(native.GetScan()), consumer_idx, *expression, bindings, types,
					                   native.columns);
				}
			}
		}
	});
	for (auto &entry : contexts) {
		auto &context = *entry.second;
		auto &info = *context.cte.filter_pushdown;
		for (auto &native_entry : context.targets) {
			auto &target = *native_entry.second;
			if (target.contributions.empty()) {
				continue;
			}
			vector<bool> restricted(context.consumer_indexes.size(), false);
			for (idx_t i = 0; i < target.static_filters.size(); i++) {
				auto &filters = target.static_filters[i];
				restricted[i] = filters && (filters->HasFilters() || filters->HasMultiColumnFilters());
			}
			for (auto &contribution : target.contributions) {
				restricted[contribution.consumer_idx] = true;
			}
			bool useful = true;
			for (auto value : restricted) {
				useful &= value;
			}
			if (!useful) {
				continue;
			}
			auto group = make_shared_ptr<DynamicTableFilterGroup>(std::move(target.static_filters));
			if (!target.get.dynamic_filters) {
				target.get.dynamic_filters = make_shared_ptr<DynamicTableFilterSet>();
			}
			info.targets.push_back({target.get.dynamic_filters, group});
			for (auto &contribution : target.contributions) {
				JoinFilterPushdownFilter filter;
				filter.dynamic_filters = group->GetConsumerFilters(contribution.consumer_idx);
				filter.columns = std::move(contribution.columns);
				filter.estimated_scan_cardinality =
				    target.get.has_estimated_cardinality ? target.get.estimated_cardinality : idx_t(0);
				contribution.join.filter_pushdown->probe_info.push_back(std::move(filter));
			}
		}
		if (info.targets.empty()) {
			context.cte.filter_pushdown.reset();
		}
	}
	VisitRuntimeFilterPlan(op, [&](LogicalOperator &current) {
		if (current.type != LogicalOperatorType::LOGICAL_COMPARISON_JOIN) {
			return;
		}
		auto &join = current.Cast<LogicalComparisonJoin>();
		if (join.filter_pushdown && join.filter_pushdown->probe_info.empty() &&
		    !ComputeAggregatesForPerfectHash(join, join.filter_pushdown->join_condition.size())) {
			join.filter_pushdown.reset();
		}
	});
}

} // namespace duckdb
