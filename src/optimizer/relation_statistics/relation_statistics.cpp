#include "duckdb/optimizer/relation_statistics/relation_statistics.hpp"

namespace duckdb {

DistinctCount::DistinctCount(idx_t distinct_count, DistinctCountSource source)
    : distinct_count(distinct_count), source(source) {
}

void CurrentDomainEvidence::TightenFilterDomainBound(idx_t bound) {
	if (!filter_domain_bound.IsValid() || bound < filter_domain_bound.GetIndex()) {
		filter_domain_bound = optional_idx(bound);
	}
}

void CurrentDomainEvidence::Invalidate() {
	filter_domain_bound = optional_idx();
	is_unique = false;
}

RelationColumnStats::RelationColumnStats(ColumnBinding binding, DistinctCount domain, Identifier name)
    : RelationColumnStats(binding, domain, domain, std::move(name)) {
}

RelationColumnStats::RelationColumnStats(ColumnBinding binding, DistinctCount total_domain,
                                         DistinctCount current_domain, Identifier name,
                                         CurrentDomainEvidence current_domain_evidence)
    : binding(binding), total_domain(total_domain), current_domain(current_domain),
      current_domain_evidence(current_domain_evidence), name(std::move(name)) {
}

optional_idx RelationColumnStats::GetSemiAntiJoinDomainSize() const {
	optional_idx result;
	if (current_domain_evidence.is_unique) {
		result = optional_idx(current_domain.distinct_count);
	}
	if (current_domain_evidence.filter_domain_bound.IsValid()) {
		auto filter_domain_bound = current_domain_evidence.filter_domain_bound.GetIndex();
		result =
		    optional_idx(result.IsValid() ? MinValue(result.GetIndex(), filter_domain_bound) : filter_domain_bound);
	}
	return result;
}

RelationStats::RelationStats() : cardinality(1), filter_strength(1), stats_initialized(false) {
}

optional_idx RelationStats::FindColumn(ColumnBinding binding) const {
	for (idx_t column_idx = 0; column_idx < columns.size(); column_idx++) {
		if (columns[column_idx].binding == binding) {
			return column_idx;
		}
	}
	return {};
}

optional_ptr<const RelationColumnStats> RelationStats::GetColumnStats(ColumnBinding binding) const {
	auto column_idx = FindColumn(binding);
	if (!column_idx.IsValid()) {
		return nullptr;
	}
	return columns[column_idx.GetIndex()];
}

bool RelationStats::MatchesBindings(const vector<ColumnBinding> &bindings) const {
	if (columns.size() != bindings.size()) {
		return false;
	}
	for (idx_t column_idx = 0; column_idx < columns.size(); column_idx++) {
		if (columns[column_idx].binding != bindings[column_idx]) {
			return false;
		}
	}
	return true;
}

void RelationStats::Verify(const vector<ColumnBinding> &bindings) const {
	D_ASSERT(!stats_initialized || MatchesBindings(bindings));
	if (!stats_initialized) {
		return;
	}
	D_ASSERT(filter_strength >= 0 && filter_strength <= 1);
	for (auto &column : columns) {
		(void)column;
		D_ASSERT(column.current_domain.distinct_count <= column.total_domain.distinct_count);
		D_ASSERT(column.current_domain.distinct_count <= cardinality);
	}
}

} // namespace duckdb
