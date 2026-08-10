#include "duckdb/optimizer/relation_statistics/relation_statistics.hpp"

namespace duckdb {

DistinctCount::DistinctCount(idx_t distinct_count, DistinctCountSource source)
    : distinct_count(distinct_count), source(source) {
}

CurrentDomainInfo::CurrentDomainInfo(CurrentDomainProvenance provenance, optional_idx direct_bound, bool is_unique)
    : provenance(provenance), direct_bound(direct_bound), is_unique(is_unique) {
}

bool CurrentDomainInfo::IsEligibleForSemiAnti() const {
	return direct_bound.IsValid() || is_unique;
}

void CurrentDomainInfo::UpdateDirectBound(idx_t bound) {
	if (!direct_bound.IsValid() || bound < direct_bound.GetIndex()) {
		direct_bound = optional_idx(bound);
	}
}

void CurrentDomainInfo::MarkModeled() {
	provenance = CurrentDomainProvenance::MODELED;
}

void CurrentDomainInfo::InvalidateStructuralEvidence() {
	MarkModeled();
	direct_bound = optional_idx();
	is_unique = false;
}

RelationColumnStats::RelationColumnStats(ColumnBinding binding, DistinctCount domain, Identifier name)
    : RelationColumnStats(binding, domain, domain, std::move(name)) {
}

RelationColumnStats::RelationColumnStats(ColumnBinding binding, DistinctCount total_domain,
                                         DistinctCount current_domain, Identifier name,
                                         CurrentDomainInfo current_domain_info)
    : binding(binding), total_domain(total_domain), current_domain(current_domain),
      current_domain_info(current_domain_info), name(std::move(name)) {
}

optional_idx RelationColumnStats::GetSemiAntiCurrentDomain() const {
	optional_idx result;
	if (current_domain_info.is_unique) {
		result = optional_idx(current_domain.distinct_count);
	}
	if (current_domain_info.direct_bound.IsValid()) {
		auto direct_bound = current_domain_info.direct_bound.GetIndex();
		result = optional_idx(result.IsValid() ? MinValue(result.GetIndex(), direct_bound) : direct_bound);
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
