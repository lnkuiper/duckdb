//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/optimizer/relation_statistics/relation_statistics.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/optional_idx.hpp"
#include "duckdb/common/optional_ptr.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/planner/column_binding.hpp"

namespace duckdb {

enum class DistinctCountSource : uint8_t { CARDINALITY, MIN_MAX, HLL, EXACT };

enum class CurrentDomainProvenance : uint8_t { UNKNOWN, BASE, MODELED };

struct DistinctCount {
public:
	DistinctCount(idx_t distinct_count, DistinctCountSource source);

public:
	idx_t distinct_count;
	DistinctCountSource source;
};

struct CurrentDomainInfo {
public:
	explicit CurrentDomainInfo(CurrentDomainProvenance provenance = CurrentDomainProvenance::UNKNOWN,
	                           optional_idx direct_bound = {}, bool is_unique = false);

public:
	bool IsEligibleForSemiAnti() const;
	void UpdateDirectBound(idx_t bound);
	void MarkModeled();
	void InvalidateStructuralEvidence();

public:
	CurrentDomainProvenance provenance;
	optional_idx direct_bound;
	bool is_unique;
};

struct RelationColumnStats {
public:
	RelationColumnStats(ColumnBinding binding, DistinctCount domain, Identifier name);
	RelationColumnStats(ColumnBinding binding, DistinctCount total_domain, DistinctCount current_domain,
	                    Identifier name, CurrentDomainInfo current_domain_info = CurrentDomainInfo());
	optional_idx GetSemiAntiCurrentDomain() const;

public:
	ColumnBinding binding;
	//! The estimated size of the complete value domain. Equality joins use this as their denominator.
	DistinctCount total_domain;
	//! The estimated number of values from the total domain that are currently present.
	DistinctCount current_domain;
	//! Describes how current_domain was derived and which structural bounds survive later row filtering.
	CurrentDomainInfo current_domain_info;
	Identifier name;
};

struct RelationStats {
public:
	RelationStats();

public:
	optional_idx FindColumn(ColumnBinding binding) const;
	optional_ptr<const RelationColumnStats> GetColumnStats(ColumnBinding binding) const;
	bool MatchesBindings(const vector<ColumnBinding> &bindings) const;
	void Verify(const vector<ColumnBinding> &bindings) const;

public:
	vector<RelationColumnStats> columns;
	idx_t cardinality;
	double filter_strength = 1;
	bool stats_initialized = false;
	Identifier table_name;
};

} // namespace duckdb
