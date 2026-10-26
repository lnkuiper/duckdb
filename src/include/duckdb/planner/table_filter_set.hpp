//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/planner/table_filter_set.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/planner/table_filter.hpp"
#include "duckdb/common/mutex.hpp"
#include "duckdb/common/reference_map.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/map.hpp"
#include "duckdb/common/atomic.hpp"

namespace duckdb {

class TableFilterSet {
public:
	void PushFilter(ProjectionIndex col_idx, unique_ptr<TableFilter> filter);
	void PushMultiColumnFilter(unique_ptr<TableFilter> filter);
	bool HasFilters() const;
	idx_t FilterCount() const;
	bool HasMultiColumnFilters() const;
	bool HasFilter(ProjectionIndex col_idx) const;
	TableFilter &GetFilterByColumnIndexMutable(ProjectionIndex col_idx);
	optional_ptr<TableFilter> TryGetFilterByColumnIndexMutable(ProjectionIndex col_idx);
	const TableFilter &GetFilterByColumnIndex(ProjectionIndex col_idx) const;
	optional_ptr<const TableFilter> TryGetFilterByColumnIndex(ProjectionIndex col_idx) const;
	void SetFilterByColumnIndex(ProjectionIndex col_idx, unique_ptr<TableFilter> filter);
	void RemoveFilterByColumnIndex(ProjectionIndex col_idx);
	void ClearFilters();
	const vector<unique_ptr<TableFilter>> &GetMultiColumnFilters() const;

	bool Equals(TableFilterSet &other);
	static bool Equals(TableFilterSet *left, TableFilterSet *right);

	unique_ptr<TableFilterSet> Copy() const;

	void Serialize(Serializer &serializer) const;
	static TableFilterSet Deserialize(Deserializer &deserializer);

	map<ProjectionIndex, unique_ptr<TableFilter>> GetTableFiltersForSerialization(Serializer &serializer) const;
	map<ProjectionIndex, unique_ptr<TableFilter>> &GetTableFiltersForDeserialization(Deserializer &deserializer);

public:
	class TableFilterIteratorEntry {
	public:
		explicit TableFilterIteratorEntry(map<ProjectionIndex, unique_ptr<TableFilter>>::iterator);

		ProjectionIndex GetIndex() const;
		TableFilter &Filter();
		const TableFilter &Filter() const;
		unique_ptr<TableFilter> TakeFilter();

	public:
		map<ProjectionIndex, unique_ptr<TableFilter>>::iterator iterator;
	};

	class ConstTableFilterIteratorEntry {
	public:
		explicit ConstTableFilterIteratorEntry(map<ProjectionIndex, unique_ptr<TableFilter>>::const_iterator);

		ProjectionIndex GetIndex() const;
		const TableFilter &Filter() const;

	public:
		map<ProjectionIndex, unique_ptr<TableFilter>>::const_iterator iterator;
	};

	// iterator
	template <class T>
	class TableFilterIterator {
	public:
		explicit TableFilterIterator(T entry_p) : entry(std::move(entry_p)) {
		}

	public:
		TableFilterIterator &operator++() {
			++entry.iterator;
			return *this;
		}
		bool operator!=(const TableFilterIterator &other) const {
			return entry.iterator != other.entry.iterator;
		}
		T &operator*() {
			return entry;
		}
		const T &operator*() const {
			return entry;
		}

	private:
		T entry;
	};

	TableFilterIterator<TableFilterIteratorEntry> begin() { // NOLINT: match stl API
		return TableFilterIterator<TableFilterIteratorEntry>(TableFilterIteratorEntry(filters.begin()));
	}
	TableFilterIterator<TableFilterIteratorEntry> end() { // NOLINT: match stl API
		return TableFilterIterator<TableFilterIteratorEntry>(TableFilterIteratorEntry(filters.end()));
	}
	TableFilterIterator<ConstTableFilterIteratorEntry> begin() const { // NOLINT: match stl API
		return TableFilterIterator<ConstTableFilterIteratorEntry>(ConstTableFilterIteratorEntry(filters.begin()));
	}
	TableFilterIterator<ConstTableFilterIteratorEntry> end() const { // NOLINT: match stl API
		return TableFilterIterator<ConstTableFilterIteratorEntry>(ConstTableFilterIteratorEntry(filters.end()));
	}

private:
	map<ProjectionIndex, unique_ptr<TableFilter>> filters;
	vector<unique_ptr<TableFilter>> multi_column_filters;
};

class DynamicTableFilterGroup;

class DynamicTableFilterSet {
public:
	explicit DynamicTableFilterSet(bool require_registration = false);

	void RegisterProducer(const PhysicalOperator &op);
	//! Reset a producer contribution and its availability for a new execution.
	void ResetFilters(const PhysicalOperator &op, bool enabled);
	void ClearFilters(const PhysicalOperator &op);
	void PushFilter(const PhysicalOperator &op, ProjectionIndex column_index, unique_ptr<TableFilter> filter);

	//! Groups are attached during planning; registered-producer leaves cannot contain groups.
	void AddGroup(shared_ptr<DynamicTableFilterGroup> group);
	vector<shared_ptr<DynamicTableFilterSet>> GetActiveGroupFilters() const;
	bool HasFilters() const;
	unique_ptr<TableFilterSet> GetFinalTableFilters(const PhysicalTableScan &scan,
	                                                optional_ptr<TableFilterSet> existing_filters) const;

private:
	friend class DynamicTableFilterGroup;

	bool HasRegisteredProducers() const;
	void ResetRegisteredProducers();
	unique_ptr<TableFilterSet> CopyFilters(optional_ptr<TableFilterSet> existing_filters = nullptr) const;
	vector<shared_ptr<DynamicTableFilterGroup>> GetGroups() const;

	struct FilterEntry {
		bool enabled = true;
		unique_ptr<TableFilterSet> filters;
	};

	const bool require_registration;
	mutable mutex lock;
	reference_map_t<const PhysicalOperator, FilterEntry> filters;
	vector<shared_ptr<DynamicTableFilterGroup>> groups;
};

//! A bounded disjunction of consumer conjunctions, shared by a scan and its CTE owner.
class DynamicTableFilterGroup {
public:
	explicit DynamicTableFilterGroup(vector<unique_ptr<TableFilterSet>> static_filters);

	idx_t ConsumerCount() const;
	shared_ptr<DynamicTableFilterSet> GetConsumerFilters(idx_t consumer_idx) const;
	bool HasPotentialFilters() const;
	bool HasFilters() const;
	bool IsActive() const;
	void SetActive(bool enabled);
	//! Reset all registered producers, including those absent from the next schedule.
	void Reset();
	void PushFilterSnapshot(TableFilterSet &result) const;

private:
	struct ConsumerFilters {
		explicit ConsumerFilters(unique_ptr<TableFilterSet> static_filters);

		unique_ptr<TableFilterSet> static_filters;
		shared_ptr<DynamicTableFilterSet> runtime_filters;
	};

	vector<ConsumerFilters> consumers;
	atomic<bool> active {false};
};

} // namespace duckdb
