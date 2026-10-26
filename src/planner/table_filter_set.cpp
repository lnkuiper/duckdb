#include "duckdb/planner/table_filter_set.hpp"
#include "duckdb/planner/filter/expression_filter.hpp"
#include "duckdb/planner/filter/conjunction_filter.hpp"
#include "duckdb/planner/filter/constant_filter.hpp"
#include "duckdb/planner/filter/dynamic_filter.hpp"
#include "duckdb/planner/filter/in_filter.hpp"
#include "duckdb/planner/filter/null_filter.hpp"
#include "duckdb/planner/filter/optional_filter.hpp"
#include "duckdb/planner/filter/struct_filter.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/expression/bound_operator_expression.hpp"
#include "duckdb/planner/expression/bound_reference_expression.hpp"
#include "duckdb/planner/expression/bound_conjunction_expression.hpp"
#include "duckdb/planner/expression_iterator.hpp"
#include "duckdb/planner/filter/table_filter_functions.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/serializer/serializer.hpp"
#include "duckdb/common/serializer/deserializer.hpp"
#include "duckdb/function/scalar/struct_utils.hpp"

namespace duckdb {

static bool ContainsInternalTableFilterFunction(const Expression &expr) {
	if (expr.GetExpressionClass() == ExpressionClass::BOUND_FUNCTION) {
		auto &func = expr.Cast<BoundFunctionExpression>();
		if (TableFilterFunctions::IsTableFilterFunction(func.Function())) {
			return true;
		}
		auto optional_child = ExpressionFilter::GetOptionalFilterChild(func);
		if (optional_child && ContainsInternalTableFilterFunction(*optional_child)) {
			return true;
		}
	}
	bool found = false;
	ExpressionIterator::EnumerateChildren(expr, [&](const Expression &child) {
		if (found) {
			return;
		}
		found = ContainsInternalTableFilterFunction(child);
	});
	return found;
}

static unique_ptr<TableFilter> SerializeExpressionToLegacyFilter(const Expression &expr);

static ExpressionType FlipComparisonType(ExpressionType type) {
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

static bool IsSupportedConstantComparison(ExpressionType type) {
	switch (type) {
	case ExpressionType::COMPARE_EQUAL:
	case ExpressionType::COMPARE_NOTEQUAL:
	case ExpressionType::COMPARE_LESSTHAN:
	case ExpressionType::COMPARE_LESSTHANOREQUALTO:
	case ExpressionType::COMPARE_GREATERTHAN:
	case ExpressionType::COMPARE_GREATERTHANOREQUALTO:
		return true;
	default:
		return false;
	}
}

static bool TryExtractLegacySubject(const Expression &expr, vector<StructExtractPathEntry> &struct_path) {
	return ExpressionFilter::IsSimpleFilterColumnRef(PeelStructExtractPath(expr, struct_path));
}

static unique_ptr<TableFilter> WrapStructFilterPath(unique_ptr<TableFilter> filter,
                                                    const vector<StructExtractPathEntry> &struct_path) {
	for (auto it = struct_path.rbegin(); it != struct_path.rend(); ++it) {
		filter = make_uniq<LegacyStructFilter>(it->child_idx, it->child_name, std::move(filter));
	}
	return filter;
}

static void NormalizeLegacyExpression(unique_ptr<Expression> &expr) {
	ExpressionIterator::VisitExpressionMutable<BoundColumnRefExpression>(
	    expr, [](BoundColumnRefExpression &col_ref, unique_ptr<Expression> &owned_expr) {
		    owned_expr = make_uniq<BoundReferenceExpression>(col_ref.GetAlias(), col_ref.GetReturnType(), 0ULL);
	    });
	ExpressionIterator::VisitExpressionMutable<BoundReferenceExpression>(
	    expr, [](BoundReferenceExpression &ref, unique_ptr<Expression> &owned_expr) { ref.IndexMutable() = 0; });
}

static unique_ptr<TableFilter> TrySerializeComparisonToLegacyFilter(const BoundFunctionExpression &comparison) {
	auto comparison_type = comparison.GetExpressionType();
	auto &left = BoundComparisonExpression::Left(comparison);
	auto &right = BoundComparisonExpression::Right(comparison);
	bool lhs_constant = left.GetExpressionClass() == ExpressionClass::BOUND_CONSTANT;
	bool rhs_constant = right.GetExpressionClass() == ExpressionClass::BOUND_CONSTANT;
	if (!lhs_constant && !rhs_constant) {
		return nullptr;
	}
	auto &subject = rhs_constant ? left : right;
	auto &constant_expr = rhs_constant ? right : left;
	const auto &constant = constant_expr.Cast<BoundConstantExpression>().GetValue();
	if (!rhs_constant) {
		comparison_type = FlipComparisonType(comparison_type);
	}

	vector<StructExtractPathEntry> struct_path;
	if (!TryExtractLegacySubject(subject, struct_path)) {
		return nullptr;
	}
	if (constant.IsNull()) {
		switch (comparison_type) {
		case ExpressionType::COMPARE_DISTINCT_FROM:
			return WrapStructFilterPath(make_uniq<LegacyIsNotNullFilter>(), struct_path);
		case ExpressionType::COMPARE_NOT_DISTINCT_FROM:
			return WrapStructFilterPath(make_uniq<LegacyIsNullFilter>(), struct_path);
		default:
			return nullptr;
		}
	}
	if (!IsSupportedConstantComparison(comparison_type)) {
		return nullptr;
	}
	return WrapStructFilterPath(make_uniq<LegacyConstantFilter>(comparison_type, constant), struct_path);
}

static unique_ptr<TableFilter> TrySerializeOperatorToLegacyFilter(const BoundOperatorExpression &op) {
	switch (op.GetExpressionType()) {
	case ExpressionType::OPERATOR_IS_NULL:
	case ExpressionType::OPERATOR_IS_NOT_NULL: {
		if (op.GetChildren().size() != 1) {
			return nullptr;
		}
		vector<StructExtractPathEntry> struct_path;
		if (!TryExtractLegacySubject(*op.GetChildren()[0], struct_path)) {
			return nullptr;
		}
		if (op.GetExpressionType() == ExpressionType::OPERATOR_IS_NULL) {
			return WrapStructFilterPath(make_uniq<LegacyIsNullFilter>(), struct_path);
		}
		return WrapStructFilterPath(make_uniq<LegacyIsNotNullFilter>(), struct_path);
	}
	case ExpressionType::COMPARE_IN: {
		if (op.GetChildren().empty()) {
			return nullptr;
		}
		vector<StructExtractPathEntry> struct_path;
		if (!TryExtractLegacySubject(*op.GetChildren()[0], struct_path)) {
			return nullptr;
		}
		vector<Value> values;
		values.reserve(op.GetChildren().size() - 1);
		for (idx_t i = 1; i < op.GetChildren().size(); i++) {
			if (op.GetChildren()[i]->GetExpressionClass() != ExpressionClass::BOUND_CONSTANT) {
				return nullptr;
			}
			auto value = op.GetChildren()[i]->Cast<BoundConstantExpression>().GetValue();
			if (value.IsNull()) {
				return nullptr;
			}
			values.push_back(std::move(value));
		}
		if (values.empty()) {
			return nullptr;
		}
		return WrapStructFilterPath(make_uniq<LegacyInFilter>(std::move(values)), struct_path);
	}
	default:
		return nullptr;
	}
}
static unique_ptr<TableFilter> SerializeOptionalChild(const optional_ptr<const Expression> child_expr) {
	if (!child_expr) {
		return nullptr;
	}
	return SerializeExpressionToLegacyFilter(*child_expr);
}

static unique_ptr<TableFilter> SerializeInternalFunctionToLegacyFilter(const BoundFunctionExpression &func_expr) {
	auto &func_name = func_expr.Function().GetName();
	if (ExpressionFilter::IsRootOptionalExpression(func_expr)) {
		auto child_filter = SerializeOptionalChild(ExpressionFilter::GetOptionalFilterChild(func_expr));
		return make_uniq<LegacyOptionalFilter>(std::move(child_filter));
	}
	if (func_name == DynamicFilterScalarFun::NAME) {
		if (!func_expr.BindInfo()) {
			return make_uniq<LegacyDynamicFilter>();
		}
		if (func_expr.GetChildren().size() != 1 ||
		    func_expr.GetChildren()[0]->GetExpressionType() != ExpressionType::BOUND_REF) {
			// the dynamic filter is evaluated on top of a reconstructed expression (e.g. a cast chain
			// over the raw scan column) - the legacy format cannot represent this, so drop the filter
			// (the optional wrapper degrades to a no-op)
			return nullptr;
		}
		auto &data = func_expr.BindInfo()->Cast<DynamicFilterFunctionData>();
		return make_uniq<LegacyDynamicFilter>(data.filter_data);
	}
	if (func_name == BloomFilterScalarFun::NAME || func_name == PrefixRangeScalarFun::NAME) {
		return make_uniq<LegacyOptionalFilter>();
	}
	throw SerializationException("Unsupported internal tablefilter function %s during serialization", func_name);
}

static unique_ptr<TableFilter> SerializeConjunctionToLegacyFilter(const BoundConjunctionExpression &conjunction) {
	unique_ptr<LegacyConjunctionFilter> result;
	if (conjunction.GetExpressionType() == ExpressionType::CONJUNCTION_AND) {
		result = make_uniq<LegacyConjunctionAndFilter>();
	} else if (conjunction.GetExpressionType() == ExpressionType::CONJUNCTION_OR) {
		result = make_uniq<LegacyConjunctionOrFilter>();
	} else {
		throw SerializationException("Unsupported conjunction type %s during table-filter serialization",
		                             EnumUtil::ToString(conjunction.GetExpressionType()));
	}
	for (auto &child : conjunction.GetChildren()) {
		auto child_filter = SerializeExpressionToLegacyFilter(*child);
		if (!child_filter) {
			return nullptr;
		}
		result->child_filters.push_back(std::move(child_filter));
	}
	return std::move(result);
}

static unique_ptr<TableFilter> SerializeExpressionToLegacyFilter(const Expression &expr) {
	if (expr.GetExpressionClass() == ExpressionClass::BOUND_CONJUNCTION) {
		return SerializeConjunctionToLegacyFilter(expr.Cast<BoundConjunctionExpression>());
	}
	if (BoundComparisonExpression::IsComparison(expr.GetExpressionType())) {
		auto result = TrySerializeComparisonToLegacyFilter(expr.Cast<BoundFunctionExpression>());
		if (result) {
			return result;
		}
	}
	if (expr.GetExpressionClass() == ExpressionClass::BOUND_OPERATOR) {
		auto result = TrySerializeOperatorToLegacyFilter(expr.Cast<BoundOperatorExpression>());
		if (result) {
			return result;
		}
	}
	if (expr.GetExpressionClass() == ExpressionClass::BOUND_FUNCTION) {
		auto &func = expr.Cast<BoundFunctionExpression>();
		if (TableFilterFunctions::IsTableFilterFunction(func.Function())) {
			return SerializeInternalFunctionToLegacyFilter(func);
		}
	}
	if (ContainsInternalTableFilterFunction(expr)) {
		return nullptr;
	}
	auto normalized_expr = expr.Copy();
	NormalizeLegacyExpression(normalized_expr);
	return make_uniq<ExpressionFilter>(std::move(normalized_expr));
}

TableFilterSet::ConstTableFilterIteratorEntry::ConstTableFilterIteratorEntry(
    map<ProjectionIndex, unique_ptr<TableFilter>>::const_iterator it)
    : iterator(it) {
}

ProjectionIndex TableFilterSet::ConstTableFilterIteratorEntry::GetIndex() const {
	return iterator->first;
}

const TableFilter &TableFilterSet::ConstTableFilterIteratorEntry::Filter() const {
	return *iterator->second;
}

TableFilterSet::TableFilterIteratorEntry::TableFilterIteratorEntry(
    map<ProjectionIndex, unique_ptr<TableFilter>>::iterator it)
    : iterator(it) {
}

ProjectionIndex TableFilterSet::TableFilterIteratorEntry::GetIndex() const {
	return iterator->first;
}

TableFilter &TableFilterSet::TableFilterIteratorEntry::Filter() {
	return *iterator->second;
}

const TableFilter &TableFilterSet::TableFilterIteratorEntry::Filter() const {
	return *iterator->second;
}

unique_ptr<TableFilter> TableFilterSet::TableFilterIteratorEntry::TakeFilter() {
	return std::move(iterator->second);
}

bool TableFilterSet::HasFilters() const {
	return !filters.empty();
}
idx_t TableFilterSet::FilterCount() const {
	return filters.size();
}
bool TableFilterSet::HasMultiColumnFilters() const {
	return !multi_column_filters.empty();
}
bool TableFilterSet::HasFilter(ProjectionIndex col_idx) const {
	return filters.find(col_idx) != filters.end();
}

const TableFilter &TableFilterSet::GetFilterByColumnIndex(ProjectionIndex col_idx) const {
	auto filter = TryGetFilterByColumnIndex(col_idx);
	if (!filter) {
		throw InternalException("Table filter set does not have a filter for column idx %d", col_idx);
	}
	return *filter;
}

optional_ptr<const TableFilter> TableFilterSet::TryGetFilterByColumnIndex(ProjectionIndex col_idx) const {
	if (!col_idx.IsValid()) {
		throw InternalException("TableFilterSet::TryGetFilterByColumnIndex called with invalid column index");
	}
	auto entry = filters.find(col_idx);
	if (entry == filters.end()) {
		return nullptr;
	}
	return *entry->second;
}

TableFilter &TableFilterSet::GetFilterByColumnIndexMutable(ProjectionIndex col_idx) {
	auto filter = TryGetFilterByColumnIndexMutable(col_idx);
	if (!filter) {
		throw InternalException("Table filter set does not have a filter for column idx %d", col_idx);
	}
	return *filter;
}

optional_ptr<TableFilter> TableFilterSet::TryGetFilterByColumnIndexMutable(ProjectionIndex col_idx) {
	auto entry = filters.find(col_idx);
	if (entry == filters.end()) {
		return nullptr;
	}
	return *entry->second;
}

void TableFilterSet::RemoveFilterByColumnIndex(ProjectionIndex col_idx) {
	filters.erase(col_idx);
}

void TableFilterSet::SetFilterByColumnIndex(ProjectionIndex col_idx, unique_ptr<TableFilter> filter) {
	filters[col_idx] = std::move(filter);
}

void TableFilterSet::ClearFilters() {
	filters.clear();
	multi_column_filters.clear();
}

bool TableFilterSet::Equals(TableFilterSet &other) {
	if (filters.size() != other.filters.size() || multi_column_filters.size() != other.multi_column_filters.size()) {
		return false;
	}
	for (auto &entry : filters) {
		auto other_entry = other.filters.find(entry.first);
		if (other_entry == other.filters.end()) {
			return false;
		}
		if (!entry.second->Cast<ExpressionFilter>().Equals(other_entry->second->Cast<ExpressionFilter>())) {
			return false;
		}
	}
	for (idx_t filter_idx = 0; filter_idx < multi_column_filters.size(); ++filter_idx) {
		const auto &filter =
		    ExpressionFilter::GetExpressionFilter(*multi_column_filters[filter_idx], "TableFilterSet::Equals");
		const auto &other_filter =
		    ExpressionFilter::GetExpressionFilter(*other.multi_column_filters[filter_idx], "TableFilterSet::Equals");
		if (!filter.Equals(other_filter)) {
			return false;
		}
	}
	return true;
}

bool TableFilterSet::Equals(TableFilterSet *left, TableFilterSet *right) {
	if (left == right) {
		return true;
	}
	if (!left || !right) {
		return false;
	}
	return left->Equals(*right);
}

unique_ptr<TableFilterSet> TableFilterSet::Copy() const {
	auto copy = make_uniq<TableFilterSet>();
	for (auto &it : filters) {
		copy->filters.emplace(it.first, it.second->Cast<ExpressionFilter>().Copy());
	}
	for (const auto &filter : multi_column_filters) {
		copy->multi_column_filters.push_back(
		    ExpressionFilter::GetExpressionFilter(*filter, "TableFilterSet::Copy").Copy());
	}
	return copy;
}

void TableFilterSet::PushFilter(ProjectionIndex col_idx, unique_ptr<TableFilter> filter) {
	if (!col_idx.IsValid()) {
		throw InternalException("Cannot push a filter over an invalid ProjectionIndex");
	}
	auto &new_filter = ExpressionFilter::GetExpressionFilter(*filter, "TableFilterSet::PushFilter");
	if (!new_filter.column_indexes.empty()) {
		throw InternalException("Cannot push a multi-column filter as a single-column filter");
	}
	auto entry = filters.find(col_idx);
	if (entry == filters.end()) {
		// no filter yet: push the filter directly
		filters[col_idx] = std::move(filter);
	} else {
		// there is already a filter: AND it together
		auto &existing = ExpressionFilter::GetExpressionFilter(*entry->second, "TableFilterSet::PushFilter");
		auto and_expr = make_uniq<BoundConjunctionExpression>(ExpressionType::CONJUNCTION_AND);
		and_expr->GetChildrenMutable().push_back(std::move(existing.expr));
		and_expr->GetChildrenMutable().push_back(std::move(new_filter.expr));
		filters[col_idx] = make_uniq<ExpressionFilter>(std::move(and_expr));
	}
}

void TableFilterSet::PushMultiColumnFilter(unique_ptr<TableFilter> filter) {
	const auto &expression_filter =
	    ExpressionFilter::GetExpressionFilter(*filter, "TableFilterSet::PushMultiColumnFilter");
	if (!expression_filter.expr || expression_filter.column_indexes.size() < 2) {
		throw InternalException("Multi-column filter must reference at least two columns");
	}
	for (const auto &column_index : expression_filter.column_indexes) {
		if (!column_index.IsValid()) {
			throw InternalException("Cannot push a multi-column filter over an invalid ProjectionIndex");
		}
	}
	for (const auto &existing_filter : multi_column_filters) {
		if (ExpressionFilter::GetExpressionFilter(*existing_filter, "TableFilterSet::PushMultiColumnFilter")
		        .Equals(expression_filter)) {
			return;
		}
	}
	multi_column_filters.push_back(std::move(filter));
}

const vector<unique_ptr<TableFilter>> &TableFilterSet::GetMultiColumnFilters() const {
	return multi_column_filters;
}

static bool HasAnyFilters(const TableFilterSet &filters) {
	return filters.HasFilters() || filters.HasMultiColumnFilters();
}

static void AppendTableFilters(TableFilterSet &result, const TableFilterSet &filters) {
	for (auto &entry : filters) {
		result.PushFilter(entry.GetIndex(), entry.Filter().Cast<ExpressionFilter>().Copy());
	}
	for (auto &filter : filters.GetMultiColumnFilters()) {
		result.PushMultiColumnFilter(filter->Cast<ExpressionFilter>().Copy());
	}
}

DynamicTableFilterSet::DynamicTableFilterSet(bool require_registration_p)
    : require_registration(require_registration_p) {
}

void DynamicTableFilterSet::RegisterProducer(const PhysicalOperator &op) {
	lock_guard<mutex> l(lock);
	auto entry = filters.find(op);
	if (entry == filters.end()) {
		filters[op].enabled = !require_registration;
	}
}

void DynamicTableFilterSet::ResetFilters(const PhysicalOperator &op, bool enabled) {
	lock_guard<mutex> l(lock);
	if (require_registration && filters.find(op) == filters.end()) {
		return;
	}
	auto &entry = filters[op];
	entry.enabled = enabled;
	entry.filters.reset();
}

bool DynamicTableFilterSet::HasRegisteredProducers() const {
	lock_guard<mutex> l(lock);
	return !filters.empty();
}

void DynamicTableFilterSet::ResetRegisteredProducers() {
	lock_guard<mutex> l(lock);
	D_ASSERT(require_registration);
	for (auto &entry : filters) {
		entry.second.enabled = false;
		entry.second.filters.reset();
	}
}

void DynamicTableFilterSet::ClearFilters(const PhysicalOperator &op) {
	lock_guard<mutex> l(lock);
	auto entry = filters.find(op);
	if (entry != filters.end()) {
		entry->second.filters.reset();
	}
}

void DynamicTableFilterSet::PushFilter(const PhysicalOperator &op, ProjectionIndex column_index,
                                       unique_ptr<TableFilter> filter) {
	lock_guard<mutex> l(lock);
	if (require_registration && filters.find(op) == filters.end()) {
		return;
	}
	auto &entry = filters[op];
	if (!entry.enabled) {
		return;
	}
	if (!entry.filters) {
		entry.filters = make_uniq<TableFilterSet>();
	}
	entry.filters->PushFilter(column_index, std::move(filter));
}

void DynamicTableFilterSet::AddGroup(shared_ptr<DynamicTableFilterGroup> group) {
	if (require_registration || !group) {
		throw InternalException("Dynamic filter groups require a root filter set");
	}
	lock_guard<mutex> l(lock);
	groups.push_back(std::move(group));
}

vector<shared_ptr<DynamicTableFilterGroup>> DynamicTableFilterSet::GetGroups() const {
	lock_guard<mutex> l(lock);
	return groups;
}

vector<shared_ptr<DynamicTableFilterSet>> DynamicTableFilterSet::GetActiveGroupFilters() const {
	vector<shared_ptr<DynamicTableFilterSet>> result;
	for (auto &group : GetGroups()) {
		if (!group->IsActive() || !group->HasPotentialFilters()) {
			continue;
		}
		for (idx_t i = 0; i < group->ConsumerCount(); i++) {
			result.push_back(group->GetConsumerFilters(i));
		}
	}
	return result;
}

bool DynamicTableFilterSet::HasFilters() const {
	{
		lock_guard<mutex> l(lock);
		for (auto &entry : filters) {
			if (entry.second.enabled && entry.second.filters && HasAnyFilters(*entry.second.filters)) {
				return true;
			}
		}
	}
	for (auto &group : GetGroups()) {
		if (group->HasFilters()) {
			return true;
		}
	}
	return false;
}

unique_ptr<TableFilterSet> DynamicTableFilterSet::CopyFilters(optional_ptr<TableFilterSet> existing_filters) const {
	auto result = make_uniq<TableFilterSet>();
	if (existing_filters) {
		AppendTableFilters(*result, *existing_filters);
	}
	lock_guard<mutex> l(lock);
	for (auto &entry : filters) {
		if (entry.second.enabled && entry.second.filters) {
			AppendTableFilters(*result, *entry.second.filters);
		}
	}
	return result;
}

unique_ptr<TableFilterSet>
DynamicTableFilterSet::GetFinalTableFilters(const PhysicalTableScan &scan,
                                            optional_ptr<TableFilterSet> existing_filters) const {
	auto result = CopyFilters(existing_filters);
	for (auto &group : GetGroups()) {
		group->PushFilterSnapshot(*result);
	}
	if (!HasAnyFilters(*result)) {
		return nullptr;
	}
	return result;
}

DynamicTableFilterGroup::ConsumerFilters::ConsumerFilters(unique_ptr<TableFilterSet> static_filters_p)
    : static_filters(std::move(static_filters_p)), runtime_filters(make_shared_ptr<DynamicTableFilterSet>(true)) {
}

DynamicTableFilterGroup::DynamicTableFilterGroup(vector<unique_ptr<TableFilterSet>> static_filters) {
	for (auto &filters : static_filters) {
		consumers.emplace_back(std::move(filters));
	}
}

idx_t DynamicTableFilterGroup::ConsumerCount() const {
	return consumers.size();
}

shared_ptr<DynamicTableFilterSet> DynamicTableFilterGroup::GetConsumerFilters(idx_t consumer_idx) const {
	return consumers.at(consumer_idx).runtime_filters;
}

bool DynamicTableFilterGroup::HasPotentialFilters() const {
	if (consumers.empty()) {
		return false;
	}
	for (auto &consumer : consumers) {
		if ((!consumer.static_filters || !HasAnyFilters(*consumer.static_filters)) &&
		    !consumer.runtime_filters->HasRegisteredProducers()) {
			return false;
		}
	}
	return true;
}

bool DynamicTableFilterGroup::HasFilters() const {
	if (!IsActive() || consumers.empty()) {
		return false;
	}
	for (auto &consumer : consumers) {
		if ((!consumer.static_filters || !HasAnyFilters(*consumer.static_filters)) &&
		    !consumer.runtime_filters->HasFilters()) {
			return false;
		}
	}
	return true;
}

bool DynamicTableFilterGroup::IsActive() const {
	return active;
}

void DynamicTableFilterGroup::SetActive(bool enabled) {
	active = enabled;
}

void DynamicTableFilterGroup::Reset() {
	active = false;
	for (auto &consumer : consumers) {
		consumer.runtime_filters->ResetRegisteredProducers();
	}
}

static void AppendFilterExpressions(const TableFilterSet &filters, vector<unique_ptr<Expression>> &expressions,
                                    vector<ProjectionIndex> &columns) {
	auto append = [&](const ExpressionFilter &filter, const vector<ProjectionIndex> &inputs) {
		vector<idx_t> indexes;
		for (auto column : inputs) {
			idx_t index = 0;
			for (; index < columns.size(); index++) {
				if (columns[index] == column) {
					break;
				}
			}
			if (index == columns.size()) {
				columns.push_back(column);
			}
			indexes.push_back(index);
		}
		auto expression = filter.expr->Copy();
		ExpressionIterator::VisitExpressionMutable<BoundReferenceExpression>(
		    expression, [&](BoundReferenceExpression &ref, unique_ptr<Expression> &) {
			    ref.IndexMutable() = indexes.at(ref.Index());
		    });
		expressions.push_back(std::move(expression));
	};
	for (auto &entry : filters) {
		append(entry.Filter().Cast<ExpressionFilter>(), {entry.GetIndex()});
	}
	for (auto &filter : filters.GetMultiColumnFilters()) {
		auto &expression_filter = filter->Cast<ExpressionFilter>();
		append(expression_filter, expression_filter.column_indexes);
	}
}

void DynamicTableFilterGroup::PushFilterSnapshot(TableFilterSet &result) const {
	if (!IsActive() || consumers.empty()) {
		return;
	}
	auto disjunction = make_uniq<BoundConjunctionExpression>(ExpressionType::CONJUNCTION_OR);
	vector<ProjectionIndex> columns;
	for (auto &consumer : consumers) {
		auto conjunction = make_uniq<BoundConjunctionExpression>(ExpressionType::CONJUNCTION_AND);
		auto &children = conjunction->GetChildrenMutable();
		if (consumer.static_filters) {
			AppendFilterExpressions(*consumer.static_filters, children, columns);
		}
		auto runtime_filters = consumer.runtime_filters->CopyFilters();
		AppendFilterExpressions(*runtime_filters, children, columns);
		if (children.empty()) {
			return;
		}
		if (children.size() == 1) {
			disjunction->GetChildrenMutable().push_back(std::move(children[0]));
		} else {
			disjunction->GetChildrenMutable().push_back(std::move(conjunction));
		}
	}
	if (columns.empty()) {
		return;
	}
	unique_ptr<Expression> expression;
	if (disjunction->GetChildren().size() == 1) {
		expression = std::move(disjunction->GetChildrenMutable()[0]);
	} else {
		expression = std::move(disjunction);
	}
	if (columns.size() == 1) {
		result.PushFilter(columns[0], make_uniq<ExpressionFilter>(std::move(expression)));
	} else {
		result.PushMultiColumnFilter(make_uniq<ExpressionFilter>(std::move(expression), std::move(columns)));
	}
}

map<ProjectionIndex, unique_ptr<TableFilter>>
TableFilterSet::GetTableFiltersForSerialization(Serializer &serializer) const {
	map<ProjectionIndex, unique_ptr<TableFilter>> result;
	for (auto &entry : filters) {
		auto &expr_filter =
		    ExpressionFilter::GetExpressionFilter(*entry.second, "TableFilterSet::GetTableFiltersForSerialization");
		auto serialized_filter = SerializeExpressionToLegacyFilter(*expr_filter.expr);
		if (!serialized_filter) {
			throw SerializationException(
			    "Could not serialize table filter for projection index %llu to the legacy format",
			    entry.first.GetIndex());
		}
		result.emplace(entry.first, std::move(serialized_filter));
	}
	return result;
}

map<ProjectionIndex, unique_ptr<TableFilter>> &
TableFilterSet::GetTableFiltersForDeserialization(Deserializer &deserializer) {
	return filters;
}

} // namespace duckdb
