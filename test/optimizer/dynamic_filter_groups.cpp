#include "catch.hpp"
#include "test_helpers.hpp"

#include "duckdb/common/vector/vector_writer.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/execution/operator/scan/physical_dummy_scan.hpp"
#include "duckdb/execution/operator/scan/physical_table_scan.hpp"
#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/expression/bound_reference_expression.hpp"
#include "duckdb/planner/expression_iterator.hpp"
#include "duckdb/planner/filter/expression_filter.hpp"
#include "duckdb/planner/filter/table_filter_functions.hpp"

using namespace duckdb;

static unique_ptr<TableFilterSet> GroupTestStaticFilter(idx_t column, int32_t value) {
	auto result = make_uniq<TableFilterSet>();
	result->PushFilter(ProjectionIndex(column),
	                   ExpressionFilter::CreateComparisonFilter(ExpressionType::COMPARE_EQUAL, Value::INTEGER(value)));
	return result;
}

static bool EvaluateGroupTestFilter(ClientContext &context, const TableFilterSet &filters,
                                    const vector<Value> &values) {
	auto evaluate = [&](const ExpressionFilter &filter, const vector<ProjectionIndex> &columns) {
		auto expression = filter.expr->Copy();
		ExpressionIterator::VisitExpressionMutable<BoundReferenceExpression>(
		    expression, [&](BoundReferenceExpression &ref, unique_ptr<Expression> &owned) {
			    owned = make_uniq<BoundConstantExpression>(values.at(columns.at(ref.Index()).GetIndex()));
		    });
		auto result = ExpressionExecutor::EvaluateScalar(context, *expression);
		return !result.IsNull() && result.GetValue<bool>();
	};
	for (auto &entry : filters) {
		if (!evaluate(entry.Filter().Cast<ExpressionFilter>(), {entry.GetIndex()})) {
			return false;
		}
	}
	for (auto &filter : filters.GetMultiColumnFilters()) {
		auto &expression_filter = filter->Cast<ExpressionFilter>();
		if (!evaluate(expression_filter, expression_filter.column_indexes)) {
			return false;
		}
	}
	return true;
}

struct FilterGroupTest {
	FilterGroupTest()
	    : db(nullptr), con(db), plan(Allocator::Get(*con.context)),
	      first(plan.Make<PhysicalDummyScan>(vector<LogicalType> {}, 1)),
	      second(plan.Make<PhysicalDummyScan>(vector<LogicalType> {}, 1)),
	      scan(plan.Make<PhysicalTableScan>(
	                   vector<LogicalType> {LogicalType::INTEGER, LogicalType::INTEGER},
	                   BoundTableFunction(TableFunction("filter_group", vector<LogicalType> {}, nullptr)), nullptr,
	                   vector<LogicalType> {LogicalType::INTEGER, LogicalType::INTEGER},
	                   vector<ColumnIndex> {ColumnIndex(0), ColumnIndex(1)}, vector<idx_t> {0, 1},
	                   vector<string> {"a", "b"}, nullptr, 100, ExtraOperatorInfo(), vector<Value> {},
	                   virtual_column_map_t {})
	               .Cast<PhysicalTableScan>()) {
	}

	void Push(DynamicTableFilterSet &filters, const PhysicalOperator &producer, idx_t column, int32_t value) {
		filters.PushFilter(
		    producer, ProjectionIndex(column),
		    ExpressionFilter::CreateComparisonFilter(ExpressionType::COMPARE_EQUAL, Value::INTEGER(value)));
	}

	DuckDB db;
	Connection con;
	PhysicalPlan plan;
	PhysicalOperator &first;
	PhysicalOperator &second;
	PhysicalTableScan &scan;
};

TEST_CASE("Consumer filter groups preserve every alternative through partial publication", "[dynamic_filter_groups]") {
	FilterGroupTest test;
	DynamicTableFilterSet root;
	vector<unique_ptr<TableFilterSet>> static_filters;
	static_filters.push_back(GroupTestStaticFilter(0, 1));
	static_filters.push_back(GroupTestStaticFilter(0, 3));
	auto group = make_shared_ptr<DynamicTableFilterGroup>(std::move(static_filters));
	root.AddGroup(group);
	auto first = group->GetConsumerFilters(0);
	auto second = group->GetConsumerFilters(1);
	first->RegisterProducer(test.first);
	second->RegisterProducer(test.second);
	REQUIRE(group->HasPotentialFilters());
	REQUIRE_FALSE(root.HasFilters());
	REQUIRE(root.GetActiveGroupFilters().empty());
	group->SetActive(true);
	REQUIRE(root.GetActiveGroupFilters().size() == 2);

	first->ResetFilters(test.first, true);
	test.Push(*first, test.first, 1, 2);
	auto partial = root.GetFinalTableFilters(test.scan, nullptr);
	REQUIRE(partial);
	REQUIRE(partial->HasMultiColumnFilters());
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *partial, {Value::INTEGER(1), Value::INTEGER(2)}));
	REQUIRE_FALSE(EvaluateGroupTestFilter(*test.con.context, *partial, {Value::INTEGER(1), Value::INTEGER(4)}));
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *partial, {Value::INTEGER(3), Value::INTEGER(2)}));
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *partial, {Value::INTEGER(3), Value::INTEGER(4)}));

	second->ResetFilters(test.second, true);
	test.Push(*second, test.second, 1, 4);
	auto complete = root.GetFinalTableFilters(test.scan, nullptr);
	REQUIRE(complete);
	for (int32_t a = 0; a <= 4; a++) {
		for (int32_t b = 0; b <= 5; b++) {
			const auto expected = (a == 1 && b == 2) || (a == 3 && b == 4);
			REQUIRE(EvaluateGroupTestFilter(*test.con.context, *complete, {Value::INTEGER(a), Value::INTEGER(b)}) ==
			        expected);
		}
	}
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *partial, {Value::INTEGER(3), Value::INTEGER(2)}));
	group->Reset();
	REQUIRE_FALSE(root.HasFilters());
	test.Push(*first, test.first, 1, 5);
	test.Push(*second, test.second, 1, 5);
	REQUIRE_FALSE(first->HasFilters());
	REQUIRE_FALSE(second->HasFilters());
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *complete, {Value::INTEGER(1), Value::INTEGER(2)}));

	group->SetActive(true);
	auto next = root.GetFinalTableFilters(test.scan, nullptr);
	REQUIRE(next);
	REQUIRE_FALSE(next->HasMultiColumnFilters());
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *next, {Value::INTEGER(1)}));
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *next, {Value::INTEGER(3)}));
	REQUIRE_FALSE(EvaluateGroupTestFilter(*test.con.context, *next, {Value::INTEGER(2)}));
}

TEST_CASE("Unrestricted consumer alternatives weaken groups without weakening ordinary filters",
          "[dynamic_filter_groups]") {
	FilterGroupTest test;
	DynamicTableFilterSet root;
	vector<unique_ptr<TableFilterSet>> static_filters(2);
	auto group = make_shared_ptr<DynamicTableFilterGroup>(std::move(static_filters));
	root.AddGroup(group);
	auto first = group->GetConsumerFilters(0);
	auto second = group->GetConsumerFilters(1);
	first->RegisterProducer(test.first);
	group->SetActive(true);
	first->ResetFilters(test.first, true);
	test.Push(*first, test.first, 0, 1);
	REQUIRE_FALSE(group->HasPotentialFilters());
	REQUIRE(root.GetActiveGroupFilters().empty());
	REQUIRE_FALSE(root.HasFilters());
	REQUIRE_FALSE(root.GetFinalTableFilters(test.scan, nullptr));

	second->ResetFilters(test.second, true);
	test.Push(*second, test.second, 0, 2);
	REQUIRE_FALSE(second->HasFilters());
	second->RegisterProducer(test.second);
	test.Push(*second, test.second, 0, 2);
	REQUIRE_FALSE(second->HasFilters());
	REQUIRE(group->HasPotentialFilters());
	REQUIRE(root.GetActiveGroupFilters().size() == 2);
	REQUIRE_FALSE(root.HasFilters());

	second->ResetFilters(test.second, true);
	test.Push(*second, test.second, 0, 2);
	REQUIRE(root.HasFilters());
	auto snapshot = root.GetFinalTableFilters(test.scan, nullptr);
	REQUIRE(snapshot);
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(1)}));
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(2)}));
	REQUIRE_FALSE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(3)}));

	test.Push(root, test.first, 0, 2);
	auto ordinary = root.GetFinalTableFilters(test.scan, nullptr);
	REQUIRE(ordinary);
	REQUIRE_FALSE(EvaluateGroupTestFilter(*test.con.context, *ordinary, {Value::INTEGER(1)}));
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *ordinary, {Value::INTEGER(2)}));
	second->ClearFilters(test.second);
	auto weakened = root.GetFinalTableFilters(test.scan, nullptr);
	REQUIRE(weakened);
	REQUIRE_FALSE(EvaluateGroupTestFilter(*test.con.context, *weakened, {Value::INTEGER(1)}));
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *weakened, {Value::INTEGER(2)}));
#ifndef DUCKDB_CRASH_ON_ASSERT
	REQUIRE_THROWS(first->AddGroup(group));
#endif
}

TEST_CASE("Static multi-column alternatives and multiple groups preserve conjunctions", "[dynamic_filter_groups]") {
	FilterGroupTest test;
	DynamicTableFilterSet root;
	vector<unique_ptr<TableFilterSet>> static_filters;
	auto paired = make_uniq<TableFilterSet>();
	paired->PushMultiColumnFilter(make_uniq<ExpressionFilter>(
	    BoundComparisonExpression::Create(ExpressionType::COMPARE_EQUAL,
	                                      make_uniq<BoundReferenceExpression>(LogicalType::INTEGER, idx_t(0)),
	                                      make_uniq<BoundReferenceExpression>(LogicalType::INTEGER, idx_t(1))),
	    vector<ProjectionIndex> {ProjectionIndex(1), ProjectionIndex(0)}));
	static_filters.push_back(std::move(paired));
	static_filters.push_back(GroupTestStaticFilter(0, 3));
	auto first = make_shared_ptr<DynamicTableFilterGroup>(std::move(static_filters));
	first->SetActive(true);
	root.AddGroup(first);
	vector<unique_ptr<TableFilterSet>> other;
	other.push_back(GroupTestStaticFilter(1, 2));
	other.push_back(GroupTestStaticFilter(1, 4));
	auto second = make_shared_ptr<DynamicTableFilterGroup>(std::move(other));
	second->SetActive(true);
	root.AddGroup(second);
	auto snapshot = root.GetFinalTableFilters(test.scan, nullptr);
	REQUIRE(snapshot);
	for (int32_t a = 0; a < 6; a++) {
		for (int32_t b = 0; b < 6; b++) {
			const auto expected = (a == b || a == 3) && (b == 2 || b == 4);
			REQUIRE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(a), Value::INTEGER(b)}) ==
			        expected);
		}
	}
}

TEST_CASE("Each available producer remains an independent consumer conjunct", "[dynamic_filter_groups]") {
	FilterGroupTest test;
	for (idx_t mask = 0; mask < 16; mask++) {
		DynamicTableFilterSet root;
		vector<unique_ptr<TableFilterSet>> static_filters(2);
		auto group = make_shared_ptr<DynamicTableFilterGroup>(std::move(static_filters));
		root.AddGroup(group);
		for (idx_t consumer_idx = 0; consumer_idx < 2; consumer_idx++) {
			auto leaf = group->GetConsumerFilters(consumer_idx);
			leaf->RegisterProducer(test.first);
			leaf->RegisterProducer(test.second);
			leaf->ResetFilters(test.first, (mask & (1ULL << (2 * consumer_idx))) != 0);
			leaf->ResetFilters(test.second, (mask & (1ULL << (2 * consumer_idx + 1))) != 0);
			test.Push(*leaf, test.first, 0, int32_t(consumer_idx + 1));
			test.Push(*leaf, test.second, 1, int32_t(consumer_idx + 3));
		}
		group->SetActive(true);
		auto snapshot = root.GetFinalTableFilters(test.scan, nullptr);
		const auto unrestricted = (mask & 3) == 0 || (mask & 12) == 0;
		REQUIRE(bool(snapshot) == !unrestricted);
		if (!snapshot) {
			continue;
		}
		for (int32_t a = 0; a <= 3; a++) {
			for (int32_t b = 0; b <= 4; b++) {
				const auto first = (!(mask & 1) || a == 1) && (!(mask & 2) || b == 3);
				const auto second = (!(mask & 4) || a == 2) && (!(mask & 8) || b == 4);
				REQUIRE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(a), Value::INTEGER(b)}) ==
				        (first || second));
			}
		}
	}
}

TEST_CASE("Constant and NULL consumer predicates remain valid filter snapshots", "[dynamic_filter_groups]") {
	FilterGroupTest test;
	idx_t consumers = 2;
	SECTION("Multiple consumers") {
	}
	SECTION("One consumer") {
		consumers = 1;
	}
	for (auto value : {Value::BOOLEAN(false), Value(LogicalType::BOOLEAN), Value::BOOLEAN(true)}) {
		DynamicTableFilterSet root;
		vector<unique_ptr<TableFilterSet>> static_filters;
		for (idx_t i = 0; i < consumers; i++) {
			auto filters = make_uniq<TableFilterSet>();
			filters->PushFilter(ProjectionIndex(0),
			                    make_uniq<ExpressionFilter>(make_uniq<BoundConstantExpression>(value)));
			static_filters.push_back(std::move(filters));
		}
		auto group = make_shared_ptr<DynamicTableFilterGroup>(std::move(static_filters));
		root.AddGroup(group);
		group->SetActive(true);
		REQUIRE(root.HasFilters());
		auto snapshot = root.GetFinalTableFilters(test.scan, nullptr);
		REQUIRE(snapshot);
		REQUIRE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(0)}) ==
		        (!value.IsNull() && value.GetValue<bool>()));
	}
}

TEST_CASE("Composed prefix-range snapshots retain their payload after owner reset", "[dynamic_filter_groups]") {
	FilterGroupTest test;
	unique_ptr<TableFilterSet> snapshot;
	{
		DynamicTableFilterSet root;
		vector<unique_ptr<TableFilterSet>> static_filters;
		static_filters.push_back(nullptr);
		static_filters.push_back(GroupTestStaticFilter(1, 7));
		auto group = make_shared_ptr<DynamicTableFilterGroup>(std::move(static_filters));
		root.AddGroup(group);
		auto leaf = group->GetConsumerFilters(0);
		leaf->RegisterProducer(test.first);
		leaf->ResetFilters(test.first, true);
		shared_ptr<PrefixRangeFilter> payload = PrefixRangeFilter::CreatePrefixRangeFilter(LogicalType::INTEGER);
		payload->Initialize(*test.con.context, 2, Value::INTEGER(1), Value::INTEGER(3), 128);
		Vector keys(LogicalType::INTEGER, 2);
		{
			auto writer = FlatVector::Writer<int32_t>(keys, 2);
			writer.WriteValue(1);
			writer.WriteValue(3);
		}
		auto state = payload->InitializeBuildState(*test.con.context);
		payload->InsertKeys(keys, *state);
		payload->MergeBuildState(*state);
		vector<unique_ptr<Expression>> children;
		children.push_back(make_uniq<BoundReferenceExpression>(LogicalType::INTEGER, idx_t(0)));
		leaf->PushFilter(
		    test.first, ProjectionIndex(0),
		    make_uniq<ExpressionFilter>(make_uniq<BoundFunctionExpression>(
		        BoundScalarFunction(PrefixRangeScalarFun::GetFunction(LogicalType::INTEGER)), std::move(children),
		        make_uniq<PrefixRangeFunctionData>(payload, true, "a", LogicalType::INTEGER, 0.0f, idx_t(0)))));
		group->SetActive(true);
		snapshot = root.GetFinalTableFilters(test.scan, nullptr);
		group->Reset();
		payload.reset();
	}
	REQUIRE(snapshot);
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(1), Value::INTEGER(0)}));
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(3), Value::INTEGER(0)}));
	REQUIRE_FALSE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(4), Value::INTEGER(0)}));
	REQUIRE(EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value::INTEGER(4), Value::INTEGER(7)}));
	REQUIRE_FALSE(
	    EvaluateGroupTestFilter(*test.con.context, *snapshot, {Value(LogicalType::INTEGER), Value::INTEGER(0)}));
}
