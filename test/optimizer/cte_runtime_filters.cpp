#include "catch.hpp"
#include "test_helpers.hpp"

#include "duckdb/execution/executor.hpp"
#include "duckdb/execution/operator/join/join_filter_pushdown.hpp"
#include "duckdb/execution/operator/scan/physical_dummy_scan.hpp"
#include "duckdb/execution/operator/set/physical_cte.hpp"
#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/main/config.hpp"
#include "duckdb/optimizer/join_filter_pushdown_optimizer.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/optimizer/optimizer_extension.hpp"
#include "duckdb/parallel/meta_pipeline.hpp"
#include "duckdb/parallel/pipeline.hpp"
#include "duckdb/parallel/pipeline_broadcast_exchange.hpp"
#include "duckdb/parser/parser.hpp"
#include "duckdb/planner/filter/expression_filter.hpp"
#include "duckdb/planner/operator/logical_cteref.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_materialized_cte.hpp"
#include "duckdb/planner/planner.hpp"

using namespace duckdb;

static unique_ptr<LogicalOperator> OptimizeCTEFilterTest(Connection &con, bool table_parameter = false) {
	Parser parser(*con.context);
	string sql = "WITH p AS MATERIALIZED (SELECT k FROM a) "
	             "SELECT p.k FROM p JOIN b ON p.k=b.k "
	             "UNION ALL SELECT p.k FROM p JOIN c ON p.k=c.k";
	parser.ParseQuery(table_parameter ? "SELECT count(*) FROM summary((" + sql + "))" : sql);
	Planner planner(*con.context);
	planner.CreatePlan(std::move(parser.statements[0]));
	Optimizer optimizer(*planner.binder, *con.context);
	auto plan = optimizer.Optimize(std::move(planner.plan));
	plan->ResolveOperatorTypes();
	return plan;
}

template <class T>
static optional_ptr<T> FindCTEFilterTestOperator(LogicalOperator &op) {
	if (op.type == T::TYPE) {
		return op.Cast<T>();
	}
	for (auto &child : op.children) {
		auto result = FindCTEFilterTestOperator<T>(*child);
		if (result) {
			return result;
		}
	}
	return nullptr;
}

static optional_ptr<PhysicalCTE> FindCTEFilterTestOwner(PhysicalOperator &op) {
	if (op.type == PhysicalOperatorType::CTE) {
		return op.Cast<PhysicalCTE>();
	}
	for (auto &child : op.children) {
		auto result = FindCTEFilterTestOwner(child);
		if (result) {
			return result;
		}
	}
	return nullptr;
}

static void InitializeCTEFilterTest(Connection &con) {
	REQUIRE_NO_FAIL(con.Query("CREATE TABLE a AS SELECT i::BIGINT k FROM range(100) t(i)"));
	REQUIRE_NO_FAIL(con.Query("CREATE TABLE b AS SELECT 1::BIGINT k"));
	REQUIRE_NO_FAIL(con.Query("CREATE TABLE c AS SELECT 3::BIGINT k"));
	REQUIRE_NO_FAIL(con.Query("SET disabled_optimizers="
	                          "'cte_filter_pusher,join_order,build_side_probe_side,statistics_propagation'"));
	con.BeginTransaction();
}

TEST_CASE("CTE filter endpoints require explicit terminal discovery", "[cte_runtime_filters]") {
	LogicalCTERef ref(TableIndex(10), TableIndex(1), {LogicalType::BIGINT, LogicalType::INTEGER}, {"a", "b"});
	vector<JoinFilterPushdownColumn> columns(3);
	columns[0].probe_column_index = ColumnBinding(TableIndex(10), ProjectionIndex(0));
	columns[1].probe_column_index = ColumnBinding(TableIndex(10), ProjectionIndex(1));
	columns[2].probe_column_index = ColumnBinding(TableIndex(20), ProjectionIndex(0));
	vector<PushdownFilterTarget> targets;
	JoinFilterPushdownOptimizer::GetPushdownFilterTargets(ref, columns, targets);
	REQUIRE(targets.empty());
	JoinFilterPushdownOptimizer::GetPushdownFilterTargets(ref, columns, targets, EqualityPropagation::DISABLED,
	                                                      CTEFilterTargetMode::INCLUDE);
	REQUIRE(targets.size() == 1);
	REQUIRE(targets[0].IsCTEReference());
	REQUIRE(&targets[0].GetCTEReference() == &ref);
	REQUIRE(targets[0].columns.size() == 2);
	REQUIRE(targets[0].columns[0].storage_type == LogicalType::BIGINT);
	REQUIRE(targets[0].columns[1].storage_type == LogicalType::INTEGER);
	ref.correlated_columns = 1;
	targets.clear();
	JoinFilterPushdownOptimizer::GetPushdownFilterTargets(ref, columns, targets, EqualityPropagation::DISABLED,
	                                                      CTEFilterTargetMode::INCLUDE);
	REQUIRE(targets.empty());
}

TEST_CASE("Final CTE filter mapping rejects stale endpoints and unsupported producer routes", "[cte_runtime_filters]") {
	DuckDB db(nullptr);
	Connection con(db);
	InitializeCTEFilterTest(con);
	auto plan = OptimizeCTEFilterTest(con);
	auto cte = FindCTEFilterTestOperator<LogicalMaterializedCTE>(*plan);
	REQUIRE(cte);
	REQUIRE(cte->filter_pushdown);
	REQUIRE(cte->filter_pushdown->consumers.size() == 2);
	bool expected = false;
	SECTION("Unchanged endpoints and routes") {
		expected = true;
	}
	SECTION("Changed consumer identity") {
		auto ref = FindCTEFilterTestOperator<LogicalCTERef>(*plan);
		REQUIRE(ref);
		ref->table_index = TableIndex(999999);
	}
	SECTION("Removed consumer from this CTE") {
		auto ref = FindCTEFilterTestOperator<LogicalCTERef>(*plan);
		REQUIRE(ref);
		ref->cte_index = TableIndex(999999);
	}
	SECTION("Changed consumer type") {
		auto ref = FindCTEFilterTestOperator<LogicalCTERef>(*plan);
		REQUIRE(ref);
		ref->chunk_types[0] = LogicalType::VARCHAR;
	}
	SECTION("Changed producer type") {
		cte->children[0]->types[0] = LogicalType::VARCHAR;
	}
	SECTION("Unsupported final scan") {
		auto get = FindCTEFilterTestOperator<LogicalGet>(*cte->children[0]);
		REQUIRE(get);
		get->function.filter_pushdown = false;
	}
	JoinFilterPushdownOptimizer::FinalizeCTEFilters(*plan);
	const auto has_targets = cte->filter_pushdown && !cte->filter_pushdown->targets.empty();
	REQUIRE(has_targets == expected);
	if (expected) {
		REQUIRE(cte->filter_pushdown->targets.size() == 1);
		REQUIRE(cte->filter_pushdown->targets[0].group->ConsumerCount() == 2);
		REQUIRE_FALSE(cte->filter_pushdown->targets[0].group->IsActive());
	}
}

TEST_CASE("Physical CTE filter groups follow resolved consumer modes", "[cte_runtime_filters]") {
	DuckDB db(nullptr);
	Connection con(db);
	InitializeCTEFilterTest(con);
	PhysicalPlanGenerator generator(*con.context);
	auto plan = generator.Plan(OptimizeCTEFilterTest(con));
	auto cte = FindCTEFilterTestOwner(plan->Root());
	REQUIRE(cte);
	REQUIRE(cte->filter_groups.size() == 1);
	REQUIRE(cte->cte_scans.size() == 2);
	REQUIRE(cte->exchange);
	REQUIRE_FALSE(cte->filter_groups[0]->IsActive());
	bool materialized = false;
	SECTION("Both consumers materialized") {
		materialized = true;
		cte->RegisterMaterializedConsumer(0);
		cte->RegisterMaterializedConsumer(1);
	}
	SECTION("Mixed consumer modes") {
		cte->RegisterMaterializedConsumer(0);
		cte->exchange->SelectBufferedConsumer(1, PipelineBroadcastExchangeScanMode::CHUNK);
	}
	SECTION("Streaming consumers") {
		cte->exchange->SelectBufferedConsumer(0, PipelineBroadcastExchangeScanMode::CHUNK);
		cte->exchange->SelectBufferedConsumer(1, PipelineBroadcastExchangeScanMode::CHUNK);
	}
	REQUIRE_FALSE(cte->filter_groups[0]->IsActive());
	cte->SetPipelineSelectionResolved();
	REQUIRE(cte->filter_groups[0]->IsActive() == materialized);
}

TEST_CASE("Post-optimizer extensions receive no captured CTE filter proof", "[cte_runtime_filters]") {
	DBConfig config;
	OptimizerExtension extension;
	extension.optimize_function = [](OptimizerExtensionInput &, unique_ptr<LogicalOperator> &plan) {
		auto cte = FindCTEFilterTestOperator<LogicalMaterializedCTE>(*plan);
		if (cte) {
			REQUIRE_FALSE(cte->filter_pushdown);
		}
	};
	OptimizerExtension::Register(config, std::move(extension));
	DuckDB db(nullptr, &config);
	Connection con(db);
	InitializeCTEFilterTest(con);
	auto plan = OptimizeCTEFilterTest(con);
	auto cte = FindCTEFilterTestOperator<LogicalMaterializedCTE>(*plan);
	REQUIRE(cte);
	REQUIRE_FALSE(cte->filter_pushdown);
	JoinFilterPushdownOptimizer::FinalizeCTEFilters(*plan);
	REQUIRE_FALSE(cte->filter_pushdown);
}

TEST_CASE("Table-parameter planning finalizes each CTE filter group once", "[cte_runtime_filters]") {
	DuckDB db(nullptr);
	Connection con(db);
	InitializeCTEFilterTest(con);
	PhysicalPlanGenerator generator(*con.context);
	auto plan = generator.Plan(OptimizeCTEFilterTest(con, true));
	auto cte = FindCTEFilterTestOwner(plan->Root());
	REQUIRE(cte);
	REQUIRE(cte->filter_groups.size() == 1);
	auto &group = *cte->filter_groups[0];
	REQUIRE(group.ConsumerCount() == 2);
	REQUIRE(group.HasPotentialFilters());
}

TEST_CASE("CTE pipeline construction resets producers absent from the schedule", "[cte_runtime_filters]") {
	DuckDB db(nullptr);
	Connection con(db);
	InitializeCTEFilterTest(con);
	PhysicalPlanGenerator generator(*con.context);
	auto plan = generator.Plan(OptimizeCTEFilterTest(con));
	auto cte = FindCTEFilterTestOwner(plan->Root());
	REQUIRE(cte);
	REQUIRE(cte->filter_groups.size() == 1);
	auto &absent = plan->Make<PhysicalDummyScan>(vector<LogicalType> {}, 1);
	auto group = cte->filter_groups[0];
	auto leaf = group->GetConsumerFilters(0);
	leaf->RegisterProducer(absent);
	for (idx_t iteration = 0; iteration < 2; iteration++) {
		group->SetActive(true);
		leaf->ResetFilters(absent, true);
		leaf->PushFilter(absent, ProjectionIndex(0),
		                 ExpressionFilter::CreateComparisonFilter(ExpressionType::COMPARE_EQUAL, Value::BIGINT(42)));
		REQUIRE(leaf->HasFilters());
		Executor executor(*con.context);
		PipelineBuildState state;
		auto meta = make_shared_ptr<MetaPipeline>(executor, state, nullptr);
		cte->BuildPipelines(*meta->GetBasePipeline(), *meta);
		REQUIRE_FALSE(leaf->HasFilters());
		leaf->PushFilter(absent, ProjectionIndex(0),
		                 ExpressionFilter::CreateComparisonFilter(ExpressionType::COMPARE_EQUAL, Value::BIGINT(43)));
		REQUIRE_FALSE(leaf->HasFilters());
	}
}
