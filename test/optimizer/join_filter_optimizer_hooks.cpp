#include "catch.hpp"
#include "test_helpers.hpp"

#include "duckdb/main/config.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/optimizer/optimizer_extension.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_cross_product.hpp"
#include "duckdb/planner/operator/logical_cteref.hpp"
#include "duckdb/planner/operator/logical_materialized_cte.hpp"

using namespace duckdb;

struct JoinFilterHookTest : OptimizerExtensionInfo {
	idx_t ordinary_filters = 0;
	bool add_consumer = false;
	bool added_consumer = false;

	static void Inspect(LogicalOperator &op, JoinFilterHookTest &info) {
		if (op.type == LogicalOperatorType::LOGICAL_COMPARISON_JOIN) {
			auto &join = op.Cast<LogicalComparisonJoin>();
			if (join.filter_pushdown) {
				info.ordinary_filters += join.filter_pushdown->probe_info.size();
			}
		}
		for (auto &child : op.children) {
			Inspect(*child, info);
		}
	}

	static bool AddConsumer(LogicalOperator &op, LogicalMaterializedCTE &cte, Binder &binder) {
		if (op.type == LogicalOperatorType::LOGICAL_AGGREGATE_AND_GROUP_BY) {
			auto ref = make_uniq<LogicalCTERef>(binder.GenerateTableIndex(), cte.table_index, cte.children[0]->types,
			                                    vector<Identifier> {"k"});
			op.children[0] = LogicalCrossProduct::Create(std::move(op.children[0]), std::move(ref));
			return true;
		}
		for (auto &child : op.children) {
			if (AddConsumer(*child, cte, binder)) {
				return true;
			}
		}
		return false;
	}

	static void Optimize(OptimizerExtensionInput &input, unique_ptr<LogicalOperator> &plan) {
		auto &info = static_cast<JoinFilterHookTest &>(*input.info);
		Inspect(*plan, info);
		if (info.add_consumer && plan->type == LogicalOperatorType::LOGICAL_MATERIALIZED_CTE) {
			auto &cte = plan->Cast<LogicalMaterializedCTE>();
			cte.children[0]->ResolveOperatorTypes();
			info.added_consumer = AddConsumer(*cte.children[1], cte, input.optimizer.binder);
		}
	}
};

TEST_CASE("Join filter generation preserves optimizer hooks and uses final CTE consumers",
          "[optimizer][dynamic_filter_groups]") {
	DuckDB db(nullptr);
	Connection con(db);
	REQUIRE_NO_FAIL(con.Query("SET disabled_optimizers='cte_filter_pusher,statistics_propagation,join_order,"
	                          "build_side_probe_side,compressed_materialization'"));
	REQUIRE_NO_FAIL(con.Query("CREATE TABLE probe AS SELECT i AS k FROM range(1000) t(i)"));
	REQUIRE_NO_FAIL(con.Query("CREATE TABLE build AS SELECT i AS k FROM range(3) t(i)"));
	auto info = make_shared_ptr<JoinFilterHookTest>();
	OptimizerExtension extension;
	extension.optimizer_info = info;
	extension.optimize_function = JoinFilterHookTest::Optimize;
	OptimizerExtension::Register(DBConfig::GetConfig(*con.context), std::move(extension));

	REQUIRE_NO_FAIL(con.Query("SELECT count(*) FROM probe JOIN build USING (k)"));
	REQUIRE(info->ordinary_filters > 0);
	info->add_consumer = true;
	auto result = con.Query("WITH p AS MATERIALIZED (SELECT k FROM probe) "
	                        "SELECT count(*) FROM p a JOIN build b ON a.k=b.k "
	                        "CROSS JOIN (p c JOIN build d ON c.k=d.k)");
	REQUIRE_NO_FAIL(*result);
	REQUIRE(info->added_consumer);
	REQUIRE(CHECK_COLUMN(result, 0, {9000}));
}
