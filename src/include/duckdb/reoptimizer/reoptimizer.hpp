//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/reoptimizer/reoptimizer.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/main/client_context.hpp"
#include "duckdb/optimizer/join_order/join_relation.hpp"
#include "duckdb/planner/binder.hpp"
#include "duckdb/planner/column_binding.hpp"
#include "duckdb/planner/logical_operator.hpp"

namespace duckdb {
class Binder;

class ReOptimizer {
public:
	ReOptimizer(ClientContext &context, Binder &binder);

	//! Re-optimization loop until only 2 joins remain
	unique_ptr<LogicalOperator> ReOptimize(unique_ptr<LogicalOperator> plan, const string query);
	//! Simulated re-optimization similar to Stonebraker paper
	unique_ptr<LogicalOperator> SimulatedReOptimize(unique_ptr<LogicalOperator> plan, const string query);

private:
	//! Baseline DuckDB
	unique_ptr<LogicalOperator> AlgorithmBaseline(unique_ptr<LogicalOperator> plan);
	//! Simple re-optimization strategy that performs all filter operations, optimizes, then executes the rest of the plan
	unique_ptr<LogicalOperator> AlgorithmFiltersOnly(unique_ptr<LogicalOperator> plan, const string temporary_table_name);
	//! Simple re-optimization strategy that joins 2 tables at a time
	unique_ptr<LogicalOperator> AlgorithmJoinsOnly(unique_ptr<LogicalOperator> plan, const string temporary_table_name);
	//! Simple re-optimization strategy that performs all filter operations, and 1 join between 2 tables at a time
	unique_ptr<LogicalOperator> AlgorithmOneStep(unique_ptr<LogicalOperator> plan, const string temporary_table_name);
	//! Simple re-optimization strategy that performs joins between X tables at a time
	unique_ptr<LogicalOperator> AlgorithmNStep(idx_t n, unique_ptr<LogicalOperator> plan, const string temporary_table_name);
	//! Set true cardinality of an operator by measuring it
	void SetTrueCardinality(LogicalOperator &plan, LogicalOperator &subquery_plan);
	//! One iteration of the re-optimization procedure
	unique_ptr<LogicalOperator> PerformPartialPlan(unique_ptr<LogicalOperator> plan, LogicalOperator *subquery_plan, const string temporary_table_name);
	//! Gets the true cost of a plan by querying for COUNT(*)
	idx_t GetTrueCost(LogicalOperator &plan);
	//! Queries for the true cardinality of a plan
	idx_t GetTrueCardinality(LogicalOperator &subquery_plan);
	//! Returns all join operators in the plan
	vector<LogicalOperator *> ExtractJoinOperators(LogicalOperator &plan);
	//! Returns all filter operators in the plan
	vector<LogicalOperator *> ExtractFilterOperators(LogicalOperator &plan);
	//! Generate left projection map for joins in the plan, and change right map accordingly
	unique_ptr<LogicalOperator> GenerateProjectionMaps(unique_ptr<LogicalOperator> plan);
	//! Fill binding_name_mapping with (binding -> alias)
	void FindAliases(LogicalOperator &plan);
	//! Creates a CREATE TEMPORARY TABLE query string for the first join to be executed in 'plan'
	string CreateSubQuery(LogicalOperator &plan, vector<string> &queried_tables, vector<string> &where_conditions, bool proj_map_filled);
	//! Reconstructs the filter conditions as strings from a LogicalFilter
	vector<string> GetFilterStrings(LogicalFilter *filter);
	//! Reconstructs an Expression to a SQL string
	string GetExpressionString(Expression *expr);
	//! Reconstructs BoundComparisonExpression to a SQL string
	string GetBoundComparisonString(BoundComparisonExpression *func);
	//! Reconstructs BoundFunctionExpression to a SQL string
	string GetBoundFunctionString(BoundFunctionExpression *func);
	//! Adjusts the original plan by replacing the join with a LogicalGet on the temporary table
	unique_ptr<LogicalOperator> AdjustPlan(unique_ptr<LogicalOperator> plan, LogicalOperator &old_op,
	                                       string temporary_table_name);
	//! Call Catalog::GetTable (works around autocommit stuff)
	TableCatalogEntry *GetTable(string schema, string table_name);
	//! Replaces the join 'old_op' in 'plan' with the given operator 'new_op'
	bool ReplaceLogicalOperator(LogicalOperator &plan, LogicalOperator &old_op, TableCatalogEntry *table,
	                            idx_t depth = 3);
	//! Fixes column bindings after replacing JOIN with GET
	void FixColumnBindings(LogicalOperator &plan);
	//! Fixes column bindings in expressions
	void FixColumnBindings(Expression *expr);
	//! Executes a query in the middle of the re-optimization process
	unique_ptr<QueryResult> ExecuteSubQuery(const string subquery, bool enable_profiling);
	//! Stores the true/estimated cardinality of a plan in saved_cardinalities
	void InjectCardinalities(LogicalOperator &plan, string temp_table_name);
	//! Gets valid subsets of temp_tables
	vector<vector<string>> TempTablePowerset();
	//! Second half of the re-optimization iteration procedure: call Optimizer::Optimize on the adjusted plan
	unique_ptr<LogicalOperator> CallOptimizer(unique_ptr<LogicalOperator> plan);
	//! Empty all left projection maps again (required by PhysicalPlanGenerator - PhysicalHashJoin assert)
	unique_ptr<LogicalOperator> ClearLeftProjectionMaps(unique_ptr<LogicalOperator> plan);

	//! Utility to join strings like in Java, Python
	string JoinStrings(vector<string> strings, string delimiter);

	//! The client context
	ClientContext &context;
	//! The binder
	Binder &binder;

	//! Binding to column alias mapping (created by CreateMaps)
	unordered_map<string, string> binding_name_mapping;
	//! The new column bindings (after replacing an operator with GET)
	unordered_map<string, ColumnBinding> rebind_mapping;

	//! Stores true/estimated cardinalities of sets of relations
	unordered_map<string, idx_t> cardinalities;
	//! Stores pairs of (set of relations, temp table name) of subqueries that were executed
	unordered_map<string, vector<string>> temp_table_relations;
	//! vector of temp tables made so far
	vector<string> temp_tables;

	//! Whether we are done re-optimizing the plan
	bool done = false;

	//! Whether to compute the cost of the plan
	bool compute_cost = false;
	//! The cost of the plan
	idx_t plan_cost = 0;

};

} // namespace duckdb
