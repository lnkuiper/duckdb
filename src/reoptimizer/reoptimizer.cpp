#include "duckdb/reoptimizer/reoptimizer.hpp"

#include "duckdb/catalog/catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/unordered_map.hpp"
#include "duckdb/common/enums/catalog_type.hpp"
#include "duckdb/common/enums/expression_type.hpp"
#include "duckdb/common/enums/logical_operator_type.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/optimizer/join_order_optimizer.hpp"
#include "duckdb/planner/expression/bound_operator_expression.hpp"
#include "duckdb/planner/expression/bound_aggregate_expression.hpp"
#include "duckdb/planner/expression/bound_between_expression.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/expression/bound_conjunction_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/planner/operator/logical_chunk_get.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"
#include "duckdb/planner/operator/logical_join.hpp"
#include "duckdb/planner/operator/logical_projection.hpp"
#include "duckdb/storage/data_table.hpp"

using namespace duckdb;
using namespace std;

ReOptimizer::ReOptimizer(ClientContext &context, Binder &binder) : context(context), binder(binder) {
}

unique_ptr<LogicalOperator> ReOptimizer::ReOptimize(unique_ptr<LogicalOperator> plan, const string query) {
	compute_cost = true;
	const string tablename_prefix = "_reopt_temp_" + to_string(hash<string>{}(query));
	// re-optimization loop
	for (int iter = 0; true; iter++) {
		const string temp_table_name = tablename_prefix + "_" + to_string(iter);
		plan = AlgorithmSmartStep(2, move(plan), temp_table_name);
		if (done) {
			break;
		}
	}
	if (compute_cost) {
		binding_name_mapping.clear();
		FindAliases(*plan);
		plan_cost += GetTrueCost(*plan);
		Printer::Print(to_string(plan_cost));
	}
	//	Printer::Print(to_string(materialize_size));
	//	Printer::Print(to_string(plan->EstimateCost()));
	return plan;
}

unique_ptr<LogicalOperator> ReOptimizer::AlgorithmBaseline(unique_ptr<LogicalOperator> plan) {
	done = true;
	return plan;
}

unique_ptr<LogicalOperator> ReOptimizer::AlgorithmFiltersOnly(unique_ptr<LogicalOperator> plan,
															  const string temporary_table_name) {
	vector<LogicalOperator *> filters = ExtractFilterOperators(*plan);
	if (filters.empty()) {
		context.profiler.StartPhase("optimizer");
		plan = CallOptimizer(move(plan));
		context.profiler.EndPhase();

		done = true;
		return plan;
	}
	return PerformPartialPlan(move(plan), filters.back(), temporary_table_name);
}

unique_ptr<LogicalOperator> ReOptimizer::AlgorithmJoinsOnly(unique_ptr<LogicalOperator> plan,
															const string temporary_table_name) {
	vector<LogicalOperator *> joins = ExtractJoinOperators(*plan);
	if (joins.size() <= 2) {
		done = true;
		return plan;
	}

	plan = PerformPartialPlan(move(plan), joins.back(), temporary_table_name);

	context.profiler.StartPhase("optimizer");
	plan = CallOptimizer(move(plan));
	context.profiler.EndPhase();

	return plan;
}

unique_ptr<LogicalOperator> ReOptimizer::AlgorithmOneStep(unique_ptr<LogicalOperator> plan,
														  const string temporary_table_name) {
	vector<LogicalOperator *> filters = ExtractFilterOperators(*plan);
	if (filters.size() == 1) {
		plan = PerformPartialPlan(move(plan), filters.back(), temporary_table_name);

		context.profiler.StartPhase("optimizer");
		plan = CallOptimizer(move(plan));
		context.profiler.EndPhase();

		return plan;
	} else if (!filters.empty()) {
		return PerformPartialPlan(move(plan), filters.back(), temporary_table_name);
	}

	vector<LogicalOperator *> joins = ExtractJoinOperators(*plan);
	if (joins.size() <= 2) {
		done = true;
		return plan;
	}

	plan = PerformPartialPlan(move(plan), joins.back(), temporary_table_name);

	context.profiler.StartPhase("optimizer");
	plan = CallOptimizer(move(plan));
	context.profiler.EndPhase();

	return plan;
}

static idx_t CountOperatorType(LogicalOperator &plan, LogicalOperatorType type) {
	idx_t count = plan.type == type;
	for (auto &child : plan.children) {
		count += CountOperatorType(*child, type);
	}
	return count;
}

unique_ptr<LogicalOperator> ReOptimizer::AlgorithmNStep(idx_t n, unique_ptr<LogicalOperator> plan,
															const string temporary_table_name) {
	// leave at least 3 tables remaining after an iteration, else do the rest in 1 go
	idx_t gets = CountOperatorType(*plan, LogicalOperatorType::GET);
	if (gets - (n - 1) < 3) {
		done = true;
		return plan;
	}

	vector<LogicalOperator *> joins = ExtractJoinOperators(*plan);
	idx_t chosen_index = joins.size();
	idx_t chosen_gets = gets;
	for (idx_t i = 1; i < joins.size(); i++) {
		LogicalOperator *join = joins[i];
		idx_t join_gets = CountOperatorType(*join, LogicalOperatorType::GET);
		if (join_gets == n) {
			// if we find a join with n tables, we go for it
			chosen_index = i;
			break;
		}
		// keep track of join with lowest # of tables in case we dont find a join with n tables
		if (join_gets < chosen_gets && join_gets >= n) {
			chosen_index = i;
			chosen_gets = join_gets;
		}
	}
	// the join we chose will leave less than 3 tables remaining, just execute the rest
	if (gets - (chosen_gets - 1) < 3) {
		done = true;
		return plan;
	}

	// we were able to find a suitable join
	plan = PerformPartialPlan(move(plan), joins[chosen_index], temporary_table_name);

	context.profiler.StartPhase("optimizer");
	plan = CallOptimizer(move(plan));
	context.profiler.EndPhase();

	return plan;
}

unique_ptr<LogicalOperator> ReOptimizer::AlgorithmSmartStep(idx_t n, unique_ptr<LogicalOperator> plan,
															const string temporary_table_name) {
	// if there are less than 'n' risky operators, execute the rest of the plan
	idx_t total_risky_count = plan->RiskyOperatorCount();
	if (total_risky_count <= n) {
		done = true;
		return plan;
	}
	// find the child of a _somewhat_ risky join
	vector<LogicalOperator *> joins = ExtractJoinOperators(*plan);
	idx_t chosen_parent_risky_count = total_risky_count;
	idx_t chosen_child_risky_count = 0;
	idx_t chosen_child_cardinality = 0;
	LogicalOperator *chosen_child = joins[0];
	for (idx_t i = 1; i < joins.size(); i++) {
		LogicalOperator *join = joins[i];
		
		// find the operator with the lowest risky operator count greater than 'n'
		idx_t join_risky_count = join->RiskyOperatorCount();
		if (join_risky_count <= n || join_risky_count > chosen_parent_risky_count)
			continue;

		// reset values if we find a lower risky count (those are for comparing with equal count)
		if (join_risky_count < chosen_parent_risky_count) {
			chosen_parent_risky_count = join_risky_count;
			chosen_child_risky_count = 0;
			chosen_child_cardinality = 0;
		}

		// find riskiest child with highest estimated cardinality
		for (auto &child : join->children) {
			idx_t child_risk = child->RiskyOperatorCount();

			// reset value if we find a higher risky count
			if (child_risk > chosen_child_risky_count)
				chosen_child_cardinality = 0;

			if (child_risk >= chosen_child_risky_count) {
				chosen_child_risky_count = child_risk;
				idx_t child_cardinality = child->EstimateCardinality();
				if (child_cardinality >= chosen_child_cardinality) {
					chosen_child_cardinality = child_cardinality;
					chosen_child = child.get();
				}
			}
		}
	}
	// no suitable join was found, execute the rest
	if (CountOperatorType(*plan, LogicalOperatorType::GET) - (CountOperatorType(*chosen_child, LogicalOperatorType::GET) - 1) < 3) {
		done = true;
		return plan;
	}
	// we were able to find a suitable join
	plan = PerformPartialPlan(move(plan), chosen_child, temporary_table_name);

	context.profiler.StartPhase("optimizer");
	plan = CallOptimizer(move(plan));
	context.profiler.EndPhase();

	return plan;
}

static bool QErrorOverThreshold(idx_t x, idx_t y, idx_t thresh) {
	if (std::min(x, y) == 0)
		return true;
	return std::max(x, y) / std::min(x, y) > thresh;
}

void ReOptimizer::SetTrueCardinality(LogicalOperator &plan, LogicalOperator &subquery_plan) {
	if (plan.children.empty())
		return;
	// search children
	for (idx_t child_i = 0; child_i < plan.children.size(); child_i++) {
		auto &child = plan.children[child_i];
		if (!(child->type == subquery_plan.type && child->ParamsToString() == subquery_plan.ParamsToString()))
			continue;
		// we found the operator, set cardinality
		plan.children[child_i]->true_cardinality = GetTrueCardinality(*child);
		return;
	}
	// enter recursion until the operator is found
	for (auto &child : plan.children)
		SetTrueCardinality(*child.get(), subquery_plan);
}

unique_ptr<LogicalOperator> ReOptimizer::SimulatedReOptimize(unique_ptr<LogicalOperator> plan, const string query) {
	// compute_cost = true;
	if (compute_cost) {
		binding_name_mapping.clear();
		FindAliases(*plan);
		plan_cost += GetTrueCost(*plan);
	}

	idx_t minimum_remaining_plan_size = 2;
	idx_t q_error_threshold = 32;

	const string tablename_prefix = "_reopt_temp_" + to_string(hash<string>{}(query));

	for (int iter = 0; true; iter++) {
		binding_name_mapping.clear();
		FindAliases(*plan);

		const string temp_table_name = tablename_prefix + "_" + to_string(iter);
		vector<LogicalOperator *> nodes = ExtractJoinOperators(*plan);
		vector<LogicalOperator *> filters = ExtractFilterOperators(*plan);
		nodes.insert(nodes.end(), filters.begin(), filters.end());

		while (nodes.size() > minimum_remaining_plan_size) {
			LogicalOperator *node = nodes.back();
			nodes.pop_back();
			context.profiler.StartPhase("collect_statistics");
			SetTrueCardinality(*plan, *node);
			context.profiler.EndPhase();
			
			if (QErrorOverThreshold(node->true_cardinality, node->EstimateCardinality(), q_error_threshold)) {
				plan = PerformPartialPlan(move(plan), node, temp_table_name);

				context.profiler.StartPhase("optimizer");
				plan = CallOptimizer(move(plan));
				context.profiler.EndPhase();
				break;
			}
		}

		if (nodes.size() <= minimum_remaining_plan_size)
			break;
	}
	if (compute_cost) {
		binding_name_mapping.clear();
		FindAliases(*plan);
		plan_cost += GetTrueCost(*plan);
		Printer::Print(to_string(plan_cost));
	}
	// Printer::Print(to_string(materialize_size));
	return plan;
}

static bool QErrorOverThreshold(idx_t x, idx_t y, double thresh) {
	if (std::min(x, y) == 0)
		return true;
	return ((double) std::max(x, y)) / ((double) std::min(x, y)) > thresh;
}

idx_t ReOptimizer::RemainingCost(LogicalOperator &plan, LogicalOperator &subquery_plan, idx_t true_cardinality) {
	if (plan.children.empty())
		return 0;
	if (plan.type == subquery_plan.type && plan.ParamsToString() == subquery_plan.ParamsToString()) {
		// we found the operator, return true cardinality, skippings its children
		return true_cardinality;
	}
	idx_t result = 0;
	if (plan.type == LogicalOperatorType::COMPARISON_JOIN)
		result += plan.EstimateCardinality();
	for (auto &child : plan.children)
		result += RemainingCost(*child.get(), subquery_plan, true_cardinality);
	return result;
}

unique_ptr<LogicalOperator> ReOptimizer::SimulatedReOptimizeCost(unique_ptr<LogicalOperator> plan, const string query, double thresh) {
	// compute_cost = true;
	idx_t minimum_remaining_plan_size = 2;

	const string tablename_prefix = "_reopt_temp_" + to_string(hash<string>{}(query));
	for (int iter = 0; true; iter++) {
		binding_name_mapping.clear();
		FindAliases(*plan);

		// loop over 'risky' nodes
		vector<LogicalOperator *> nodes = ExtractJoinOperators(*plan);
		vector<LogicalOperator *> filters = ExtractFilterOperators(*plan);
		nodes.insert(nodes.end(), filters.begin(), filters.end());
		while (nodes.size() > minimum_remaining_plan_size) {
			LogicalOperator *node = nodes.back();
			nodes.pop_back();

			// save estimated cost of remaining plan (if 'node' was materialized)
			idx_t estimated_plan_cost = RemainingCost(*plan, *node, node->EstimateCardinality());

			// set true cardinalities of the node and its children
			context.profiler.StartPhase("collect_statistics");
			idx_t true_cardinality = GetTrueCardinality(*node);
			context.profiler.EndPhase();

			// re-compute estimated cost of remaining plan given the true cardinality of 'node'
			idx_t updated_plan_cost = RemainingCost(*plan, *node, true_cardinality);
			
			if (QErrorOverThreshold(updated_plan_cost, estimated_plan_cost, thresh)) {
				const string temp_table_name = tablename_prefix + "_" + to_string(iter);
				plan = PerformPartialPlan(move(plan), node, temp_table_name);

				context.profiler.StartPhase("optimizer");
				plan = CallOptimizer(move(plan));
				context.profiler.EndPhase();
				break;
			}
		}

		if (nodes.size() <= minimum_remaining_plan_size)
			break;
	}
	if (compute_cost) {
		binding_name_mapping.clear();
		FindAliases(*plan);
		plan_cost += GetTrueCost(*plan);
		Printer::Print(to_string(plan_cost));
	}
	// Printer::Print(to_string(materialize_size));
	return plan;
}

idx_t ReOptimizer::GetTrueCardinality(LogicalOperator &subquery_plan) {
	vector<string> queried_tables;
	vector<string> where_conditions;
	string subquery = "SELECT COUNT(*) FROM (" + CreateSubQuery(subquery_plan, queried_tables, where_conditions, false) + ") AS subquery;";
	auto result = ExecuteSubQuery(subquery, false);
	MaterializedQueryResult *mqr = static_cast<MaterializedQueryResult *>(result.get());
	Value count = mqr->collection.GetValue(0, 0);
	// Printer::Print(subquery_plan.ParamsToString() + " - " + count.ToString());
	return stoi(count.ToString());
}

unique_ptr<LogicalOperator> ReOptimizer::PerformPartialPlan(unique_ptr<LogicalOperator> plan,
															LogicalOperator *subquery_plan,
                                                            const string temporary_table_name) {
	// Printer::Print("\n----------------------------- before");
	// plan->Print();
	// Printer::Print("-----------------------------\n");

	if (compute_cost) {
		binding_name_mapping.clear();
		FindAliases(*plan);
		plan_cost += GetTrueCost(*subquery_plan);
	}

	context.profiler.StartPhase("reopt_pre_tooling");
	plan = GenerateProjectionMaps(move(plan));
	binding_name_mapping.clear();
	FindAliases(*plan);
	context.profiler.EndPhase();

	context.profiler.StartPhase("subquery");
	vector<string> queried_tables;
	vector<string> where_conditions;
	string subquery = "CREATE TEMPORARY TABLE main." + temporary_table_name + " AS (" +
					  CreateSubQuery(*subquery_plan, queried_tables, where_conditions, true) + ");";
	ExecuteSubQuery(subquery, true);
	context.profiler.EndPhase();

	// materialize_size += GetTable("main", temporary_table_name)->storage->info->cardinality;

	// Printer::Print(subquery);

	context.profiler.StartPhase("reopt_post_tooling");
	plan = AdjustPlan(move(plan), *subquery_plan, temporary_table_name);
	FixColumnBindings(*plan);
	rebind_mapping.clear();
	plan = ClearLeftProjectionMaps(move(plan));
	context.profiler.EndPhase();

	// Printer::Print("\n----------------------------- after");
	// plan->Print();
	// Printer::Print("-----------------------------\n\n");

	return plan;
}

idx_t ReOptimizer::GetTrueCost(LogicalOperator &plan) {
	idx_t cost = 0;
	vector<LogicalOperator *> joins = ExtractJoinOperators(plan);
	for (LogicalOperator *join : joins) {
		cost += GetTrueCardinality(*join);
	}
	return cost;
}

vector<LogicalOperator *> ReOptimizer::ExtractJoinOperators(LogicalOperator &plan) {
	// the deepest join operators appear at the end of the vector
	vector<LogicalOperator *> joins;
	if (plan.children.empty())
		return joins;
	bool recurse = true;
	switch (plan.type) {
	case LogicalOperatorType::FILTER:
		// special case: IN operator, has a join down the line, but should not be counted
		for (auto &expr : plan.expressions) {
			if (expr->type == ExpressionType::COMPARE_IN) {
				recurse = false;
				break;
			}
		}
		break;
	// case LogicalOperatorType::JOIN:
	// case LogicalOperatorType::ANY_JOIN:
	// case LogicalOperatorType::DELIM_JOIN:
	// case LogicalOperatorType::CROSS_PRODUCT:
	case LogicalOperatorType::COMPARISON_JOIN: // Pretty sure we only get COMPARISON_JOIN in JOB
		if (plan.children[1]->type == LogicalOperatorType::CHUNK_GET) {
			recurse = false;
			break;
		}
		joins.push_back(&plan);
		break;
	default:
		// fall through
		break;
	}
	if (recurse) {
		for (auto &child : plan.children) {
			vector<LogicalOperator *> child_joins = ExtractJoinOperators(*child);
			joins.insert(joins.end(), child_joins.begin(), child_joins.end());
		}
	}
	return joins;
}

vector<LogicalOperator *> ReOptimizer::ExtractFilterOperators(LogicalOperator &plan) {
	vector<LogicalOperator *> filters;
	if (plan.children.empty())
		return filters;
	if (plan.type == LogicalOperatorType::FILTER) {
		if (!plan.expressions.empty())
			filters.push_back(&plan);
		return filters;
	} else if (plan.type == LogicalOperatorType::COMPARISON_JOIN && plan.children[1]->type == LogicalOperatorType::CHUNK_GET) {
		filters.push_back(&plan);
		return filters;
	}
	for (auto &child : plan.children) {
		vector<LogicalOperator *> child_filters = ExtractFilterOperators(*child);
		filters.insert(filters.end(), child_filters.begin(), child_filters.end());
	}
	return filters;
}

//! Gets column bindings from a join, but assumes projection maps are actually filled
static vector<ColumnBinding> GetColumnBindings(LogicalComparisonJoin &join) {
	vector<ColumnBinding> cbs;
	if (!join.left_projection_map.empty()) {
		vector<ColumnBinding> left_cbs =
		    LogicalOperator::MapBindings(join.children[0]->GetColumnBindings(), join.left_projection_map);
		for (auto cb : left_cbs)
			cbs.push_back(cb);
	}
	if (!join.right_projection_map.empty()) {
		vector<ColumnBinding> right_cbs =
		    LogicalOperator::MapBindings(join.children[1]->GetColumnBindings(), join.right_projection_map);
		for (auto cb : right_cbs)
			cbs.push_back(cb);
	}
	return cbs;
}

unique_ptr<LogicalOperator> ReOptimizer::GenerateProjectionMaps(unique_ptr<LogicalOperator> plan) {
	if (CountOperatorType(*plan, LogicalOperatorType::COMPARISON_JOIN) == 0)
		return plan;
	// store the column bindings that are needed by the parents of the child join
	vector<ColumnBinding> column_bindings;
	if (plan->type == LogicalOperatorType::COMPARISON_JOIN &&
	    (plan->children[0]->type == LogicalOperatorType::COMPARISON_JOIN ||
	     plan->children[1]->type == LogicalOperatorType::COMPARISON_JOIN)) {
		auto *join = static_cast<LogicalComparisonJoin *>(plan.get());
		for (idx_t cond_i = 0; cond_i < join->conditions.size(); cond_i++) {
			// from conditions
			JoinCondition &jc = join->conditions[cond_i];
			auto l = (BoundColumnRefExpression &)*jc.left.get();
			auto r = (BoundColumnRefExpression &)*jc.right.get();
			column_bindings.push_back(l.binding);
			column_bindings.push_back(r.binding);
		}
		// from bindings of this plan
		for (auto new_cb : GetColumnBindings(*join)) {
			bool found = false;
			for (auto existing_cb : column_bindings) {
				if (existing_cb == new_cb) {
					found = true;
					break;
				}
			}
			if (!found)
				column_bindings.push_back(new_cb);
		}

		// update child join projection maps
		for (idx_t i = 0; i < 2; i++) {
			if (plan->children[i]->type != LogicalOperatorType::COMPARISON_JOIN)
				continue;
			// if (plan->children[i]->children[1]->type == LogicalOperatorType::CHUNK_GET)
			// 	continue;
			auto *child_join = static_cast<LogicalComparisonJoin *>(join->children[i].get());
			vector<ColumnBinding> child_left_cbs = LogicalOperator::MapBindings(child_join->children[0]->GetColumnBindings(), child_join->left_projection_map);
			idx_t n_kept = 0;
			vector<idx_t> removed_columns;
			for (idx_t cb_i = 0; cb_i < child_left_cbs.size(); cb_i++) {
				bool keep = false;
				for (auto cb : column_bindings) {
					if (child_left_cbs[cb_i] == cb) {
						keep = true;
						n_kept++;
						break;
					}
				}
				if (keep) {
					if (find(child_join->left_projection_map.begin(), child_join->left_projection_map.end(), cb_i) == child_join->left_projection_map.end()) {
						child_join->left_projection_map.push_back(cb_i);
					}
				} else {
					removed_columns.push_back(cb_i);
				}
			}

			// Need to only remove columns that were kept in the first place!
			// somewhat of a hack fix ...
			if (n_kept == 0)
				continue;

			// update left/right projection maps of this join accordingly
			if (i == 0) {
				for (idx_t lpi = 0; lpi < join->left_projection_map.size(); lpi++) {
					idx_t decr = 0;
					for (idx_t removed_col : removed_columns) {
						if (removed_col < join->left_projection_map[lpi]) {
							decr++;
						}
					}
					join->left_projection_map[lpi] -= decr;
				}
			} else {
				for (idx_t rpi = 0; rpi < join->right_projection_map.size(); rpi++) {
					idx_t decr = 0;
					for (idx_t removed_col : removed_columns) {
						if (removed_col < join->right_projection_map[rpi]) {
							decr++;
						}
					}
					join->right_projection_map[rpi] -= decr;
				}
			}
		}
	} else if (plan->children[0]->type == LogicalOperatorType::COMPARISON_JOIN) {
		for (idx_t expr_i = 0; expr_i < plan->expressions.size(); expr_i++) {
			auto &expr = *plan->expressions[expr_i];
			if (expr.GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF)
				continue;
			auto bcre = (BoundColumnRefExpression &)expr;
			column_bindings.push_back(bcre.binding);
		}
		switch (plan->type) {
		case LogicalOperatorType::AGGREGATE_AND_GROUP_BY: {
			for (idx_t i = 0; i < plan->expressions.size(); i++) {
				auto &bae = (BoundAggregateExpression &)*plan->expressions[i];
				for (auto &child : bae.children) {
					auto bcre = (BoundColumnRefExpression &)*child.get();
					column_bindings.push_back(bcre.binding);
				}
			}
			break;
		}
		case LogicalOperatorType::FILTER: {
			for (auto cb : plan->GetColumnBindings()) {
				column_bindings.push_back(cb);
			}
			break;
		}
		default:
			// fall through
			break;
		}
		// get the child join, modify its left projection map
		LogicalComparisonJoin *child_join = static_cast<LogicalComparisonJoin *>(plan->children[0].get());
		auto left_cbs = child_join->children[0]->GetColumnBindings();
		// add index of needed bindings
		for (idx_t cb_i = 0; cb_i < left_cbs.size(); cb_i++) {
			for (auto cb : column_bindings) {
				if (left_cbs[cb_i] == cb) {
					child_join->left_projection_map.push_back(cb_i);
					break;
				}
			}
		}
	}
	// recursively propagate operator tree
	for (idx_t child_i = 0; child_i < plan->children.size(); child_i++) {
		if (plan->children[child_i]->type != LogicalOperatorType::FILTER)
			plan->children[child_i] = GenerateProjectionMaps(move(plan->children[child_i]));
	}
	return plan;
}

void ReOptimizer::FindAliases(LogicalOperator &plan) {
	if (plan.type == LogicalOperatorType::GET) {
		auto *get = static_cast<LogicalGet *>(&plan);
		for (idx_t i = 0; i < get->column_ids.size(); i++) {
			binding_name_mapping["#[" + std::to_string(get->table_index) + "." + std::to_string(i) + "]"] =
			    get->table->columns[get->column_ids[i]].name;
		}
	}
	// recursively propagate operator tree
	for (auto &child : plan.children) {
		FindAliases(*child);
	}
}

string ReOptimizer::CreateSubQuery(LogicalOperator &plan, vector<string> &queried_tables, vector<string> &where_conditions, bool proj_map_filled) {
	switch (plan.type) {
	case LogicalOperatorType::COMPARISON_JOIN: {
		LogicalComparisonJoin *join = static_cast<LogicalComparisonJoin *>(&plan);
		if (plan.children[1]->type == LogicalOperatorType::CHUNK_GET) {
			// special case: IN operator
			JoinCondition &join_condition = join->conditions[0];
			auto l_bind = ((BoundColumnRefExpression &)*join_condition.left.get()).binding;
			LogicalChunkGet *chunk_get = static_cast<LogicalChunkGet *>(join->children[1].get());
			idx_t count = chunk_get->collection->chunks[0]->size();
			Vector *vals = &chunk_get->collection->chunks[0]->data[0];
			bool is_string = vals->type == TypeId::VARCHAR;
			vector<string> in_vals;
			for (idx_t i = 0; i < count; i++) {
				if (is_string)
					in_vals.push_back("'" + vals->GetValue(i).ToString() + "'");
				else
					in_vals.push_back(vals->GetValue(i).ToString());
			}
			where_conditions.push_back("t" + to_string(l_bind.table_index) + "." + binding_name_mapping[l_bind.ToString()] + " IN (" +
			                           JoinStrings(in_vals, ",") + ")");
		} else {
			for (idx_t cond_i = 0; cond_i < join->conditions.size(); cond_i++) {
				JoinCondition &join_condition = join->conditions[cond_i];
				auto l_bind = ((BoundColumnRefExpression &)*join_condition.left.get()).binding;
				auto r_bind = ((BoundColumnRefExpression &)*join_condition.right.get()).binding;
				where_conditions.push_back("t" + to_string(l_bind.table_index) + "." + binding_name_mapping[l_bind.ToString()] + " = " +
				                           "t" + to_string(r_bind.table_index) + "." + binding_name_mapping[r_bind.ToString()]);
			}
		}
		break;
	}
	case LogicalOperatorType::FILTER: {
		LogicalFilter *filter = static_cast<LogicalFilter *>(&plan);
		auto filter_conditions = GetFilterStrings(filter);
		where_conditions.insert(where_conditions.end(), filter_conditions.begin(), filter_conditions.end());
		break;
	}
	case LogicalOperatorType::GET: {
		LogicalGet *get = static_cast<LogicalGet *>(&plan);
		queried_tables.push_back(get->table->schema->name + "." + get->table->name + " AS t" +
		                         to_string(get->table_index));
		break;
	}
	case LogicalOperatorType::CHUNK_GET: {
		// handled by special COMPARISON_JOIN case above
		break;
	}
	default:
		Printer::Print("Exception 1 " + LogicalOperatorToString(plan.type));
		throw new ReOptimizerException("Unexpected operator in query plan: '%s'", LogicalOperatorToString(plan.type));
	}
	// recursively propagate operator tree to fill vectors
	for (auto &child : plan.children) {
		CreateSubQuery(*child, queried_tables, where_conditions, proj_map_filled);
	}
	// grab selected columns
	vector<ColumnBinding> selected_column_bindings;
	if (plan.type == LogicalOperatorType::COMPARISON_JOIN) {
		LogicalComparisonJoin *join = static_cast<LogicalComparisonJoin *>(&plan);
		if (proj_map_filled)
			selected_column_bindings = GetColumnBindings(*join);
		else
			selected_column_bindings = join->GetColumnBindings();
	} else {
		selected_column_bindings = plan.GetColumnBindings();
	}
	vector<string> selected_columns;
	for (auto cb : selected_column_bindings)
		selected_columns.push_back("t" + to_string(cb.table_index) + "." + binding_name_mapping[cb.ToString()] + " AS " + binding_name_mapping[cb.ToString()] + to_string(cb.table_index));

	// create and return query
	return "SELECT " + JoinStrings(selected_columns, ", ") + " " +
		   "FROM " + JoinStrings(queried_tables, ", ") + " " +
		   "WHERE " + JoinStrings(where_conditions, " AND ");
}

vector<string> ReOptimizer::GetFilterStrings(LogicalFilter *filter) {
	vector<string> conditions;
	for (auto &expr : filter->expressions) {
		string expr_string = GetExpressionString(expr.get());
		if (expr_string != "")
			conditions.push_back(expr_string);
	}
	return conditions;
}

string ReOptimizer::GetExpressionString(Expression *expr) {
	switch (expr->GetExpressionClass()) {
	case ExpressionClass::BOUND_COLUMN_REF:
		return "";
	case ExpressionClass::BOUND_OPERATOR: {
		// special case: IN operator, has a join down the line, processed differently
		switch (expr->type) {
		case ExpressionType::COMPARE_IN:
			return "";
		case ExpressionType::OPERATOR_IS_NOT_NULL: {
			auto *op = static_cast<BoundOperatorExpression *>(expr);
			auto binding = static_cast<BoundColumnRefExpression *>(op->children[0].get())->binding;
			return "t" + to_string(binding.table_index) + "." + binding_name_mapping[binding.ToString()] + " IS NOT NULL";
		}
		case ExpressionType::OPERATOR_IS_NULL: {
			auto *op = static_cast<BoundOperatorExpression *>(expr);
			auto binding = static_cast<BoundColumnRefExpression *>(op->children[0].get())->binding;
			return "t" + to_string(binding.table_index) + "." + binding_name_mapping[binding.ToString()] + " IS NULL";
		}
		default:
			Printer::Print("Exception 4: " + expr->ToString() + " - " + ExpressionTypeToString(expr->type));
			throw new ReOptimizerException(
			    "Expected expression type of class BOUND_OPERATOR to be COMPARE_IN or IS_NOT_NULL, got '%s'",
			    ExpressionTypeToString(expr->type));
		}
	}
	case ExpressionClass::BOUND_COMPARISON: {
		auto *comparison = static_cast<BoundComparisonExpression *>(expr);
		return GetBoundComparisonString(comparison);
	}
	case ExpressionClass::BOUND_FUNCTION: {
		auto *func = static_cast<BoundFunctionExpression *>(expr);
		return GetBoundFunctionString(func);
	}
	case ExpressionClass::BOUND_CONJUNCTION: {
		auto *conjunction = static_cast<BoundConjunctionExpression *>(expr);
		vector<string> child_conditions;
		for (auto &child : conjunction->children) {
			child_conditions.push_back(GetExpressionString(child.get()));
		}
		return "(" +
		       JoinStrings(child_conditions, " " + ExpressionTypeToOperator(conjunction->GetExpressionType()) + " ") +
		       ")";
	}
	case ExpressionClass::BOUND_BETWEEN: {
		auto *between = static_cast<BoundBetweenExpression *>(expr);
		auto binding = static_cast<BoundColumnRefExpression *>(between->input.get())->binding;
		string condition = "t" + to_string(binding.table_index) + "." + binding_name_mapping[binding.ToString()] + " BETWEEN ";
		auto lower = static_cast<BoundConstantExpression *>(between->lower.get());
		auto upper = static_cast<BoundConstantExpression *>(between->upper.get());
		if (lower->value.type == TypeId::VARCHAR)
			condition += "'" + lower->value.ToString() + "'";
		else
			condition += lower->value.ToString();
		condition += " AND ";
		if (upper->value.type == TypeId::VARCHAR)
			condition += "'" + upper->value.ToString() + "'";
		else
			condition += upper->value.ToString();
		return condition;
	}
	default:
		Printer::Print("Exception 2: " + expr->ToString() + " - " +
		               ExpressionClassToString(expr->GetExpressionClass()) + ", " + ExpressionTypeToString(expr->type));
		throw new ReOptimizerException(
		    "Expected filter class to be BOUND_COMPARISON, BOUND_FUNCTION or BOUND_CONJUNCTION, got '%s'",
		    ExpressionClassToString(expr->GetExpressionClass()));
	}
}

string ReOptimizer::GetBoundComparisonString(BoundComparisonExpression *comparison) {
	// filter conditions - table reference is 'always' placed on the left, constant on the right
	auto binding = static_cast<BoundColumnRefExpression *>(comparison->left.get())->binding;
	string condition = "t" + to_string(binding.table_index) + "." + binding_name_mapping[binding.ToString()] + " " +
	                   ExpressionTypeToOperator(comparison->type) + " ";
	if (comparison->right->return_type == TypeId::VARCHAR) {
		// strings need to be escaped with quotes
		condition += "'" + comparison->right->ToString() + "'";
	} else {
		condition += comparison->right->ToString();
	}
	return condition;
}

string ReOptimizer::GetBoundFunctionString(BoundFunctionExpression *func) {
	// filter conditions - table reference is 'always' placed on the left, constant on the right
	auto binding = static_cast<BoundColumnRefExpression *>(func->children[0].get())->binding;
	string condition = "t" + to_string(binding.table_index) + "." + binding_name_mapping[binding.ToString()];
	if (func->function.name == "!~~")
		condition += " NOT LIKE ";
	else if (func->function.name == "contains" || func->function.name == "~~" || func->function.name == "prefix")
		condition += " LIKE ";
	if (func->function.name == "contains") {
		// % have gone missing with this expression
		condition += "'%" + func->children[1]->ToString() + "%'";
	} else if (func->function.name == "prefix") {
		// trailing % has gone missing
		condition += "'" + func->children[1]->ToString() + "%'";
	} else if (func->children[0]->return_type == TypeId::VARCHAR) {
		// strings need to be escaped with quotes
		condition += "'" + func->children[1]->ToString() + "'";
	} else {
		condition += func->children[1]->ToString();
	}
	return condition;
}

//! Extracts table names from GET operators
static vector<string> GetRelationSet(LogicalOperator &plan) {
	vector<string> relations;
	if (plan.type == LogicalOperatorType::GET) {
		auto *get = static_cast<LogicalGet *>(&plan);
		relations.push_back(get->table->schema->name + "." + get->table->name);
	} else {
		for (auto &child : plan.children) {
			vector<string> child_indices = GetRelationSet(*child);
			relations.insert(relations.end(), child_indices.begin(), child_indices.end());
		}
	}
	sort(relations.begin(), relations.end());
	return relations;
}

unique_ptr<LogicalOperator> ReOptimizer::AdjustPlan(unique_ptr<LogicalOperator> plan, LogicalOperator &old_op,
                                                    const string temporary_table_name) {
	TableCatalogEntry *table = GetTable("main", temporary_table_name);
	// Create a LogicalGet for the newly made temporary table (empty column_ids - filled in "ReplaceLogicalOperator")
	unique_ptr<LogicalGet> temp_table_get = make_unique<LogicalGet>(table, 0, vector<idx_t>());
	// replace 'subplan' with 'temp_table_get' in 'plan'
	rebind_mapping = {}; // reset
	ReplaceLogicalOperator(*plan, old_op, table);
	return plan;
}

TableCatalogEntry *ReOptimizer::GetTable(string schema, string table_name) {
	// Catalog::GetTable can only be called if there is an active transaction - else segfault
	CatalogEntry *entry;
	TableCatalogEntry *table;
	if (!context.transaction.HasActiveTransaction()) {
		context.transaction.BeginTransaction();
		entry = context.catalog.GetEntry(context, CatalogType::TABLE, schema, table_name);
		context.transaction.Commit();
	} else {
		entry = context.catalog.GetEntry(context, CatalogType::TABLE, schema, table_name);
	}
	table = static_cast<TableCatalogEntry *>(entry);
	return table;
}

bool ReOptimizer::ReplaceLogicalOperator(LogicalOperator &plan, LogicalOperator &old_op, TableCatalogEntry *table,
                                         idx_t depth) {
	if (plan.children.empty())
		return false;
	// search children
	for (idx_t child_i = 0; child_i < plan.children.size(); child_i++) {
		auto &child = plan.children[child_i];
		if (!(child->type == old_op.type && child->ParamsToString() == old_op.ParamsToString()))
			continue;
		vector<ColumnBinding> old_op_cbs;
		if (child->type == LogicalOperatorType::COMPARISON_JOIN) {
			auto *join = static_cast<LogicalComparisonJoin *>(child.get());	
			old_op_cbs = GetColumnBindings(*join);
		} else {
			old_op_cbs = child->GetColumnBindings();
		}
		idx_t new_table_index = 0;
		for (auto cb : old_op_cbs) {
			if (cb.table_index > new_table_index)
				new_table_index = cb.table_index;
		}
		for (idx_t new_index = 0; new_index < old_op_cbs.size(); new_index++) {
			rebind_mapping[old_op_cbs[new_index].ToString()] = ColumnBinding(new_table_index, new_index);
		}
		vector<idx_t> column_ids;
		for (idx_t column_id = 0; column_id < table->columns.size(); column_id++)
			column_ids.push_back(column_id);
		unique_ptr<LogicalGet> replacement_operator = make_unique<LogicalGet>(table, new_table_index, column_ids);
		plan.children[child_i] = move(replacement_operator);
		return true;
	}
	// enter recursion until the operator is found
	bool recursion_done;
	for (auto &child : plan.children) {
		recursion_done = ReplaceLogicalOperator(*child.get(), old_op, table, depth + 1);
		if (recursion_done)
			break;
	}
	return recursion_done;
}

void ReOptimizer::FixColumnBindings(LogicalOperator &plan) {
	// fix BoundColumRefs in JoinConditions, aggregates, filters
	switch (plan.type) {
	case LogicalOperatorType::COMPARISON_JOIN: {
		LogicalComparisonJoin *join = static_cast<LogicalComparisonJoin *>(&plan);
		for (idx_t condition_i = 0; condition_i < join->conditions.size(); condition_i++) {
			JoinCondition &jc = join->conditions[condition_i];
			auto l = (BoundColumnRefExpression &)*jc.left.get();
			auto r = (BoundColumnRefExpression &)*jc.right.get();
			if (rebind_mapping.find(l.ToString()) != rebind_mapping.end()) {
				unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
				    l.alias, l.return_type, rebind_mapping[l.ToString()], l.depth);
				join->conditions[condition_i].left = move(fixed_bcre);
			}
			if (rebind_mapping.find(r.ToString()) != rebind_mapping.end()) {
				unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
				    r.alias, r.return_type, rebind_mapping[r.ToString()], r.depth);
				join->conditions[condition_i].right = move(fixed_bcre);
			}
		}
		break;
	}
	// case LogicalOperatorType::AGGREGATE_AND_GROUP_BY: {
	// 	for (idx_t expr_i = 0; expr_i < plan.expressions.size(); expr_i++) {
	// 		auto &bae = (BoundAggregateExpression &)*plan.expressions[expr_i];
	// 		for (idx_t child_i = 0; child_i < bae.children.size(); child_i++) {
	// 			auto &child = bae.children[child_i];
	// 			auto e = (BoundColumnRefExpression &)*child.get();
	// 			if (rebind_mapping.find	(e.ToString()) != rebind_mapping.end()) {
	// 				unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
	// 				    e.alias, e.return_type, rebind_mapping[e.ToString()], e.depth);
	// 				bae.children[child_i] = move(fixed_bcre);
	// 			}
	// 		}
	// 	}
	// 	break;
	// }
	default:
		// fix other BoundColumnRefExpressions found in expressions (for e.g. LogicalProjection)
		for (idx_t expr_i = 0; expr_i < plan.expressions.size(); expr_i++) {
			auto &expr = *plan.expressions[expr_i];
			switch (expr.GetExpressionClass()) {
			case ExpressionClass::BOUND_COLUMN_REF: {
				auto bcre = (BoundColumnRefExpression &)expr;
				if (rebind_mapping.find(bcre.ToString()) != rebind_mapping.end()) {
					unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
						bcre.alias, bcre.return_type, rebind_mapping[bcre.ToString()], bcre.depth);
					plan.expressions[expr_i] = move(fixed_bcre);
				}
				break;
			}
			default:
				FixColumnBindings(&expr);
			}
		}
	}
	
	// recursively propagate operator tree
	for (auto &child : plan.children) {
		FixColumnBindings(*child);
	}
}

void ReOptimizer::FixColumnBindings(Expression *expr) {
	switch (expr->GetExpressionClass()) {
	case ExpressionClass::BOUND_OPERATOR: {
		auto *op = static_cast<BoundOperatorExpression *>(expr);
		for (idx_t child_i = 0; child_i < op->children.size(); child_i++) {
			if (op->children[child_i]->type == ExpressionType::BOUND_COLUMN_REF) {
				auto *bcre = static_cast<BoundColumnRefExpression *>(op->children[child_i].get());
				if (rebind_mapping.find(bcre->ToString()) != rebind_mapping.end()) {
					unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(bcre->alias, bcre->return_type, rebind_mapping[bcre->ToString()], bcre->depth);
					op->children[child_i] = move(fixed_bcre);
				}
			} else {
				auto &child = *op->children[child_i];
				FixColumnBindings(&child);
			}
		}
		break;
	}
	case ExpressionClass::BOUND_COMPARISON: {
		auto *comparison = static_cast<BoundComparisonExpression *>(expr);
		auto *bcre = static_cast<BoundColumnRefExpression *>(comparison->left.get());
		if (rebind_mapping.find(bcre->ToString()) != rebind_mapping.end()) {
			unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
			    bcre->alias, bcre->return_type, rebind_mapping[bcre->ToString()], bcre->depth);
			comparison->left = move(fixed_bcre);
		}
		break;
	}
	case ExpressionClass::BOUND_FUNCTION: {
		auto *func = static_cast<BoundFunctionExpression *>(expr);
		auto *bcre = static_cast<BoundColumnRefExpression *>(func->children[0].get());
		if (rebind_mapping.find(bcre->ToString()) != rebind_mapping.end()) {
			unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
			    bcre->alias, bcre->return_type, rebind_mapping[bcre->ToString()], bcre->depth);
			func->children[0] = move(fixed_bcre);
		}
		break;
	}
	case ExpressionClass::BOUND_CONJUNCTION: {
		auto *conjunction = static_cast<BoundConjunctionExpression *>(expr);
		for (idx_t child_i = 0; child_i < conjunction->children.size(); child_i++) {
			auto &expr = *conjunction->children[child_i];
			FixColumnBindings(&expr);
		}
		break;
	}
	case ExpressionClass::BOUND_BETWEEN: {
		auto *between = static_cast<BoundBetweenExpression *>(expr);
		auto bcre = static_cast<BoundColumnRefExpression *>(between->input.get());
		if (rebind_mapping.find(bcre->ToString()) != rebind_mapping.end()) {
			unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
			    bcre->alias, bcre->return_type, rebind_mapping[bcre->ToString()], bcre->depth);
			between->input = move(fixed_bcre);
		}
		break;
	}
	case ExpressionClass::BOUND_AGGREGATE: {
		auto *aggregate = static_cast<BoundAggregateExpression *>(expr);
		for (idx_t child_i = 0; child_i < aggregate->children.size(); child_i++) {
			auto *bcre = static_cast<BoundColumnRefExpression *>(aggregate->children[child_i].get());
			if (rebind_mapping.find(bcre->ToString()) != rebind_mapping.end()) {
				unique_ptr<BoundColumnRefExpression> fixed_bcre = make_unique<BoundColumnRefExpression>(
				    bcre->alias, bcre->return_type, rebind_mapping[bcre->ToString()], bcre->depth);
				aggregate->children[child_i] = move(fixed_bcre);
			}
		}
		break;
	}
	default:
		Printer::Print("Exception 3: " + expr->ToString() + " - " +
		               ExpressionClassToString(expr->GetExpressionClass()) + ", " + ExpressionTypeToString(expr->type));
		throw new ReOptimizerException(
		    "Expected filter class to be BOUND_COMPARISON, BOUND_FUNCTION or BOUND_CONJUNCTION, got '%s'",
		    ExpressionClassToString(expr->GetExpressionClass()));
	}
}

unique_ptr<QueryResult> ReOptimizer::ExecuteSubQuery(const string subquery, bool enable_profiler) {
	context.subquery = true;
	context.transaction.Commit();
	// store the state of the profiler
	bool profiler_was_enabled = context.profiler.IsEnabled();
	if (!enable_profiler)
		context.profiler.Disable();
	auto result = context.QueryWithoutLock(subquery, false);
	// restore the state of the profiler
	context.transaction.BeginTransaction();
	if (profiler_was_enabled)
		context.profiler.Enable(); 
	context.subquery = false;
	return result;
}

unique_ptr<LogicalOperator> ReOptimizer::CallOptimizer(unique_ptr<LogicalOperator> plan) {
	// not interested in optimizer sub-phases
	bool enabled = context.profiler.IsEnabled();
	context.profiler.Disable();
	Optimizer optimizer(binder, context);
	plan = optimizer.Optimize(move(plan));
	if (enabled)
		context.profiler.Enable();
	return plan;
}

unique_ptr<LogicalOperator> ReOptimizer::ClearLeftProjectionMaps(unique_ptr<LogicalOperator> plan) {
	if (plan->children.empty())
		return plan;
	if (plan->type == LogicalOperatorType::COMPARISON_JOIN) {
		auto *join = static_cast<LogicalComparisonJoin *>(plan.get());
		join->left_projection_map.clear();
	}
	for (idx_t i = 0; i < plan->children.size(); i++) {
		plan->children[i] = ClearLeftProjectionMaps(move(plan->children[i]));
	}
	return plan;
}

string ReOptimizer::JoinStrings(vector<string> strings, string delimiter) {
	string joined_strings = "";
	for (idx_t i = 0; i < strings.size(); i++) {
		joined_strings += strings[i];
		if (i < strings.size() - 1)
			joined_strings += delimiter;
	}
	return joined_strings;
}
