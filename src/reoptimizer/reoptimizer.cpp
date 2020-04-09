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

// TODO: write tests that figure out if this shit actually works
ReOptimizer::ReOptimizer(ClientContext &context, Binder &binder) : context(context), binder(binder) {
}

unique_ptr<LogicalOperator> ReOptimizer::ReOptimize(unique_ptr<LogicalOperator> plan, const string query) {
	const string tablename_prefix = "_reopt_temp_" + to_string(hash<string>{}(query));
	// re-optimization loop
	for (int iter = 0; true; iter++) {
		// profiling for each iteration
		context.profiler.StartPhase(to_string(iter));
		// create and execute subquery, adjust plan accordingly
		const string temp_table_name = tablename_prefix + "_" + to_string(iter);
		plan = PerformPartialPlan(move(plan), temp_table_name);
		if (done) {
			context.profiler.EndPhase();
			break;
		}

		// TODO: if q-error (or some other criterion) is low we should not call optimizer
		context.profiler.StartPhase("optimizer");
		plan = CallOptimizer(move(plan));
		context.profiler.EndPhase();

		// if (ExtractJoinOperators(*plan).size() <= 1)
		// 	break;
		// end iteration profiling phase
		context.profiler.EndPhase();
	}
	context.profiler.StartPhase("clear_left_proj_maps");
	plan = ClearLeftProjectionMaps(move(plan));
	context.profiler.EndPhase();

	Printer::Print("OUTPUT");

	return plan;
}

unique_ptr<LogicalOperator> ReOptimizer::PerformPartialPlan(unique_ptr<LogicalOperator> plan,
                                                            const string temporary_table_name) {
	context.profiler.StartPhase("decide_subquery");
	auto *subquery_plan = DecideSubQueryPlan(*plan);
	context.profiler.EndPhase();

	if (done)
		return plan;
	// auto *join_subquery_plan = static_cast<LogicalComparisonJoin *>(subquery_plan);

	Printer::Print("\n----------------------------- before");
	plan->Print();
	Printer::Print("-----------------------------\n");

	context.profiler.StartPhase("generate_projection_maps");
	plan = GenerateProjectionMaps(move(plan));
	context.profiler.EndPhase();

	context.profiler.StartPhase("map_binding_names");
	bta.clear();
	CreateMaps(*plan);
	context.profiler.EndPhase();

	context.profiler.StartPhase("create_subquery");
	vector<string> queried_tables;
	vector<string> where_conditions;
	string subquery = CreateSubQuery(*subquery_plan, temporary_table_name, queried_tables, where_conditions);
	context.profiler.EndPhase();

	Printer::Print(subquery);

	context.profiler.StartPhase("execute_subquery");
	ExecuteSubQuery(subquery);
	context.profiler.EndPhase();

	// context.profiler.StartPhase("inject_cardinalities");
	// InjectCardinalities(*subquery_plan, temporary_table_name);
	// context.profiler.EndPhase();

	context.profiler.StartPhase("adjust_plan");
	plan = AdjustPlan(move(plan), *subquery_plan, temporary_table_name);
	context.profiler.EndPhase();

	context.profiler.StartPhase("fix_bindings");
	FixColumnBindings(*plan);
	rebind_mapping.clear();
	context.profiler.EndPhase();

	Printer::Print("\n----------------------------- after");
	plan->Print();
	Printer::Print("-----------------------------\n\n");

	return plan;
}

LogicalOperator *ReOptimizer::DecideSubQueryPlan(LogicalOperator &plan) {
	vector<LogicalOperator *> joins = ExtractJoinOperators(plan);
	if (joins.size() <= 1) {
		done = true;
		return &plan;
	}

	vector<LogicalOperator *> filters = ExtractFilterOperators(plan);
	if (!filters.empty()) {
		idx_t selected_filter_i = filters.size() - 1;
		return filters[selected_filter_i];
	}

	// TODO: implement selection procedure
	// select last one for now - deciding this is basically the whole re-optimization strategy
	// the selection procedure might also decide that we are done, instead of selecting an index
	idx_t selected_join_i = joins.size() - 1;
	return joins[selected_join_i];
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
			if (expr->type == ExpressionType::COMPARE_IN)
				recurse = false;
		}
		break;
	// case LogicalOperatorType::JOIN:
	// case LogicalOperatorType::ANY_JOIN:
	// case LogicalOperatorType::DELIM_JOIN:
	// case LogicalOperatorType::CROSS_PRODUCT:
	case LogicalOperatorType::COMPARISON_JOIN: // Pretty sure we only get COMPARISON_JOIN in JOB
		if (plan.children[1]->type != LogicalOperatorType::CHUNK_GET)
			joins.push_back(&plan);
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

static idx_t CountJoinOperators(LogicalOperator &plan) {
	idx_t count = plan.type == LogicalOperatorType::COMPARISON_JOIN;
	for (auto &child : plan.children) {
		count += CountJoinOperators(*child);
	}
	return count;
}

unique_ptr<LogicalOperator> ReOptimizer::GenerateProjectionMaps(unique_ptr<LogicalOperator> plan) {
	// idx_t n_child_joins = 0;
	// for (auto &child : plan->children) {
	// 	n_child_joins += ExtractJoinOperators(*child).size();
	// }
	// Printer::Print(plan->ParamsToString() + " - " + to_string(n_child_joins));
	if (CountJoinOperators(*plan) == 0)
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
			auto *child_join = static_cast<LogicalComparisonJoin *>(plan->children[i].get());
			vector<ColumnBinding> left_cbs = child_join->children[0]->GetColumnBindings();
			idx_t n_kept = 0;
			vector<idx_t> removed_columns;
			for (idx_t cb_i = 0; cb_i < left_cbs.size(); cb_i++) {
				bool keep = false;
				for (auto cb : column_bindings) {
					if (left_cbs[cb_i] == cb) {
						keep = true;
						n_kept++;
						break;
					}
				}
				if (keep) {
					// child_join->left_projection_map.push_back(cb_i);
					// Printer::Print("1: " + to_string(cb_i));
					if (find(child_join->left_projection_map.begin(), child_join->left_projection_map.end(), cb_i) == child_join->left_projection_map.end()) {
						child_join->left_projection_map.push_back(cb_i);
					}
				} else {
					removed_columns.push_back(cb_i);
				}
			}

			// FIXME: something wrong here ...
			// Need to only remove columns that were kept in the first place!
			if (n_kept == 0) // hack fix for this FIXME
				continue;

			// update left/right projection maps of this join accordingly
			if (i == 0) {
				for (idx_t rpi = 0; rpi < join->left_projection_map.size(); rpi++) {
					idx_t decr = 0;
					for (idx_t removed_col : removed_columns) {
						if (removed_col < join->left_projection_map[rpi]) {
							decr++;
						}
					}
					join->left_projection_map[rpi] -= decr;
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
			// default:
		}
		// get the child join, modify its left projection map
		LogicalComparisonJoin *child_join = child_join = static_cast<LogicalComparisonJoin *>(plan->children[0].get());
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

void ReOptimizer::CreateMaps(LogicalOperator &plan) {
	if (plan.type == LogicalOperatorType::GET) {
		auto *get = static_cast<LogicalGet *>(&plan);
		for (idx_t i = 0; i < get->column_ids.size(); i++) {
			bta["#[" + std::to_string(get->table_index) + "." + std::to_string(i) + "]"] =
			    get->table->columns[get->column_ids[i]].name;
		}
	}
	// recursively propagate operator tree
	for (auto &child : plan.children) {
		CreateMaps(*child);
	}
}

string ReOptimizer::CreateSubQuery(LogicalOperator &plan, const string temporary_table_name,
                                   vector<string> &queried_tables, vector<string> &where_conditions) {
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
			where_conditions.push_back("t" + to_string(l_bind.table_index) + "." + bta[l_bind.ToString()] + " IN (" +
			                           JoinStrings(in_vals, ",") + ")");
		} else {
			for (idx_t cond_i = 0; cond_i < join->conditions.size(); cond_i++) {
				JoinCondition &join_condition = join->conditions[cond_i];
				auto l_bind = ((BoundColumnRefExpression &)*join_condition.left.get()).binding;
				auto r_bind = ((BoundColumnRefExpression &)*join_condition.right.get()).binding;
				where_conditions.push_back("t" + to_string(l_bind.table_index) + "." + bta[l_bind.ToString()] + " = " +
				                           "t" + to_string(r_bind.table_index) + "." + bta[r_bind.ToString()]);
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
		CreateSubQuery(*child, temporary_table_name, queried_tables, where_conditions);
	}
	// grab selected columns
	vector<ColumnBinding> selected_column_bindings;
	if (plan.type == LogicalOperatorType::COMPARISON_JOIN) {
		LogicalComparisonJoin *join = static_cast<LogicalComparisonJoin *>(&plan);
		selected_column_bindings = GetColumnBindings(*join);
	} else {
		selected_column_bindings = plan.GetColumnBindings();
	}
	vector<string> selected_columns;
	for (auto cb : selected_column_bindings)
		selected_columns.push_back("t" + to_string(cb.table_index) + "." + bta[cb.ToString()] + " AS " + bta[cb.ToString()] + to_string(cb.table_index));
	// create query and return
	return "CREATE TEMPORARY TABLE main." + temporary_table_name + " AS (" + "SELECT " +
	       JoinStrings(selected_columns, ", ") + " " + "FROM " + JoinStrings(queried_tables, ", ") + " " + "WHERE " +
	       JoinStrings(where_conditions, " AND ") + ");";
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
			return "t" + to_string(binding.table_index) + "." + bta[binding.ToString()] + " IS NOT NULL";
		}
		case ExpressionType::OPERATOR_IS_NULL: {
			auto *op = static_cast<BoundOperatorExpression *>(expr);
			auto binding = static_cast<BoundColumnRefExpression *>(op->children[0].get())->binding;
			return "t" + to_string(binding.table_index) + "." + bta[binding.ToString()] + " IS NULL";
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
		string condition = "t" + to_string(binding.table_index) + "." + bta[binding.ToString()] + " BETWEEN ";
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
	string condition = "t" + to_string(binding.table_index) + "." + bta[binding.ToString()] + " " +
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
	string condition = "t" + to_string(binding.table_index) + "." + bta[binding.ToString()];
	if (func->function.name[0] == '!')
		condition += " NOT";
	if (func->function.name.find("~~") != string::npos)
		condition += " LIKE ";
	if (func->children[0]->return_type == TypeId::VARCHAR) {
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

void ReOptimizer::InjectCardinalities(LogicalOperator &plan, string temp_table_name) {
	// TODO: derive information about the cardinality of the intermediate results, implement
	// TODO: create cost function that uses injected cardinalities
	cardinalities[JoinStrings(GetRelationSet(plan), ",")] = GetTable("main", temp_table_name)->storage->cardinality;
}

vector<vector<string>> ReOptimizer::TempTablePowerset() {
	idx_t powerset_size = pow(2, temp_tables.size());
	vector<vector<string>> powerset(powerset_size);
	for (idx_t i = 0; i < powerset_size; i++) {
		vector<string> subset;
		for (idx_t j = 0; j < temp_tables.size(); j++) {
			if (i & (1 << j)) {
				subset.push_back(temp_tables[j]);
			}
		}
		powerset.push_back(subset);
	}
	return powerset;
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

void ReOptimizer::ReplaceLogicalOperator(LogicalOperator &plan, LogicalOperator &old_op, TableCatalogEntry *table,
                                         idx_t depth) {
	if (plan.children.empty())
		return;
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
		return;
	}
	// enter recursion until the operator is found
	for (auto &child : plan.children)
		ReplaceLogicalOperator(*child.get(), old_op, table, depth + 1);
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

// FIXME: disabling optimizer indeed causes less overhead (ONLY TESTED FOR SMALL QUERIES), but it's not very clean
// we already have an optimized subquery plan at this point
// slapping a create table and project on top of it should be doable... right?
// also don't know if not having filter pushdown sucks - join order should be fixed tho (by how query is made)

// FIXME: getting cardinality from data_storage only works on autocommit, assume no transactions for now
// this might be because of the fact that cardinality is atomic?
void ReOptimizer::ExecuteSubQuery(const string subquery) {
	// store state of autocommit, profiler and optimizer
	// bool auto_commit = context.transaction.IsAutoCommit();
	context.transaction.Commit();
	bool profiler = context.profiler.IsEnabled();
	// bool optimizer = context.enable_optimizer;
	// disable them
	// context.transaction.SetAutoCommit(false);
	context.profiler.Disable();
	// context.enable_optimizer = false;
	// hack in a query - without autocommit, profiling, optimizer or lock
	// this prevents segfault / assertion fail / unnecessary overhead / deadlock respectively
	context.QueryWithoutLock(subquery, false);
	// restore the state of autocommit, profiler and optimizer
	// context.transaction.SetAutoCommit(auto_commit);
	context.transaction.BeginTransaction();
	// if (profiler)
	context.profiler.Enable();
	// context.enable_optimizer = optimizer;
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
