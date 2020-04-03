#include "duckdb/planner/operator/logical_comparison_join.hpp"
#include "duckdb/common/string_util.hpp"

#include "duckdb/common/printer.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"

using namespace duckdb;
using namespace std;

LogicalComparisonJoin::LogicalComparisonJoin(JoinType join_type, LogicalOperatorType logical_type)
    : LogicalJoin(join_type, logical_type) {
}

vector<ColumnBinding> LogicalComparisonJoin::GetColumnBindings() {
	auto left_bindings = MapBindings(children[0]->GetColumnBindings(), left_projection_map);
	if (join_type == JoinType::SEMI || join_type == JoinType::ANTI) {
		// for SEMI and ANTI join we only project the left hand side
		return left_bindings;
	}
	if (join_type == JoinType::MARK) {
		// for MARK join we project the left hand side plus the MARK column
		left_bindings.push_back(ColumnBinding(mark_index, 0));
		return left_bindings;
	}
	// string testing = "COMPARISON JOIN " + ParamsToString();
	// for (JoinCondition &jc : conditions) {
	// 	BoundColumnRefExpression &l = (BoundColumnRefExpression &)*jc.left.get();
	// 	BoundColumnRefExpression &r = (BoundColumnRefExpression &)*jc.right.get();
	// 	testing += l.ToString();
	// 	testing += r.ToString();
	// 	testing += " ";
	// }

	// testing += "Bindings: ";
	// for (ColumnBinding cb : children[0]->GetColumnBindings()) {
	// 	testing += cb.ToString();
	// }
	// testing += " - ";
	// for (ColumnBinding cb : children[1]->GetColumnBindings()) {
	// 	testing += cb.ToString();
	// }
	// testing += ", l:";
	// for (idx_t i : left_projection_map) {
	// 	testing += to_string(i);
	// }
	// testing += ", r:";
	// for (idx_t i : right_projection_map) {
	// 	testing += to_string(i);
	// }

	// for other join types we project both the LHS and the RHS
	auto right_bindings = MapBindings(children[1]->GetColumnBindings(), right_projection_map);
	left_bindings.insert(left_bindings.end(), right_bindings.begin(), right_bindings.end());

	// testing += " || ";
	// for (ColumnBinding cb : left_bindings) {
	// 	testing += cb.ToString();
	// }
	// Printer::Print(testing);

	return left_bindings;
}

string LogicalComparisonJoin::ParamsToString() const {
	string result = "[" + JoinTypeToString(join_type);
	if (conditions.size() > 0) {
		result += " ";
		result += StringUtil::Join(conditions, conditions.size(), ", ", [](const JoinCondition &condition) {
			return ExpressionTypeToString(condition.comparison) + "(" + condition.left->GetName() + ", " +
			       condition.right->GetName() + ")";
		});
		result += "]";
	}

	return result;
}
