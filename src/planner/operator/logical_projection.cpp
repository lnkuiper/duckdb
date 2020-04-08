#include "duckdb/planner/operator/logical_projection.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"

using namespace duckdb;
using namespace std;

LogicalProjection::LogicalProjection(idx_t table_index, vector<unique_ptr<Expression>> select_list)
    : LogicalOperator(LogicalOperatorType::PROJECTION, move(select_list)), table_index(table_index) {
}

vector<ColumnBinding> LogicalProjection::GetColumnBindings() {
	vector<ColumnBinding> cbs = GenerateColumnBindings(table_index, expressions.size());
	return GenerateColumnBindings(table_index, expressions.size());
}

void LogicalProjection::ResolveTypes() {
	for (auto &expr : expressions) {
		types.push_back(expr->return_type);
	}
}

string LogicalProjection::ParamsToString() const {
	string result = "";
	if (expressions.size() > 0) {
		result += "[";
		result += StringUtil::Join(expressions, expressions.size(), ", ", [](const unique_ptr<Expression> &expression) {
			return expression->GetName() + ((BoundColumnRefExpression &)*expression.get()).ToString();
		});
		result += "]";
	}

	return result;
}
