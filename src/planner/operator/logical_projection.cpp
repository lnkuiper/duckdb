#include "duckdb/planner/operator/logical_projection.hpp"

#include "duckdb/common/printer.hpp"

using namespace duckdb;
using namespace std;

LogicalProjection::LogicalProjection(idx_t table_index, vector<unique_ptr<Expression>> select_list)
    : LogicalOperator(LogicalOperatorType::PROJECTION, move(select_list)), table_index(table_index) {
}

vector<ColumnBinding> LogicalProjection::GetColumnBindings() {
	vector<ColumnBinding> cbs = GenerateColumnBindings(table_index, expressions.size());
	// string testing = "PROJECTION " + ParamsToString() + " Bindings: ";
	// for (ColumnBinding cb : cbs) {
	// 	testing += cb.ToString();
	// }
	// testing += " || ";
	// for (auto &expr : expressions) {
	// 	testing += expr->alias + "-" + expr->ToString();
	// }
	// Printer::Print(testing);
	return GenerateColumnBindings(table_index, expressions.size());
}

void LogicalProjection::ResolveTypes() {
	for (auto &expr : expressions) {
		types.push_back(expr->return_type);
	}
}
