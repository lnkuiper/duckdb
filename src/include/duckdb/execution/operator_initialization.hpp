//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/execution/operator_initialization.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/enums/state_initialization.hpp"
#include "duckdb/common/optional_ptr.hpp"
#include "duckdb/planner/expression.hpp"

namespace duckdb {

struct ExpressionExecutorGlobalState;

struct OperatorInitializationInfo {
	StateInitialization effect = StateInitialization::UNKNOWN;
	optional_ptr<const vector<unique_ptr<Expression>>> expressions;

	static OperatorInitializationInfo Independent() {
		OperatorInitializationInfo result;
		result.effect = StateInitialization::INDEPENDENT;
		return result;
	}
	static OperatorInitializationInfo Shared(const vector<unique_ptr<Expression>> &expressions) {
		OperatorInitializationInfo result;
		result.effect = StateInitialization::SHARED;
		result.expressions = expressions;
		return result;
	}
};

struct OperatorStateInitInput {
	//! Query-owned states prepared using this operator's negotiated input guarantees.
	optional_ptr<const ExpressionExecutorGlobalState> expressions;
};

} // namespace duckdb
