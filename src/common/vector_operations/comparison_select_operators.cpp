//===--------------------------------------------------------------------===//
// comparison_select_operators.cpp
// Description: This file contains the implementation of the comparison select
// operations == != >= <= > <
//===--------------------------------------------------------------------===//

#include "duckdb/common/operator/comparison_operators.hpp"

#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/common/vector_operations/binary_executor.hpp"

using namespace duckdb;
using namespace std;

template <class OP> static index_t templated_select_operation(Vector &left, Vector &right, sel_t result[]) {
	// the inplace loops take the result as the last parameter
	switch (left.type) {
	case TypeId::BOOL:
	case TypeId::INT8:
		return BinaryExecutor::Select<int8_t, int8_t, OP>(left, right, result);
	case TypeId::INT16:
		return BinaryExecutor::Select<int16_t, int16_t, OP>(left, right, result);
	case TypeId::INT32:
		return BinaryExecutor::Select<int32_t, int32_t, OP>(left, right, result);
	case TypeId::INT64:
		return BinaryExecutor::Select<int64_t, int64_t, OP>(left, right, result);
	case TypeId::POINTER:
		return BinaryExecutor::Select<uint64_t, uint64_t, OP>(left, right, result);
	case TypeId::FLOAT:
		return BinaryExecutor::Select<float, float, OP>(left, right, result);
	case TypeId::DOUBLE:
		return BinaryExecutor::Select<double, double, OP>(left, right, result);
	case TypeId::VARCHAR:
		return BinaryExecutor::Select<const char *, const char *, OP>(left, right, result);
	default:
		throw InvalidTypeException(left.type, "Invalid type for comparison");
	}
}

index_t VectorOperations::SelectEquals(Vector &left, Vector &right, sel_t result[]) {
	return templated_select_operation<duckdb::Equals>(left, right, result);
}

index_t VectorOperations::SelectNotEquals(Vector &left, Vector &right, sel_t result[]) {
	return templated_select_operation<duckdb::NotEquals>(left, right, result);
}

index_t VectorOperations::SelectGreaterThanEquals(Vector &left, Vector &right, sel_t result[]) {
	return templated_select_operation<duckdb::GreaterThanEquals>(left, right, result);
}

index_t VectorOperations::SelectLessThanEquals(Vector &left, Vector &right, sel_t result[]) {
	return templated_select_operation<duckdb::LessThanEquals>(left, right, result);
}

index_t VectorOperations::SelectGreaterThan(Vector &left, Vector &right, sel_t result[]) {
	return templated_select_operation<duckdb::GreaterThan>(left, right, result);
}

index_t VectorOperations::SelectLessThan(Vector &left, Vector &right, sel_t result[]) {
	return templated_select_operation<duckdb::LessThan>(left, right, result);
}
