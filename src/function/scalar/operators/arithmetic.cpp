#include "duckdb/function/scalar/operators.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/common/operator/numeric_binary_operators.hpp"

using namespace std;

namespace duckdb {

//===--------------------------------------------------------------------===//
// + [add]
//===--------------------------------------------------------------------===//
void AddFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions("+");
	// binary add function adds two numbers together
	for (auto &type : SQLType::NUMERIC) {
		functions.AddFunction(
		    ScalarFunction({type, type}, type, ScalarFunction::GetScalarBinaryFunction<AddOperator>(type)));
	}
	// we can add integers to dates
	functions.AddFunction(ScalarFunction({SQLType::DATE, SQLType::INTEGER}, SQLType::DATE,
	                                     ScalarFunction::GetScalarBinaryFunction<AddOperator>(SQLType::INTEGER)));
	functions.AddFunction(ScalarFunction({SQLType::INTEGER, SQLType::DATE}, SQLType::DATE,
	                                     ScalarFunction::GetScalarBinaryFunction<AddOperator>(SQLType::INTEGER)));
	// unary add function is a nop, but only exists for numeric types
	for (auto &type : SQLType::NUMERIC) {
		functions.AddFunction(ScalarFunction({type}, type, ScalarFunction::NopFunction));
	}
	set.AddFunction(functions);
}

//===--------------------------------------------------------------------===//
// - [subtract]
//===--------------------------------------------------------------------===//
void SubtractFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions("-");
	// binary subtract function "a - b", subtracts b from a
	for (auto &type : SQLType::NUMERIC) {
		functions.AddFunction(
		    ScalarFunction({type, type}, type, ScalarFunction::GetScalarBinaryFunction<SubtractOperator>(type)));
	}
	functions.AddFunction(ScalarFunction({SQLType::DATE, SQLType::DATE}, SQLType::INTEGER,
	                                     ScalarFunction::GetScalarBinaryFunction<SubtractOperator>(SQLType::INTEGER)));
	functions.AddFunction(ScalarFunction({SQLType::DATE, SQLType::INTEGER}, SQLType::DATE,
	                                     ScalarFunction::GetScalarBinaryFunction<SubtractOperator>(SQLType::INTEGER)));
	// unary subtract function, negates the input (i.e. multiplies by -1)
	for (auto &type : SQLType::NUMERIC) {
		functions.AddFunction(
		    ScalarFunction({type}, type, ScalarFunction::GetScalarUnaryFunction<NegateOperator>(type)));
	}
	set.AddFunction(functions);
}

//===--------------------------------------------------------------------===//
// * [multiply]
//===--------------------------------------------------------------------===//
void MultiplyFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions("*");
	for (auto &type : SQLType::NUMERIC) {
		functions.AddFunction(
		    ScalarFunction({type, type}, type, ScalarFunction::GetScalarBinaryFunction<MultiplyOperator>(type)));
	}
	set.AddFunction(functions);
}

//===--------------------------------------------------------------------===//
// / [divide]
//===--------------------------------------------------------------------===//
struct BinaryZeroIsNullWrapper {
	template <class FUNC, class OP, class LEFT_TYPE, class RIGHT_TYPE, class RESULT_TYPE>
	static inline RESULT_TYPE Operation(FUNC fun, LEFT_TYPE left, RIGHT_TYPE right, nullmask_t &nullmask, idx_t idx) {
		if (right == 0) {
			nullmask[idx] = true;
			return 0;
		} else {
			return OP::template Operation<LEFT_TYPE, RIGHT_TYPE, RESULT_TYPE>(left, right);
		}
	}
};

template <class T, class OP>
static void BinaryScalarFunctionIgnoreZero(DataChunk &input, ExpressionState &state, Vector &result) {
	BinaryExecutor::Execute<T, T, T, OP, true, BinaryZeroIsNullWrapper>(input.data[0], input.data[1], result,
	                                                                    input.size());
}

template <class OP> static scalar_function_t GetBinaryFunctionIgnoreZero(SQLType type) {
	switch (type.id) {
	case SQLTypeId::TINYINT:
		return BinaryScalarFunctionIgnoreZero<int8_t, OP>;
	case SQLTypeId::SMALLINT:
		return BinaryScalarFunctionIgnoreZero<int16_t, OP>;
	case SQLTypeId::INTEGER:
		return BinaryScalarFunctionIgnoreZero<int32_t, OP>;
	case SQLTypeId::BIGINT:
		return BinaryScalarFunctionIgnoreZero<int64_t, OP>;
	case SQLTypeId::FLOAT:
		return BinaryScalarFunctionIgnoreZero<float, OP>;
	case SQLTypeId::DOUBLE:
		return BinaryScalarFunctionIgnoreZero<double, OP>;
	case SQLTypeId::DECIMAL:
		return BinaryScalarFunctionIgnoreZero<double, OP>;
	default:
		throw NotImplementedException("Unimplemented type for GetScalarUnaryFunction");
	}
}

void DivideFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions("/");
	for (auto &type : SQLType::NUMERIC) {
		functions.AddFunction(ScalarFunction({type, type}, type, GetBinaryFunctionIgnoreZero<DivideOperator>(type)));
	}
	set.AddFunction(functions);
}

//===--------------------------------------------------------------------===//
// % [modulo]
//===--------------------------------------------------------------------===//
template <> float ModuloOperator::Operation(float left, float right) {
	assert(right != 0);
	return fmod(left, right);
}

template <> double ModuloOperator::Operation(double left, double right) {
	assert(right != 0);
	return fmod(left, right);
}

void ModFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions("%");
	for (auto &type : SQLType::NUMERIC) {
		functions.AddFunction(ScalarFunction({type, type}, type, GetBinaryFunctionIgnoreZero<ModuloOperator>(type)));
	}
	set.AddFunction(functions);
	functions.name = "mod";
	set.AddFunction(functions);
}

} // namespace duckdb
