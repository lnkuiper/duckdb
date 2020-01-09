#include "duckdb/function/scalar/operators.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"

using namespace duckdb;
using namespace std;

namespace duckdb {

//===--------------------------------------------------------------------===//
// & [bitwise_and]
//===--------------------------------------------------------------------===//
struct BitwiseANDOperator {
	template <class T> static inline T Operation(T left, T right) {
		return left & right;
	}
};

void BitwiseAndFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions("&");
	for (auto &type : SQLType::INTEGRAL) {
		functions.AddFunction(ScalarFunction({type, type}, type,
		                                     ScalarFunction::GetScalarIntegerBinaryFunction<BitwiseANDOperator>(type)));
	}
	set.AddFunction(functions);
}

//===--------------------------------------------------------------------===//
// | [bitwise_or]
//===--------------------------------------------------------------------===//
struct BitwiseOROperator {
	template <class T> static inline T Operation(T left, T right) {
		return left | right;
	}
};

void BitwiseOrFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions("|");
	for (auto &type : SQLType::INTEGRAL) {
		functions.AddFunction(ScalarFunction({type, type}, type,
		                                     ScalarFunction::GetScalarIntegerBinaryFunction<BitwiseOROperator>(type)));
	}
	set.AddFunction(functions);
}

//===--------------------------------------------------------------------===//
// # [bitwise_xor]
//===--------------------------------------------------------------------===//
struct BitwiseXOROperator {
	template <class T> static inline T Operation(T left, T right) {
		return left ^ right;
	}
};

void BitwiseXorFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions("#");
	for (auto &type : SQLType::INTEGRAL) {
		functions.AddFunction(ScalarFunction({type, type}, type,
		                                     ScalarFunction::GetScalarIntegerBinaryFunction<BitwiseXOROperator>(type)));
	}
	set.AddFunction(functions);
}

//===--------------------------------------------------------------------===//
// << [bitwise_left_shift]
//===--------------------------------------------------------------------===//
struct BitwiseShiftLeftOperator {
	template <class T> static inline T Operation(T left, T right) {
		return left << right;
	}
};

void LeftShiftFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions("<<");
	for (auto &type : SQLType::INTEGRAL) {
		functions.AddFunction(ScalarFunction(
		    {type, type}, type, ScalarFunction::GetScalarIntegerBinaryFunction<BitwiseShiftLeftOperator>(type)));
	}
	set.AddFunction(functions);
}

//===--------------------------------------------------------------------===//
// >> [bitwise_right_shift]
//===--------------------------------------------------------------------===//
struct BitwiseShiftRightOperator {
	template <class T> static inline T Operation(T left, T right) {
		return left >> right;
	}
};

void RightShiftFun::RegisterFunction(BuiltinFunctions &set) {
	ScalarFunctionSet functions(">>");
	for (auto &type : SQLType::INTEGRAL) {
		functions.AddFunction(ScalarFunction(
		    {type, type}, type, ScalarFunction::GetScalarIntegerBinaryFunction<BitwiseShiftRightOperator>(type)));
	}
	set.AddFunction(functions);
}

} // namespace duckdb
