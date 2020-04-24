#include "duckdb/function/scalar/string_functions.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "utf8proc.hpp"

using namespace std;

namespace duckdb {

// length returns the size in characters
struct StringLengthOperator {
	template <class TA, class TR> static inline TR Operation(TA input) {
		auto input_data = input.GetData();
		auto input_length = input.GetSize();
		for (idx_t i = 0; i < input_length; i++) {
			if (input_data[i] & 0x80) {
				int64_t length = 0;
				// non-ascii character: use grapheme iterator on remainder of string
				utf8proc_grapheme_callback(input_data, input_length, [&](size_t start, size_t end) {
					length++;
					return true;
				});
				return length;
			}
		}
		return input_length;
	}
};

// strlen returns the size in bytes
struct StrLenOperator {
	template <class TA, class TR> static inline TR Operation(TA input) {
		return input.GetSize();
	}
};

void LengthFun::RegisterFunction(BuiltinFunctions &set) {
	set.AddFunction({"length", "len"}, ScalarFunction({SQLType::VARCHAR}, SQLType::BIGINT,
	                               ScalarFunction::UnaryFunction<string_t, int64_t, StringLengthOperator, true>));
	set.AddFunction(ScalarFunction("strlen", {SQLType::VARCHAR}, SQLType::BIGINT,
	                               ScalarFunction::UnaryFunction<string_t, int64_t, StrLenOperator, true>));
}

} // namespace duckdb
