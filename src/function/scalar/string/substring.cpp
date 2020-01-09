#include "duckdb/function/scalar/string_functions.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"

using namespace std;

namespace duckdb {

static const char *substring_scalar_function(const char *input_string, int offset, int length,
                                             unique_ptr<char[]> &output, index_t &current_len) {
	// reduce offset by one because SQL starts counting at 1
	offset--;

	if (offset < 0 || length < 0) {
		throw Exception("SUBSTRING cannot handle negative offsets");
	}

	index_t required_len = strlen(input_string) + 1;
	if (required_len > current_len) {
		// need a resize
		current_len = required_len;
		output = unique_ptr<char[]>{new char[required_len]};
	}

	// UTF8 chars can use more than one byte
	index_t input_char_offset = 0;
	index_t input_byte_offset = 0;
	index_t output_byte_offset = 0;

	while (input_string[input_byte_offset]) {
		char b = input_string[input_byte_offset++];
		input_char_offset += (b & 0xC0) != 0x80;
		if (input_char_offset > (index_t)(offset + length)) {
			break;
		}
		if (input_char_offset > (index_t)offset) {
			output[output_byte_offset++] = b;
		}
	}
	// terminate output
	output[output_byte_offset] = '\0';
	return output.get();
}

static void substring_function(DataChunk &args, ExpressionState &state, Vector &result) {
	assert(args.column_count == 3 && args.data[0].type == TypeId::VARCHAR && args.data[1].type == TypeId::INTEGER &&
	       args.data[2].type == TypeId::INTEGER);
	auto &input_vector = args.data[0];
	auto &offset_vector = args.data[1];
	auto &length_vector = args.data[2];

	index_t current_len = 0;
	unique_ptr<char[]> output;
	VectorOperations::TernaryExec<const char *, int, int, const char *>(
	    input_vector, offset_vector, length_vector, result,
	    [&](const char *input_string, int offset, int length, index_t result_index) {
		    return result.string_heap.AddString(
		        substring_scalar_function(input_string, offset, length, output, current_len));
	    });
}

void SubstringFun::RegisterFunction(BuiltinFunctions &set) {
	set.AddFunction(ScalarFunction("substring",       // name of function
	                               {SQLType::VARCHAR, // argument list
	                                SQLType::INTEGER, SQLType::INTEGER},
	                               SQLType::VARCHAR,     // return type
	                               substring_function)); // pointer to function implementation
}

} // namespace duckdb
