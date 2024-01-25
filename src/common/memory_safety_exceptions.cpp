#include "duckdb/common/memory_safety_exceptions.hpp"

#include "duckdb/common/exception.hpp"

namespace duckdb {

std::exception VectorOutOfBoundsException(idx_t index, idx_t size) {
	return InternalException("Attempted to access index %ld within vector of size %ld", index, size);
}

std::exception VectorBackOnEmptyException() {
	return InternalException("'back' called on an empty vector!");
}

std::exception UniquePtrNullException() {
	return InternalException("Attempted to dereference unique_ptr that is NULL!");
}

std::exception OptionalPtrNullException() {
	return InternalException("Attempting to dereference an optional pointer that is not set");
}

} // namespace duckdb
