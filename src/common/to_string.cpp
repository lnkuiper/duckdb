#include "duckdb/common/to_string.hpp"

#include "duckdb/common/stringstream.hpp"

namespace duckdb {

template <typename T>
string to_string_internal(const T &t) {
	ostringstream o;
	o << t;
	return o.str();
}

string to_string(int __val) {
	return to_string_internal(__val);
}
string to_string(unsigned __val) {
	return to_string_internal(__val);
}
string to_string(long __val) {
	return to_string_internal(__val);
}
string to_string(unsigned long __val) {
	return to_string_internal(__val);
}
string to_string(long long __val) {
	return to_string_internal(__val);
}
string to_string(unsigned long long __val) {
	return to_string_internal(__val);
}
string to_string(float __val) {
	return to_string_internal(__val);
}
string to_string(double __val) {
	return to_string_internal(__val);
}
string to_string(long double __val) {
	return to_string_internal(__val);
}

} // namespace duckdb
