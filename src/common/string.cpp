#include "duckdb/common/string.hpp"

#include "duckdb/common/types/hash.hpp"

namespace duckdb {

int stoi(const string &__str) {
	return std::stoi(__str.c_str());
}
long stol(const string &__str) {
	return std::stol(__str.c_str());
}
unsigned long stoul(const string &__str) {
	return std::stoul(__str.c_str());
}
long long stoll(const string &__str) {
	return std::stoll(__str.c_str());
}
unsigned long long stoull(const string &__str) {
	return std::stoull(__str.c_str());
}
float stof(const string &__str) {
	return std::stof(__str.c_str());
}
double stod(const string &__str) {
	return std::stod(__str.c_str());
}
long double stold(const string &__str) {
	return std::stold(__str.c_str());
}

size_t DuckStringHash(const duckdb::string &val) {
	return Hash(val.c_str(), val.length());
}

} // namespace duckdb
