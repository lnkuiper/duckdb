//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/string.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/container_allocator.hpp"

#include <string>

namespace duckdb {

template <class _CharT = char, class _Traits = std::char_traits<_CharT>, class _Alloc = container_allocator<_CharT>>
using basic_string = typename std::basic_string<_CharT, _Traits, _Alloc>;
using string = basic_string<char>;

int stoi(const string &__str);
long stol(const string &__str);
unsigned long stoul(const string &__str);
long long stoll(const string &__str);
unsigned long long stoull(const string &__str);
float stof(const string &__str);
double stod(const string &__str);
long double stold(const string &__str);

size_t DuckStringHash(const duckdb::string &val);

} // namespace duckdb

namespace std {

template <>
struct hash<duckdb::string> {
	size_t operator()(const duckdb::string &__val) const {
		return duckdb::DuckStringHash(__val);
	}
};

} // namespace std
