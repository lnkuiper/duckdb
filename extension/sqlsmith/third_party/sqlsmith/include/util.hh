#ifndef UTIL_HH
#define UTIL_HH
#include <typeinfo>

#include "duckdb/common/stringstream.hpp"

using namespace std;

/* TODO: The strings are implementation-defined.  How do they look in
   clang? */

inline duckdb::string pretty_type(const char *raw) {
	duckdb::ostringstream os;
	os << raw;
	duckdb::string s = os.str();
	while (s[0] <= '9')
		s.erase(s.begin());
	return s;
}

inline duckdb::string pretty_type(struct prod *p) {
	return pretty_type(typeid(*p).name());
}

#endif
