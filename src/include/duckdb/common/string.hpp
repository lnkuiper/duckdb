//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/string.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include <string>

namespace duckdb {

// FIXME: we would like to change this to container_allocator<_CharT> instead, so it uses our allocator.
//  However, this is a big pain as it will no longer be compatible with std::string, which we use absolutely everywhere.
//  So, for now, we use std::allocator<_CharT>, and our string is exactly the same as std::string.
template <class _CharT = char>
using string_allocator = std::allocator<_CharT>;

template <class _CharT = char, class _Traits = std::char_traits<_CharT>, class _Alloc = string_allocator<_CharT>>
using basic_string = typename std::basic_string<_CharT, _Traits, _Alloc>;
using string = basic_string<char>;

// If we were to use container_allocator instead, we would have to re-define these for our custom string type.
using std::stod;
using std::stof;
using std::stoi;
using std::stol;
using std::stold;
using std::stoll;
using std::stoul;
using std::stoull;

} // namespace duckdb
