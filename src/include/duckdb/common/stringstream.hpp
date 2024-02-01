//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/stringstream.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/string.hpp"

#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>

namespace duckdb {

template <class _CharT = char, class _Traits = std::char_traits<_CharT>, class _Alloc = string_allocator<_CharT>>
using basic_stringstream = typename std::basic_stringstream<_CharT, _Traits, _Alloc>;
using stringstream = basic_stringstream<char>;

template <class _CharT = char, class _Traits = std::char_traits<_CharT>, class _Alloc = string_allocator<_CharT>>
using basic_istringstream = typename std::basic_istringstream<_CharT, _Traits, _Alloc>;
using istringstream = basic_istringstream<char>;

template <class _CharT = char, class _Traits = std::char_traits<_CharT>, class _Alloc = string_allocator<_CharT>>
using basic_ostringstream = typename std::basic_ostringstream<_CharT, _Traits, _Alloc>;
using ostringstream = basic_ostringstream<char>;

using std::istream;
using std::ostream;

} // namespace duckdb
