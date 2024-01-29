// Copyright 2016 The RE2 Authors.  All Rights Reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef UTIL_STRUTIL_H_
#define UTIL_STRUTIL_H_

#include "duckdb/common/string.hpp"

#include "re2/stringpiece.h"
#include "util/util.h"

namespace duckdb_re2 {

duckdb::string CEscape(const StringPiece& src);
void PrefixSuccessor(duckdb::string* prefix);
duckdb::string StringPrintf(const char* format, ...);
void SStringPrintf(duckdb::string* dst, const char* format, ...);
void StringAppendF(duckdb::string* dst, const char* format, ...);

}  // namespace duckdb_re2

#endif  // UTIL_STRUTIL_H_
