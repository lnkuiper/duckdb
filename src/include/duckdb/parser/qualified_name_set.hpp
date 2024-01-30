//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/parser/qualified_name_set.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/types/hash.hpp"
#include "duckdb/common/unordered_set.hpp"
#include "duckdb/parser/qualified_name.hpp"

namespace duckdb {

struct QualifiedColumnHashFunction {
	uint64_t operator()(const QualifiedColumnName &a) const {
		return CombineHash(CombineHash(Hash(a.schema), Hash(a.table)), Hash(a.column));
	}
};

struct QualifiedColumnEquality {
	bool operator()(const QualifiedColumnName &a, const QualifiedColumnName &b) const {
		return a.schema == b.schema && a.table == b.table && a.column == b.column;
	}
};

using qualified_column_set_t = unordered_set<QualifiedColumnName, QualifiedColumnHashFunction, QualifiedColumnEquality>;

} // namespace duckdb
