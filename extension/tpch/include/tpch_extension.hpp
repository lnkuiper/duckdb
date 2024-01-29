//===----------------------------------------------------------------------===//
//                         DuckDB
//
// tpch_extension.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb.hpp"

namespace duckdb {

class TpchExtension : public Extension {
public:
	void Load(DuckDB &db) override;
	string Name() override;

	//! Gets the specified TPC-H Query number as a string
	static string GetQuery(int query);
	//! Returns the CSV answer of a TPC-H query
	static string GetAnswer(double sf, int query);
};

} // namespace duckdb
