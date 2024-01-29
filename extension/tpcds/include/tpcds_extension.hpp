//===----------------------------------------------------------------------===//
//                         DuckDB
//
// tpcds_extension.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb.hpp"
#include "duckdb/main/client_context.hpp"

namespace duckdb {

class TpcdsExtension : public Extension {
public:
	void Load(DuckDB &db) override;
	string Name() override;

	//! Gets the specified TPC-DS Query number as a string
	static string GetQuery(int query);
	//! Returns the CSV answer of a TPC-DS query
	static string GetAnswer(double sf, int query);
};

} // namespace duckdb
