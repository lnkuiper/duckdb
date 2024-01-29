//===----------------------------------------------------------------------===//
//
//                         DuckDB
//
// dbgen.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb.hpp"
#ifndef DUCKDB_AMALGAMATION
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/common/string.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#endif

namespace duckdb {
class ClientContext;
}

namespace tpch {

struct DBGenWrapper {
	//! Create the TPC-H tables in the given catalog / schema with the given suffix
	static void CreateTPCHSchema(duckdb::ClientContext &context, duckdb::string catalog, duckdb::string schema,
	                             duckdb::string suffix);
	//! Load the TPC-H data of the given scale factor
	static void LoadTPCHData(duckdb::ClientContext &context, double sf, duckdb::string catalog, duckdb::string schema,
	                         duckdb::string suffix, int children, int step);

	//! Gets the specified TPC-H Query number as a duckdb::string
	static duckdb::string GetQuery(int query);
	//! Returns the CSV answer of a TPC-H query
	static duckdb::string GetAnswer(double sf, int query);
};

} // namespace tpch
