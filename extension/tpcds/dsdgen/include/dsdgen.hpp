//===----------------------------------------------------------------------===//
//
//                         DuckDB
//
// dsdgen.hpp
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

namespace tpcds {

struct DSDGenWrapper {
	//! Create the TPC-DS tables in the given schema with the given suffix
	static void CreateTPCDSSchema(duckdb::ClientContext &context, duckdb::string catalog, duckdb::string schema,
	                              duckdb::string suffix, bool keys, bool overwrite);
	//! Generate the TPC-DS data of the given scale factor
	static void DSDGen(double scale, duckdb::ClientContext &context, duckdb::string catalog, duckdb::string schema,
	                   duckdb::string suffix);

	static uint32_t QueriesCount();
	//! Gets the specified TPC-DS Query number as a duckdb::string
	static duckdb::string GetQuery(int query);
	//! Returns the CSV answer of a TPC-DS query
	static duckdb::string GetAnswer(double sf, int query);
};

} // namespace tpcds
