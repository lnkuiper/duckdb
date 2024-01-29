//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/main/database_path_and_type.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/string.hpp"
#include "duckdb/main/config.hpp"

namespace duckdb {

struct DBPathAndType {
	//! Parse database extension type and rest of path from combined form (type:path)
	static void ExtractExtensionPrefix(string &path, string &db_type);
	//! Check the magic bytes of a file and set the database type based on that
	static void CheckMagicBytes(string &path, string &db_type, const DBConfig &config);

	//! Run ExtractExtensionPrefix followed by CheckMagicBytes
	static void ResolveDatabaseType(string &path, string &db_type, const DBConfig &config);
};
} // namespace duckdb
