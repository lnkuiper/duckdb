//===----------------------------------------------------------------------===//
//                         DuckDB
//
// postgres_parser.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/string.hpp"
#include "duckdb/common/vector.hpp"
#include "nodes/pg_list.hpp"
#include "pg_simplified_token.hpp"

namespace duckdb {
class PostgresParser {
public:
	PostgresParser();
	~PostgresParser();

	bool success;
	duckdb_libpgquery::PGList *parse_tree;
	string error_message;
	int error_location;

public:
	void Parse(const string &query);
	static vector<duckdb_libpgquery::PGSimplifiedToken> Tokenize(const string &query);

	static bool IsKeyword(const string &text);
	static vector<duckdb_libpgquery::PGKeyword> KeywordList();

	static void SetPreserveIdentifierCase(bool downcase);
};

} // namespace duckdb
