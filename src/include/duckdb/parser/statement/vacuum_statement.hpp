//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/parser/statement/explain_statement.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/parser/parsed_expression.hpp"
#include "duckdb/parser/sql_statement.hpp"
#include "duckdb/parser/parsed_data/vacuum_info.hpp"

namespace duckdb {

class VacuumStatement : public SQLStatement {
public:
	VacuumStatement() : SQLStatement(StatementType::VACUUM){};
	unique_ptr<VacuumInfo> info;

};

} // namespace duckdb
