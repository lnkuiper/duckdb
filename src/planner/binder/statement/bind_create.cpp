#include "duckdb/parser/statement/create_statement.hpp"
#include "duckdb/planner/statement/bound_create_statement.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/schema_catalog_entry.hpp"
#include "duckdb/planner/binder.hpp"

using namespace duckdb;
using namespace std;

unique_ptr<BoundSQLStatement> Binder::Bind(CreateStatement &stmt) {
	auto result = make_unique<BoundCreateStatement>();
	result->info = BindCreateInfo(move(stmt.info));
	return move(result);
}

unique_ptr<BoundCreateInfo> Binder::BindCreateInfo(unique_ptr<CreateInfo> info) {
	unique_ptr<BoundCreateInfo> result;
	if (info->schema == INVALID_SCHEMA) {
		info->schema = info->temporary ? TEMP_SCHEMA : DEFAULT_SCHEMA;
	}

	SchemaCatalogEntry *bound_schema = nullptr;
	if (!info->temporary) {
		// non-temporary create: not read only
		if (info->schema == TEMP_SCHEMA) {
			throw ParserException("Only TEMPORARY table names can use the \"temp\" schema");
		}
		this->read_only = false;
	} else {
		if (info->schema != TEMP_SCHEMA) {
			throw ParserException("TEMPORARY table names can *only* use the \"%s\" schema", TEMP_SCHEMA);
		}
	}
	if (info->type != CatalogType::SCHEMA) {
		// fetch the schema in which we want to create the object
		bound_schema = Catalog::GetCatalog(context).GetSchema(context, info->schema);
		info->schema = bound_schema->name;
	}
	switch (info->type) {
	case CatalogType::INDEX:
		result = BindCreateIndexInfo(move(info));
		break;
	case CatalogType::TABLE:
		result = BindCreateTableInfo(move(info));
		break;
	case CatalogType::VIEW:
		result = BindCreateViewInfo(move(info));
		break;
	default:
		result = make_unique<BoundCreateInfo>(move(info));
		break;
	}
	result->schema = bound_schema;
	return result;
}
