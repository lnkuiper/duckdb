#include "duckdb/main/appender.hpp"

#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/storage/data_table.hpp"

#include "duckdb/common/operator/cast_operators.hpp"

using namespace duckdb;
using namespace std;

Appender::Appender(Connection &con, string schema_name, string table_name) : con(con), column(0) {
	description = con.TableInfo(schema_name, table_name);
	if (!description) {
		// table could not be found
		throw CatalogException(
		    StringUtil::Format("Table \"%s.%s\" could not be found", schema_name.c_str(), table_name.c_str()));
	} else {
		vector<TypeId> types;
		for (auto &column : description->columns) {
			types.push_back(GetInternalType(column.type));
		}
		chunk.Initialize(types);
		con.context->RegisterAppender(this);
	}
}

Appender::Appender(Connection &con, string table_name) : Appender(con, DEFAULT_SCHEMA, table_name) {
}

Appender::~Appender() {
	Close();
}

void Appender::CheckInvalidated() {
	if (!invalidated_msg.empty()) {
		throw Exception("Invalid appender: " + invalidated_msg);
	}
}

void Appender::BeginRow() {
	CheckInvalidated();
	column = 0;
}

void Appender::EndRow() {
	CheckInvalidated();
	// check that all rows have been appended to
	if (column != chunk.column_count) {
		string msg = "Call to EndRow before all rows have been appended to!";
		Invalidate(msg, true);
		throw Exception(msg);
	}
	if (chunk.size() >= STANDARD_VECTOR_SIZE) {
		Flush();
	}
}

template <class SRC, class DST> void Appender::AppendValueInternal(Vector &col, SRC input) {
	((DST *)col.GetData())[col.count++] = Cast::Operation<SRC, DST>(input);
}

template <class T> void Appender::AppendValueInternal(T input) {
	CheckInvalidated();
	if (column >= chunk.column_count) {
		throw Exception("Too many appends for chunk!");
	}
	auto &col = chunk.data[column];
	switch (col.type) {
	case TypeId::BOOL:
		AppendValueInternal<T, bool>(col, input);
		break;
	case TypeId::INT8:
		AppendValueInternal<T, int8_t>(col, input);
		break;
	case TypeId::INT16:
		AppendValueInternal<T, int16_t>(col, input);
		break;
	case TypeId::INT32:
		AppendValueInternal<T, int32_t>(col, input);
		break;
	case TypeId::INT64:
		AppendValueInternal<T, int64_t>(col, input);
		break;
	case TypeId::FLOAT:
		AppendValueInternal<T, float>(col, input);
		break;
	case TypeId::DOUBLE:
		AppendValueInternal<T, double>(col, input);
		break;
	default:
		AppendValue(Value::CreateValue<T>(input));
		return;
	}
	column++;
}

template <> void Appender::Append(bool value) {
	AppendValueInternal<bool>(value);
}

template <> void Appender::Append(int8_t value) {
	AppendValueInternal<int8_t>(value);
}

template <> void Appender::Append(int16_t value) {
	AppendValueInternal<int16_t>(value);
}

template <> void Appender::Append(int32_t value) {
	AppendValueInternal<int32_t>(value);
}

template <> void Appender::Append(int64_t value) {
	AppendValueInternal<int64_t>(value);
}

template <> void Appender::Append(const char *value) {
	AppendValueInternal<const char *>(value);
}

template <> void Appender::Append(double value) {
	AppendValueInternal<double>(value);
}

template <> void Appender::Append(Value value) {
	if (column >= chunk.column_count) {
		throw Exception("Too many appends for chunk!");
	}
	AppendValue(move(value));
}

template <> void Appender::Append(nullptr_t value) {
	if (column >= chunk.column_count) {
		throw Exception("Too many appends for chunk!");
	}
	auto &col = chunk.data[column++];
	col.nullmask[col.count++] = true;
}

void Appender::AppendValue(Value value) {
	chunk.data[column].SetValue(chunk.data[column].count++, value);
	column++;
}

Vector &Appender::GetAppendVector(index_t col_idx) {
	if (col_idx >= chunk.column_count) {
		throw Exception("Column index out of range!");
	}
	return chunk.data[col_idx];
}

void Appender::Flush() {
	if (chunk.size() == 0) {
		return;
	}
	CheckInvalidated();
	try {
		// check that all vectors have the same length before appending
		index_t chunk_size = chunk.size();
		for (index_t i = 1; i < chunk.column_count; i++) {
			if (chunk.data[i].count != chunk_size) {
				throw Exception("Failed to Flush appender: vectors have different number of rows");
			}
		}

		con.Append(*description, chunk);
	} catch (Exception &ex) {
		con.context->RemoveAppender(this);
		Invalidate(ex.what());
		throw ex;
	}
	chunk.Reset();
	column = 0;
}

void Appender::Close() {
	if (!invalidated_msg.empty()) {
		return;
	}
	if (column == 0 || column == chunk.column_count) {
		Flush();
	}
	Invalidate("The appender has been closed!", true);
}

void Appender::Invalidate(string msg, bool close) {
	if (!invalidated_msg.empty()) {
		return;
	}
	if (close) {
		con.context->RemoveAppender(this);
	}
	assert(!msg.empty());
	invalidated_msg = msg;
}
