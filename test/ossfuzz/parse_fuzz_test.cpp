#include "duckdb.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/unordered_set.hpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
	duckdb::string input(reinterpret_cast<const char *>(data), size);
	duckdb::DuckDB db(nullptr);
	duckdb::Connection con(db);

	duckdb::unordered_set<duckdb::string> internal_error_messages = {"Unoptimized Result differs from original result!",
	                                                                 "INTERNAL"};
	con.Query("PRAGMA enable_verification");
	try {
		auto result = con.Query(input);
		if (result->HasError()) {
			for (auto &internal_error : internal_error_messages) {
				if (duckdb::StringUtil::Contains(result->GetError(), internal_error)) {
					return 1;
				}
			}
		}
	} catch (std::exception &e) {
	}
	return 0;
}
