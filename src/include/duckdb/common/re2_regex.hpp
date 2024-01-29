// RE2 compatibility layer with std::regex

#pragma once

#include "duckdb/common/string.hpp"
#include "duckdb/common/vector.hpp"
#include "duckdb/common/winapi.hpp"

#include <stdexcept>

namespace duckdb_re2 {
class RE2;

enum class RegexOptions : uint8_t { NONE, CASE_INSENSITIVE };

class Regex {
public:
	DUCKDB_API Regex(const duckdb::string &pattern, RegexOptions options = RegexOptions::NONE);
	Regex(const char *pattern, RegexOptions options = RegexOptions::NONE) : Regex(duckdb::string(pattern)) {
	}
	const duckdb_re2::RE2 &GetRegex() const {
		return *regex;
	}

private:
	std::shared_ptr<duckdb_re2::RE2> regex;
};

struct GroupMatch {
	duckdb::string text;
	uint32_t position;

	const duckdb::string &str() const {
		return text;
	}
	operator duckdb::string() const {
		return text;
	}
};

struct Match {
	duckdb::vector<GroupMatch> groups;

	GroupMatch &GetGroup(uint64_t index) {
		if (index >= groups.size()) {
			throw std::runtime_error("RE2: Match index is out of range");
		}
		return groups[index];
	}

	duckdb::string str(uint64_t index) {
		return GetGroup(index).text;
	}

	uint64_t position(uint64_t index) {
		return GetGroup(index).position;
	}

	uint64_t length(uint64_t index) {
		return GetGroup(index).text.size();
	}

	GroupMatch &operator[](uint64_t i) {
		return GetGroup(i);
	}
};

DUCKDB_API bool RegexSearch(const duckdb::string &input, Match &match, const Regex &regex);
DUCKDB_API bool RegexMatch(const duckdb::string &input, Match &match, const Regex &regex);
DUCKDB_API bool RegexMatch(const char *start, const char *end, Match &match, const Regex &regex);
DUCKDB_API bool RegexMatch(const duckdb::string &input, const Regex &regex);
DUCKDB_API duckdb::vector<Match> RegexFindAll(const duckdb::string &input, const Regex &regex);

} // namespace duckdb_re2
