#include "catch.hpp"
#include "duckdb/common/random_engine.hpp"

using namespace duckdb;

TEST_CASE("Bounded random integers respect 64-bit half-open ranges and reseeding", "[random_engine]") {
	struct Range {
		uint64_t min;
		uint64_t max;
	};
	const auto max_value = NumericLimits<uint64_t>::Maximum();
	const auto limit32 = uint64_t(1) << 32;
	const auto limit63 = uint64_t(1) << 63;
	const Range ranges[] = {{0, 1},
	                        {17, 18},
	                        {0, 3},
	                        {0, limit32 + 1},
	                        {limit32 - 1, limit32 + 2},
	                        {limit63, max_value},
	                        {max_value - 1, max_value},
	                        {max_value - 7, max_value},
	                        {0, max_value}};
	for (auto range : ranges) {
		CAPTURE(range.min, range.max);
		RandomEngine random(42);
		uint64_t values[128];
		bool has_large_value = false;
		for (auto &value : values) {
			value = random.NextRandomInteger64(range.min, range.max);
			REQUIRE(value >= range.min);
			REQUIRE(value < range.max);
			has_large_value |= value > NumericLimits<uint32_t>::Maximum();
		}
		if (range.max == max_value) {
			REQUIRE(has_large_value);
		}
		random.SetSeed(42);
		for (auto value : values) {
			REQUIRE(random.NextRandomInteger64(range.min, range.max) == value);
		}
	}
}

TEST_CASE("Bounded random integers avoid modulo bias", "[random_engine]") {
	RandomEngine random(42);
	const auto third = uint64_t(1) << 62;
	const idx_t sample_count = 8192;
	idx_t counts[3] = {};
	for (idx_t i = 0; i < sample_count; i++) {
		auto value = random.NextRandomInteger64(0, 3 * third);
		REQUIRE(value < 3 * third);
		counts[value / third]++;
	}
	// Naive modulo maps half the draws into the first third instead of one third.
	for (auto count : counts) {
		REQUIRE(count > sample_count / 3 - sample_count / 20);
		REQUIRE(count < sample_count / 3 + sample_count / 20);
	}
}
