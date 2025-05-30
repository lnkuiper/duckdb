#include "duckdb/common/random_engine.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "pcg_random.hpp"

#ifdef __linux__
#include <sys/syscall.h>
#include <unistd.h>
#else
#include <random>
#endif
namespace duckdb {

struct RandomState {
	RandomState() {
	}

	pcg32 pcg;
};

RandomEngine::RandomEngine(int64_t seed) : random_state(make_uniq<RandomState>()) {
	if (seed < 0) {
#ifdef __linux__
		idx_t random_seed = 0;
		int result = -1;
#if defined(SYS_getrandom)
		result = static_cast<int>(syscall(SYS_getrandom, &random_seed, sizeof(random_seed), 0));
#endif
		if (result == -1) {
			// Something went wrong with the syscall, we use chrono
			const auto now = std::chrono::high_resolution_clock::now();
			random_seed = now.time_since_epoch().count();
		}
		random_state->pcg.seed(random_seed);
#else
		random_state->pcg.seed(pcg_extras::seed_seq_from<std::random_device>());
#endif
	} else {
		random_state->pcg.seed(NumericCast<uint64_t>(seed));
	}
}

RandomEngine::~RandomEngine() {
}

double RandomEngine::NextRandom(double min, double max) {
	D_ASSERT(max >= min);
	return min + (NextRandom() * (max - min));
}

double RandomEngine::NextRandom() {
	auto uint64 = NextRandomInteger64();
	return std::ldexp(uint64, -64);
}

double RandomEngine::NextRandom32(double min, double max) {
	D_ASSERT(max >= min);
	return min + (NextRandom32() * (max - min));
}

double RandomEngine::NextRandom32() {
	auto uint32 = NextRandomInteger();
	return std::ldexp(uint32, -32);
}

uint32_t RandomEngine::NextRandomInteger() {
	return random_state->pcg();
}

uint64_t RandomEngine::NextRandomInteger64() {
	return (static_cast<uint64_t>(NextRandomInteger()) << UINT64_C(32)) | static_cast<uint64_t>(NextRandomInteger());
}

uint32_t RandomEngine::NextRandomInteger(uint32_t min, uint32_t max) {
	return min + static_cast<uint32_t>(NextRandom() * double(max - min));
}

uint32_t RandomEngine::NextRandomInteger32(uint32_t min, uint32_t max) {
	return min + static_cast<uint32_t>(NextRandom32() * double(max - min));
}

void RandomEngine::SetSeed(uint64_t seed) {
	random_state->pcg.seed(seed);
}

} // namespace duckdb
