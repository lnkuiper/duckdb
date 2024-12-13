#include "duckdb/common/types/hash.hpp"

#include "duckdb/common/helper.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/common/types/interval.hpp"
#include "duckdb/common/types/uhugeint.hpp"

#include <functional>
#include <cmath>

namespace duckdb {

template <>
hash_t Hash(uint64_t val) {
	return MurmurHash64(val);
}

template <>
hash_t Hash(int64_t val) {
	return MurmurHash64((uint64_t)val);
}

template <>
hash_t Hash(hugeint_t val) {
	return MurmurHash64(val.lower) ^ MurmurHash64(static_cast<uint64_t>(val.upper));
}

template <>
hash_t Hash(uhugeint_t val) {
	return MurmurHash64(val.lower) ^ MurmurHash64(val.upper);
}

template <class T>
struct FloatingPointEqualityTransform {
	static void OP(T &val) {
		if (val == (T)0.0) {
			// Turn negative zero into positive zero
			val = (T)0.0;
		} else if (std::isnan(val)) {
			val = std::numeric_limits<T>::quiet_NaN();
		}
	}
};

template <>
hash_t Hash(float val) {
	static_assert(sizeof(float) == sizeof(uint32_t), "");
	FloatingPointEqualityTransform<float>::OP(val);
	uint32_t uval = Load<uint32_t>(const_data_ptr_cast(&val));
	return MurmurHash64(uval);
}

template <>
hash_t Hash(double val) {
	static_assert(sizeof(double) == sizeof(uint64_t), "");
	FloatingPointEqualityTransform<double>::OP(val);
	uint64_t uval = Load<uint64_t>(const_data_ptr_cast(&val));
	return MurmurHash64(uval);
}

template <>
hash_t Hash(interval_t val) {
	int64_t months, days, micros;
	val.Normalize(months, days, micros);
	return Hash(days) ^ Hash(months) ^ Hash(micros);
}

template <>
hash_t Hash(const char *str) {
	return Hash(str, strlen(str));
}

template <>
hash_t Hash(string_t val) {
	return Hash(val.GetData(), val.GetSize());
}

template <>
hash_t Hash(char *val) {
	return Hash<const char *>(val);
}

hash_t HashBytes(data_ptr_t ptr, const idx_t len) noexcept {
	const auto end = ptr + len;

	hash_t h = 0;
	while (ptr + sizeof(uint64_t) < end) {
		h ^= Hash(Load<uint64_t>(ptr));
		ptr += sizeof(uint64_t);
	}

	hash_t hr = 0;
	memcpy(&hr, ptr, end - ptr);

	return h ^ hr;
}

hash_t Hash(const char *val, size_t size) {
	return HashBytes(data_ptr_cast(val), size);
}

hash_t Hash(uint8_t *val, size_t size) {
	return HashBytes(data_ptr_cast(val), size);
}

} // namespace duckdb
