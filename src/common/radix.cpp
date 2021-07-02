#include "duckdb/common/radix.hpp"

#include <cfloat>
#include <cstring> // strlen() on Solaris
#include <limits.h>

namespace duckdb {

bool IsLittleEndian() {
	int n = 1;
	if (*(char *)&n == 1) {
		return true;
	} else {
		return false;
	}
}

uint8_t FlipSign(uint8_t key_byte) {
	return key_byte ^ 128;
}

uint32_t EncodeFloat(float x) {
	uint64_t buff;

	//! zero
	if (x == 0) {
		buff = 0;
		buff |= (1u << 31);
		return buff;
	}
	//! infinity
	if (x > FLT_MAX) {
		return UINT_MAX;
	}
	//! -infinity
	if (x < -FLT_MAX) {
		return 0;
	}
	buff = Load<uint32_t>((const_data_ptr_t)&x);
	if ((buff & (1u << 31)) == 0) { //! +0 and positive numbers
		buff |= (1u << 31);
	} else {          //! negative numbers
		buff = ~buff; //! complement 1
	}

	return buff;
}

uint64_t EncodeDouble(double x) {
	uint64_t buff;
	//! zero
	if (x == 0) {
		buff = 0;
		buff += (1ull << 63);
		return buff;
	}
	//! infinity
	if (x > DBL_MAX) {
		return ULLONG_MAX;
	}
	//! -infinity
	if (x < -DBL_MAX) {
		return 0;
	}
	buff = Load<uint64_t>((const_data_ptr_t)&x);
	if (buff < (1ull << 63)) { //! +0 and positive numbers
		buff += (1ull << 63);
	} else {          //! negative numbers
		buff = ~buff; //! complement 1
	}
	return buff;
}

template <>
void EncodeData(data_ptr_t dataptr, bool value, bool is_little_endian) {
	Store<uint8_t>(value ? 1 : 0, dataptr);
}

template <>
void EncodeData(data_ptr_t dataptr, int8_t value, bool is_little_endian) {
	Store<uint8_t>(value, dataptr);
	dataptr[0] = FlipSign(dataptr[0]);
}

template <>
void EncodeData(data_ptr_t dataptr, int16_t value, bool is_little_endian) {
	Store<uint16_t>(is_little_endian ? BSWAP16(value) : value, dataptr);
	dataptr[0] = FlipSign(dataptr[0]);
}

template <>
void EncodeData(data_ptr_t dataptr, int32_t value, bool is_little_endian) {
	Store<uint32_t>(is_little_endian ? BSWAP32(value) : value, dataptr);
	dataptr[0] = FlipSign(dataptr[0]);
}

template <>
void EncodeData(data_ptr_t dataptr, int64_t value, bool is_little_endian) {
	Store<uint64_t>(is_little_endian ? BSWAP64(value) : value, dataptr);
	dataptr[0] = FlipSign(dataptr[0]);
}

template <>
void EncodeData(data_ptr_t dataptr, uint8_t value, bool is_little_endian) {
	Store<uint8_t>(value, dataptr);
}

template <>
void EncodeData(data_ptr_t dataptr, uint16_t value, bool is_little_endian) {
	Store<uint16_t>(is_little_endian ? BSWAP16(value) : value, dataptr);
}

template <>
void EncodeData(data_ptr_t dataptr, uint32_t value, bool is_little_endian) {
	Store<uint32_t>(is_little_endian ? BSWAP32(value) : value, dataptr);
}

template <>
void EncodeData(data_ptr_t dataptr, uint64_t value, bool is_little_endian) {
	Store<uint64_t>(is_little_endian ? BSWAP64(value) : value, dataptr);
}

template <>
void EncodeData(data_ptr_t dataptr, hugeint_t value, bool is_little_endian) {
	EncodeData<int64_t>(dataptr, value.upper, is_little_endian);
	EncodeData<uint64_t>(dataptr + sizeof(value.upper), value.lower, is_little_endian);
}

template <>
void EncodeData(data_ptr_t dataptr, float value, bool is_little_endian) {
	uint32_t converted_value = EncodeFloat(value);
	Store<uint32_t>(is_little_endian ? BSWAP32(converted_value) : converted_value, dataptr);
}

template <>
void EncodeData(data_ptr_t dataptr, double value, bool is_little_endian) {
	uint64_t converted_value = EncodeDouble(value);
	Store<uint64_t>(is_little_endian ? BSWAP64(converted_value) : converted_value, dataptr);
}

template <>
void EncodeData(data_ptr_t dataptr, interval_t value, bool is_little_endian) {
	EncodeData<int32_t>(dataptr, value.months, is_little_endian);
	dataptr += sizeof(value.months);
	EncodeData<int32_t>(dataptr, value.days, is_little_endian);
	dataptr += sizeof(value.days);
	EncodeData<int64_t>(dataptr, value.micros, is_little_endian);
}

void EncodeStringData(data_ptr_t dataptr, string_t value, idx_t prefix_len) {
	auto len = value.GetSize();
	memcpy(dataptr, value.GetDataUnsafe(), MinValue(len, prefix_len));
	if (len < prefix_len) {
		memset(dataptr + len, '\0', prefix_len - len);
	}
}

float DecodeFloat(uint32_t buff) {
	//! zero
	if (buff == (1u << 31)) {
		return 0;
	}
	//! infinity
	if (buff == UINT_MAX) {
		return FLT_MAX + FLT_EPSILON;
	}
	//! -infinity
	if (buff == 0) {
		return -(FLT_MAX + FLT_EPSILON);
	}
	if ((buff & (1u << 31)) == 1) { //! +0 and positive numbers
		buff &= ~(1u << 31);
	} else {          //! negative numbers
		buff = ~buff; //! complement 1
	}
	return Load<uint32_t>((const_data_ptr_t)&buff);
}

double DecodeDouble(uint64_t buff) {
	//! zero
	if (buff == (1ull << 63)) {
		return 0;
	}
	//! infinity
	if (buff == ULLONG_MAX) {
		return DBL_MAX + DBL_EPSILON;
	}
	//! -infinity
	if (buff == 0) {
		return -(DBL_MAX + DBL_EPSILON);
	}
	if (buff >= (1ull << 63)) { //! +0 and positive numbers
		buff -= (1ull << 63);
	} else {          //! negative numbers
		buff = ~buff; //! complement 1
	}
	return Load<uint64_t>((const_data_ptr_t)&buff);
}

template <>
bool DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	return Load<uint8_t>(dataptr);
}

template <>
int8_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	dataptr[0] = FlipSign(dataptr[0]);
	return Load<uint8_t>(dataptr);
}

template <>
int16_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	dataptr[0] = FlipSign(dataptr[0]);
	auto value = Load<uint16_t>(dataptr);
	return is_little_endian ? BSWAP16(value) : value;
}

template <>
int32_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	dataptr[0] = FlipSign(dataptr[0]);
	auto value = Load<uint32_t>(dataptr);
	return is_little_endian ? BSWAP32(value) : value;
}

template <>
int64_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	dataptr[0] = FlipSign(dataptr[0]);
	auto value = Load<uint64_t>(dataptr);
	return is_little_endian ? BSWAP64(value) : value;
}

template <>
uint8_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	return Load<uint8_t>(dataptr);
}

template <>
uint16_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	auto value = Load<uint16_t>(dataptr);
	return is_little_endian ? BSWAP16(value) : value;
}

template <>
uint32_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	auto value = Load<uint32_t>(dataptr);
	return is_little_endian ? BSWAP32(value) : value;
}

template <>
uint64_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	auto value = Load<uint64_t>(dataptr);
	return is_little_endian ? BSWAP64(value) : value;
}

template <>
hugeint_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	hugeint_t value;
	value.upper = DecodeData<int64_t>(dataptr, is_little_endian);
	value.lower = DecodeData<uint64_t>(dataptr + sizeof(value.upper), is_little_endian);
	return value;
}

template <>
float DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	auto value = Load<uint32_t>(dataptr);
	return DecodeFloat(is_little_endian ? BSWAP32(value) : value);
}

template <>
double DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	auto value = Load<uint64_t>(dataptr);
	return DecodeDouble(is_little_endian ? BSWAP64(value) : value);
}

template <>
interval_t DecodeData(data_ptr_t dataptr, bool is_little_endian) {
	interval_t value;
	value.months = DecodeData<int32_t>(dataptr, is_little_endian);
	dataptr += sizeof(value.months);
	value.days = DecodeData<int32_t>(dataptr, is_little_endian);
	dataptr += sizeof(value.days);
	value.micros = DecodeData<int64_t>(dataptr, is_little_endian);
	return value;
}

string_t DecodeStringData(data_ptr_t dataptr, idx_t prefix_len) {
	uint32_t len;
	for (len = 0; len < prefix_len; len++) {
		if (dataptr[0] == '\0') {
			break;
		}
	}
	return string_t((const char *)dataptr, len);
}

} // namespace duckdb
