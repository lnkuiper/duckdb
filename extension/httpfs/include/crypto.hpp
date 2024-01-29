#pragma once

#include "duckdb/common/string.hpp"

#include <stddef.h>

namespace duckdb {

typedef unsigned char hash_bytes[32];
typedef unsigned char hash_str[64];

void sha256(const char *in, size_t in_len, hash_bytes &out);

void hmac256(const string &message, const char *secret, size_t secret_len, hash_bytes &out);

void hmac256(string message, hash_bytes secret, hash_bytes &out);

void hex256(hash_bytes &in, hash_str &out);

} // namespace duckdb
