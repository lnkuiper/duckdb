//===----------------------------------------------------------------------===//
//                         DuckDB
//
// mbedtls_wrapper.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/optional_ptr.hpp"
#include "duckdb/common/common.hpp"
#include "duckdb/common/string.hpp"

namespace duckdb_mbedtls {

class MbedTlsWrapper {
public:
	static void ComputeSha256Hash(const char *in, size_t in_len, char *out);
	static duckdb::string ComputeSha256Hash(const duckdb::string &file_content);
	static bool IsValidSha256Signature(const duckdb::string &pubkey, const duckdb::string &signature,
	                                   const duckdb::string &sha256_hash);
	static void Hmac256(const char *key, size_t key_len, const char *message, size_t message_len, char *out);
	static void ToBase16(char *in, char *out, size_t len);

	DUCKDB_API static void GenerateRandomData(duckdb::data_ptr_t data, duckdb::idx_t len);

	static constexpr size_t SHA256_HASH_LENGTH_BYTES = 32;
	static constexpr size_t SHA256_HASH_LENGTH_TEXT = 64;

	class SHA256State {
	public:
		SHA256State();
		~SHA256State();
		void AddString(const duckdb::string &str);
		duckdb::string Finalize();
		void FinishHex(char *out);

	private:
		void *sha_context;
	};

	class AESGCMState {
	public:
		DUCKDB_API AESGCMState(const duckdb::string &key);
		DUCKDB_API ~AESGCMState();

	public:
		DUCKDB_API static bool ValidKey(const duckdb::string &key);
		DUCKDB_API void InitializeEncryption(duckdb::const_data_ptr_t iv, duckdb::idx_t iv_len);
		DUCKDB_API void InitializeDecryption(duckdb::const_data_ptr_t iv, duckdb::idx_t iv_len);
		DUCKDB_API size_t Process(duckdb::const_data_ptr_t in, duckdb::idx_t in_len, duckdb::data_ptr_t out,
		               duckdb::idx_t out_len);
		DUCKDB_API size_t Finalize(duckdb::data_ptr_t out, duckdb::idx_t out_len, duckdb::data_ptr_t tag, duckdb::idx_t tag_len);

	public:
		static constexpr size_t BLOCK_SIZE = 16;

	private:
		void *gcm_context;
	};
};

} // namespace duckdb_mbedtls
