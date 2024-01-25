//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/common/serializer/serialization_traits.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/common/map.hpp"
#include "duckdb/common/optional_ptr.hpp"
#include "duckdb/common/set.hpp"
#include "duckdb/common/shared_ptr.hpp"
#include "duckdb/common/unique_ptr.hpp"
#include "duckdb/common/unordered_map.hpp"
#include "duckdb/common/unordered_set.hpp"
#include "duckdb/common/vector.hpp"

#include <atomic>
#include <cstdint>
#include <type_traits>

namespace duckdb {

class Serializer;   // Forward declare
class Deserializer; // Forward declare

typedef uint16_t field_id_t;
const field_id_t MESSAGE_TERMINATOR_FIELD_ID = 0xFFFF;

// Backport to c++11
template <class...>
using void_t = void;

// Check for anything implementing a `void Serialize(Serializer &Serializer)` method
template <typename T, typename = T>
struct has_serialize : std::false_type {};
template <typename T>
struct has_serialize<
    T, typename std::enable_if<
           std::is_same<decltype(std::declval<T>().Serialize(std::declval<Serializer &>())), void>::value, T>::type>
    : std::true_type {};

template <typename T, typename = T>
struct has_deserialize : std::false_type {};

// Accept `static unique_ptr<T> Deserialize(Deserializer& deserializer)`
template <typename T>
struct has_deserialize<
    T, typename std::enable_if<std::is_same<decltype(T::Deserialize), unique_ptr<T>(Deserializer &)>::value, T>::type>
    : std::true_type {};

// Accept `static shared_ptr<T> Deserialize(Deserializer& deserializer)`
template <typename T>
struct has_deserialize<
    T, typename std::enable_if<std::is_same<decltype(T::Deserialize), shared_ptr<T>(Deserializer &)>::value, T>::type>
    : std::true_type {};

// Accept `static T Deserialize(Deserializer& deserializer)`
template <typename T>
struct has_deserialize<
    T, typename std::enable_if<std::is_same<decltype(T::Deserialize), T(Deserializer &)>::value, T>::type>
    : std::true_type {};

// Check if T is a vector, and provide access to the inner type
template <typename T>
struct is_vector : std::false_type {};
template <typename T>
struct is_vector<typename duckdb::vector<T>> : std::true_type {
	typedef T ELEMENT_TYPE;
};

template <typename T>
struct is_unsafe_vector : std::false_type {};
template <typename T>
struct is_unsafe_vector<typename duckdb::unsafe_vector<T>> : std::true_type {
	typedef T ELEMENT_TYPE;
};

// Check if T is a unordered map, and provide access to the inner type
template <typename T>
struct is_unordered_map : std::false_type {};

template <typename K, typename V>
struct is_unordered_map<typename duckdb::unordered_map<K, V>> : std::true_type {
	typedef K KEY_TYPE;
	typedef V VALUE_TYPE;
	typedef typename std::hash<KEY_TYPE> HASH_TYPE;
	typedef typename std::equal_to<KEY_TYPE> EQUAL_TYPE;
};

template <typename K, typename V, typename H, typename E>
struct is_unordered_map<typename duckdb::unordered_map<K, V, H, E>> : std::true_type {
	typedef K KEY_TYPE;
	typedef V VALUE_TYPE;
	typedef H HASH_TYPE;
	typedef E EQUAL_TYPE;
};

template <typename T>
struct is_map : std::false_type {};

template <typename K, typename V>
struct is_map<typename duckdb::map<K, V>> : std::true_type {
	typedef K KEY_TYPE;
	typedef V VALUE_TYPE;
	typedef typename std::less<KEY_TYPE> COMPARE_TYPE;
};

template <typename K, typename V, typename C>
struct is_map<typename duckdb::map<K, V, C>> : std::true_type {
	typedef K KEY_TYPE;
	typedef V VALUE_TYPE;
	typedef C COMPARE_TYPE;
};

template <typename T>
struct is_unique_ptr : std::false_type {};
template <typename T>
struct is_unique_ptr<unique_ptr<T>> : std::true_type {
	typedef T ELEMENT_TYPE;
};

template <typename T>
struct is_shared_ptr : std::false_type {};
template <typename T>
struct is_shared_ptr<shared_ptr<T>> : std::true_type {
	typedef T ELEMENT_TYPE;
};

template <typename T>
struct is_optional_ptr : std::false_type {};
template <typename T>
struct is_optional_ptr<optional_ptr<T>> : std::true_type {
	typedef T ELEMENT_TYPE;
};

template <typename T>
struct is_pair : std::false_type {};
template <typename T, typename U>
struct is_pair<std::pair<T, U>> : std::true_type {
	typedef T FIRST_TYPE;
	typedef U SECOND_TYPE;
};

template <typename T>
struct is_unordered_set : std::false_type {};

template <typename T>
struct is_unordered_set<typename duckdb::unordered_set<T>> : std::true_type {
	using ELEMENT_TYPE = T;
	typedef typename std::hash<ELEMENT_TYPE> HASH_TYPE;
	typedef typename std::equal_to<ELEMENT_TYPE> EQUAL_TYPE;
};

template <typename T, typename H, typename E>
struct is_unordered_set<typename duckdb::unordered_set<T, H, E>> : std::true_type {
	typedef T ELEMENT_TYPE;
	typedef H HASH_TYPE;
	typedef E EQUAL_TYPE;
};

template <typename T>
struct is_set : std::false_type {};

template <typename T>
struct is_set<typename duckdb::set<T>> : std::true_type {
	typedef T ELEMENT_TYPE;
	typedef typename std::less<T> COMPARE_TYPE;
};

template <typename T, typename C>
struct is_set<typename duckdb::set<T, C>> : std::true_type {
	typedef T ELEMENT_TYPE;
	typedef C COMPARE_TYPE;
};

template <typename T>
struct is_atomic : std::false_type {};

template <typename T>
struct is_atomic<std::atomic<T>> : std::true_type {
	typedef T TYPE;
};

struct SerializationDefaultValue {
	template <typename T = void>
	static inline typename std::enable_if<is_atomic<T>::value, T>::type GetDefault() {
		using INNER = typename is_atomic<T>::TYPE;
		return static_cast<T>(GetDefault<INNER>());
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<is_atomic<T>::value, T>::type &value) {
		using INNER = typename is_atomic<T>::TYPE;
		return value == GetDefault<INNER>();
	}

	template <typename T = void>
	static inline typename std::enable_if<std::is_arithmetic<T>::value, T>::type GetDefault() {
		return static_cast<T>(0);
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<std::is_arithmetic<T>::value, T>::type &value) {
		return value == static_cast<T>(0);
	}

	template <typename T = void>
	static inline typename std::enable_if<is_unique_ptr<T>::value, T>::type GetDefault() {
		return T();
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<is_unique_ptr<T>::value, T>::type &value) {
		return !value;
	}

	template <typename T = void>
	static inline typename std::enable_if<is_optional_ptr<T>::value, T>::type GetDefault() {
		return T();
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<is_optional_ptr<T>::value, T>::type &value) {
		return !value;
	}

	template <typename T = void>
	static inline typename std::enable_if<is_shared_ptr<T>::value, T>::type GetDefault() {
		return T();
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<is_shared_ptr<T>::value, T>::type &value) {
		return !value;
	}

	template <typename T = void>
	static inline typename std::enable_if<is_vector<T>::value, T>::type GetDefault() {
		return T();
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<is_vector<T>::value, T>::type &value) {
		return value.empty();
	}

	template <typename T = void>
	static inline typename std::enable_if<is_unsafe_vector<T>::value, T>::type GetDefault() {
		return T();
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<is_unsafe_vector<T>::value, T>::type &value) {
		return value.empty();
	}

	template <typename T = void>
	static inline typename std::enable_if<is_unordered_set<T>::value, T>::type GetDefault() {
		return T();
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<is_unordered_set<T>::value, T>::type &value) {
		return value.empty();
	}

	template <typename T = void>
	static inline typename std::enable_if<is_unordered_map<T>::value, T>::type GetDefault() {
		return T();
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<is_unordered_map<T>::value, T>::type &value) {
		return value.empty();
	}

	template <typename T = void>
	static inline typename std::enable_if<is_map<T>::value, T>::type GetDefault() {
		return T();
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<is_map<T>::value, T>::type &value) {
		return value.empty();
	}

	template <typename T = void>
	static inline typename std::enable_if<std::is_same<T, string>::value, T>::type GetDefault() {
		return T();
	}

	template <typename T = void>
	static inline bool IsDefault(const typename std::enable_if<std::is_same<T, string>::value, T>::type &value) {
		return value.empty();
	}
};

} // namespace duckdb
