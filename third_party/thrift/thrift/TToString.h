/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements. See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership. The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

#ifndef _DUCKDB_THRIFT_TOSTRING_H_
#define _DUCKDB_THRIFT_TOSTRING_H_ 1

#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include "duckdb/common/stringstream.hpp"
#include "duckdb/common/vector.hpp"

namespace duckdb_apache {
namespace thrift {

template <typename T>
duckdb::string to_string(const T& t) {
  duckdb::ostringstream o;
  o << t;
  return o.str();
}

// TODO: replace the computations below with std::numeric_limits::max_digits10 once C++11
// is enabled.
inline duckdb::string to_string(const float& t) {
	duckdb::ostringstream o;
  o.precision(static_cast<std::streamsize>(std::ceil(static_cast<double>(std::numeric_limits<float>::digits * std::log10(2.0f) + 1))));
  o << t;
  return o.str();
}

inline duckdb::string to_string(const double& t) {
  duckdb::ostringstream o;
  o.precision(static_cast<std::streamsize>(std::ceil(static_cast<double>(std::numeric_limits<double>::digits * std::log10(2.0f) + 1))));
  o << t;
  return o.str();
}

inline duckdb::string to_string(const long double& t) {
  duckdb::ostringstream o;
  o.precision(static_cast<std::streamsize>(std::ceil(static_cast<double>(std::numeric_limits<long double>::digits * std::log10(2.0f) + 1))));
  o << t;
  return o.str();
}

template <typename K, typename V>
duckdb::string to_string(const std::map<K, V>& m);

template <typename T>
duckdb::string to_string(const std::set<T>& s);

template <typename T>
duckdb::string to_string(const duckdb::vector<T>& t);

template <typename K, typename V>
duckdb::string to_string(const typename std::pair<K, V>& v) {
  duckdb::ostringstream o;
  o << to_string(v.first) << ": " << to_string(v.second);
  return o.str();
}

template <typename T>
duckdb::string to_string(const T& beg, const T& end) {
  duckdb::ostringstream o;
  for (T it = beg; it != end; ++it) {
    if (it != beg)
      o << ", ";
    o << to_string(*it);
  }
  return o.str();
}

template <typename T>
duckdb::string to_string(const duckdb::vector<T>& t) {
  duckdb::ostringstream o;
  o << "[" << to_string(t.begin(), t.end()) << "]";
  return o.str();
}

template <typename K, typename V>
duckdb::string to_string(const std::map<K, V>& m) {
  duckdb::ostringstream o;
  o << "{" << to_string(m.begin(), m.end()) << "}";
  return o.str();
}

template <typename T>
duckdb::string to_string(const std::set<T>& s) {
  duckdb::ostringstream o;
  o << "{" << to_string(s.begin(), s.end()) << "}";
  return o.str();
}
}
} // duckdb_apache::thrift

#endif // _DUCKDB_THRIFT_TOSTRING_H_
