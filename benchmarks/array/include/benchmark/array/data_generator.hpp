/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <algorithm>
#include <random>
#include <vector>

namespace algebra {

/// Fill an @c std::array based vector with random values
template <typename vector_t>
inline void fill_random(std::vector<vector_t> &collection) {

  // Generate a vector of the right type with random values
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<typename vector_t::value_type> dist(0.f, 1.f);

  auto rand_obj = [&]() { return vector_t{dist(mt), dist(mt), dist(mt)}; };

  collection.resize(collection.capacity());
  std::generate(collection.begin(), collection.end(), rand_obj);
}

}  // namespace algebra