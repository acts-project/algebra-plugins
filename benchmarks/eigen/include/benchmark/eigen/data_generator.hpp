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

namespace algebra {

/// Fill a @c Eigen3 based vector with random values
template <typename vector_t>
inline void fill_random(std::vector<vector_t> &collection) {

  auto rand_obj = []() { return vector_t::Random(); };

  collection.resize(collection.capacity());
  std::generate(collection.begin(), collection.end(), rand_obj);
}

}  // namespace algebra