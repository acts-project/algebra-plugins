/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <algorithm>
#include <vector>

namespace algebra {

// @todo Leave for later

/// Fill a @c Vc::SimdArray based vector with random values
/*template <typename vector_t>
inline void fill_random(
    std::vector<vector3_s<value_t>, allocator_t<vector3_s<value_t>>>
        &collection) {

  using vector_t = vector3_s<value_t>;
  // Generate a vector of the right type with random values
  auto rand_obj = [&]() { return vector_t{vector_t::array_type::Random()}; };

  std::generate(collection.begin(), collection.end(), rand_obj);
}*/

/// Fill a @c Vc::Vector based vector with random values
template <typename vector_soa_t>
inline void fill_random(std::vector<vector_soa_t> &collection) {
  // Generate a vector of the right type with random values
  auto rand_obj = []() {
    using simd_vector_t = typename vector_soa_t::value_type;
    vector_soa_t tmp{};
    tmp[0] = simd_vector_t::Random();
    tmp[1] = simd_vector_t::Random();
    tmp[2] = simd_vector_t::Random();
    return tmp;
  };

  collection.resize(collection.capacity());
  std::generate(collection.begin(), collection.end(), rand_obj);
}

}  // namespace algebra