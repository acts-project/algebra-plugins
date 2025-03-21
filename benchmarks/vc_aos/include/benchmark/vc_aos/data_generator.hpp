/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/vc_aos.hpp"

// System include(s)
#include <algorithm>
#include <vector>

namespace algebra {

/// Fill a @c Vc::SimdArray based vector with random values
template <concepts::vector vector_aos_t>
inline void fill_random_vec(std::vector<vector_aos_t> &collection) {

  // Generate a vector of the right type with random values
  auto rand_obj = [&]() {
    return vector_aos_t{vector_aos_t::array_type::Random()};
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

/// Fill a @c Vc::SimdArray based transform3 with random values
template <concepts::transform3D transform3_t>
inline void fill_random_trf(std::vector<transform3_t> &collection) {
  // Generate a random, but valid affine transformation
  auto rand_obj = []() {
    using vector_t = typename transform3_t::vector3;
    vector_t x_axis;
    vector_t z_axis;
    vector_t t;

    x_axis = vector_t{vector_t::array_type::Random()};
    x_axis = vector::normalize(x_axis);

    z_axis = vector_t{vector_t::array_type::Random()};
    z_axis = vector::normalize(z_axis);

    t = vector_t{vector_t::array_type::Random()};
    t = vector::normalize(t);

    // Gram-Schmidt projection
    typename transform3_t::scalar_type coeff =
        vector::dot(x_axis, z_axis) / vector::norm(x_axis);
    z_axis = x_axis - coeff * z_axis;

    return transform3_t{t, x_axis, vector::normalize(z_axis)};
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

/// Fill a @c Vc::SimdArray based matrix with random values
template <concepts::matrix matrix_t>
inline void fill_random_matrix(std::vector<matrix_t> &collection) {
  // Generate a random, but valid affine transformation
  auto rand_obj = []() {
    using vector_t = typename matrix_t::vector_type;
    vector_t x_axis;
    vector_t z_axis;
    vector_t t;

    matrix_t m;
    for (std::size_t j = 0u; j < matrix_t::columns(); ++j) {
      m[j] = vector_t::array_type::Random();
    }

    return m;
  };

  collection.resize(collection.capacity());
  std::ranges::generate(collection, rand_obj);
}

}  // namespace algebra
