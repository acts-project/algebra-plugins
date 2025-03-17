/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/fastor_fastor.hpp"

// System include(s)
#include <algorithm>
#include <vector>

namespace algebra {

/// Fill an @c Fastor based vector with random values
template <concepts::vector vector_t>
inline void fill_random_vec(std::vector<vector_t>& collection) {

  auto rand_obj = [](vector_t& v) { v.random(); };

  collection.resize(collection.capacity());
  std::ranges::for_each(collection, rand_obj);
}

/// Fill a @c Fastor based transform3 with random values
template <concepts::transform3D transform3_t>
inline void fill_random_trf(std::vector<transform3_t>& collection) {

  using vector_t = typename transform3_t::vector3;

  auto rand_obj = [](transform3_t& trf) {
    vector_t x_axis;
    vector_t z_axis;
    vector_t t;

    x_axis.random();
    x_axis = vector::normalize(x_axis);
    z_axis.random();
    t.random();
    t = vector::normalize(t);

    // Gram-Schmidt projection
    typename transform3_t::scalar_type coeff =
        vector::dot(x_axis, z_axis) / vector::norm(x_axis);
    z_axis = x_axis - coeff * z_axis;

    trf = transform3_t{t, x_axis, vector::normalize(z_axis)};
  };

  collection.resize(collection.capacity());
  std::ranges::for_each(collection, rand_obj);
}

/// Fill a @c Fastor based matrix with random values
template <concepts::matrix matrix_t>
inline void fill_random_matrix(std::vector<matrix_t>& collection) {

  auto rand_obj = [](matrix_t& m) { m.random(); };

  collection.resize(collection.capacity());
  std::ranges::for_each(collection, rand_obj);
}

}  // namespace algebra
