/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"
#include "algebra/storage/vc.hpp"

// System include(s).
#include <cmath>

namespace algebra::vc::math {

/** Dot product between two input vectors
 *
 * @tparam vector_type generic input vector type
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST_DEVICE inline scalar dot(const storage_type<scalar_t, N> &a,
                                      const storage_type<scalar_t, N> &b) {

  return (a * b).sum();
}

/** Get a normalized version of the input vector
 *
 * @tparam vector_type generic input vector type
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST_DEVICE inline storage_type<scalar_t, N> normalize(
    const storage_type<scalar_t, N> &v) {

  return v / std::sqrt(dot(v, v));
}

/** Cross product between two input vectors - 3 Dim
 *
 * @tparam vector_type generic input vector type
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector representing the cross product
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST_DEVICE inline storage_type<scalar_t, N> cross(
    const storage_type<scalar_t, N> &a, const storage_type<scalar_t, N> &b) {
  return {a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0],
          a[0] * b[1] - b[0] * a[1]};
}

}  // namespace algebra::vc::math
