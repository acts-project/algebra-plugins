/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// System include(s).
#include <cmath>

namespace algebra::cmath {

/** Cross product between two input vectors - 3 Dim
 *
 * @tparam derived_type_lhs is the first matrix (epresseion) template
 * @tparam derived_type_rhs is the second matrix (epresseion) template
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector (expression) representing the cross product
 **/
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t, size_type N, std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, N> cross(
    const array_t<scalar_t, N> &a, const array_t<scalar_t, N> &b) {

  return {a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0],
          a[0] * b[1] - b[0] * a[1]};
}

/** Dot product between two input vectors - 2 Dim
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t, size_type N, std::enable_if_t<N == 2, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t dot(const array_t<scalar_t, N> &a,
                                        const array_t<scalar_t, N> &b) {

  return a[0] * b[0] + a[1] * b[1];
}

/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t, size_type N, std::enable_if_t<N == 2, bool> = true>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, N> normalize(
    const array_t<scalar_t, N> &v) {

  scalar_t oon = 1. / std::sqrt(dot(v, v));
  return {v[0] * oon, v[1] * oon};
}

/** Dot product between two input vectors - 3 Dim
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t, size_type N1, size_type N2,
          std::enable_if_t<N1 >= 3, bool> = true,
          std::enable_if_t<N2 >= 3, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t dot(const array_t<scalar_t, N1> &a,
                                        const array_t<scalar_t, N2> &b) {

  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t, size_type N, std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> normalize(
    const array_t<scalar_t, N> &v) {

  scalar_t oon = 1. / std::sqrt(dot(v, v));
  return {v[0] * oon, v[1] * oon, v[2] * oon};
}

}  // namespace algebra::cmath
