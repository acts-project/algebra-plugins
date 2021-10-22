/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// Vc include(s).
#include <Vc/Vc>

// System include(s).
#include <cmath>
#include <type_traits>

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
template <template <typename, auto> class array_t, typename scalar_t, auto N,
          std::enable_if_t<std::is_same<array_t<scalar_t, N>,
                                        Vc::SimdArray<scalar_t, N> >::value,
                           bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t dot(const array_t<scalar_t, N> &a,
                                        const array_t<scalar_t, N> &b) {

  return (a * b).sum();
}

/** Get a normalized version of the input vector
 *
 * @tparam vector_type generic input vector type
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t, auto N>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, N> normalize(
    const array_t<scalar_t, N> &v) {

  return v / std::sqrt(dot<array_t>(v, v));
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
template <template <typename, auto> class array_t, typename scalar_t, auto N>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, N> cross(
    const array_t<scalar_t, N> &a, const array_t<scalar_t, N> &b) {
  return {a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0],
          a[0] * b[1] - b[0] * a[1]};
}

}  // namespace algebra::vc::math
