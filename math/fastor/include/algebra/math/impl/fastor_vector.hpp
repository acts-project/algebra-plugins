/** Algebra plugins, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// Fastor include(s).
#include <Fastor/Fastor.h>

namespace algebra::fastor::math {

/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST inline Fastor::Tensor<scalar_t, N> normalize(
    const Fastor::Tensor<scalar_t, N> &v) {

  return Fastor::norm(v);
}

/** Dot product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST_DEVICE inline auto dot(const Fastor::Tensor<scalar_t, N> &a,
                                    const Fastor::Tensor<scalar_t, N> &b) {
  return Fastor::inner(a, b);
}

/** Cross product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector (expression) representing the cross product
 **/
template <typename scalar_t>
ALGEBRA_HOST_DEVICE inline auto cross(const Fastor::Tensor<scalar_t, 3> &a,
                                      const Fastor::Tensor<scalar_t, 3> &b) {
  return Fastor::cross(a, b);
}

}  // namespace algebra::fastor::math
