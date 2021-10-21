/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// System include(s).
#include <cmath>
#include <cstddef>
#include <type_traits>

namespace algebra::cmath {

/** This method retrieves phi from a vector, vector base with rows >= 2
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t, auto N,
          std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t phi(
    const array_t<scalar_t, N> &v) noexcept {

  return std::atan2(v[1], v[0]);
}

/** This method retrieves theta from a vector, vector base with rows >= 3
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t theta(
    const array_t<scalar_t, 3> &v) noexcept {

  return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
}

/** This method retrieves the perpenticular magnitude of a vector with rows >= 2
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t, auto N,
          std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t perp(
    const array_t<scalar_t, N> &v) noexcept {

  return std::sqrt(v[0] * v[0] + v[1] * v[1]);
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t norm(const array_t<scalar_t, 2> &v) {

  return perp<array_t>(v);
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t norm(const array_t<scalar_t, 3> &v) {

  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/** This method retrieves the pseudo-rapidity from a vector or vector base with
 *rows >= 3
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t eta(
    const array_t<scalar_t, 3> &v) noexcept {

  return std::atanh(v[2] / norm<array_t>(v));
}

}  // namespace algebra::cmath
