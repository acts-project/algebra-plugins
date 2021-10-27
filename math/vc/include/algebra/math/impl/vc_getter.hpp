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

namespace algebra::vc::math {

/** This method retrieves phi from a vector, vector base with rows >= 2
 *
 * @param v the input vector
 **/
template <typename vector_type>
ALGEBRA_HOST_DEVICE inline auto phi(const vector_type &v) noexcept {

  return std::atan2(v[1], v[0]);
}

/** This method retrieves theta from a vector, vector base with rows >= 3
 *
 * @param v the input vector
 **/
template <typename vector_type>
ALGEBRA_HOST_DEVICE inline auto theta(const vector_type &v) noexcept {

  return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
}

/** This method retrieves the perpenticular magnitude of a vector with rows >= 2
 *
 * @param v the input vector
 **/
template <typename vector_type>
ALGEBRA_HOST_DEVICE inline auto perp(const vector_type &v) noexcept {

  return std::sqrt(v[0] * v[0] + v[1] * v[1]);
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
template <typename vector_type>
ALGEBRA_HOST_DEVICE inline auto norm(const vector_type &v) {

  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/** This method retrieves the pseudo-rapidity from a vector or vector base with
 *rows >= 3
 *
 * @param v the input vector
 **/
template <typename vector_type>
ALGEBRA_HOST_DEVICE inline auto eta(const vector_type &v) noexcept {

  return std::atanh(v[2] / norm(v));
}

}  // namespace algebra::vc::math
