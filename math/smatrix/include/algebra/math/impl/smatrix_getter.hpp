/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"
#include "algebra/storage/smatrix.hpp"

// ROOT/Smatrix include(s).
#include <Math/Functions.h>

// System include(s).
#include <cmath>

namespace algebra::smatrix::math {

/** This method retrieves phi from a vector, vector base with rows >= 2
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N, std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST inline scalar_t phi(const storage_type<scalar_t, N> &v) noexcept {

  scalar_t element0 = v.apply(0);
  scalar_t element1 = v.apply(1);
  return std::atan2(element1, element0);
}

/** This method retrieves theta from a vector, vector base with rows >= 3
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N, std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST inline scalar_t theta(
    const storage_type<scalar_t, N> &v) noexcept {

  return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST inline scalar_t norm(const storage_type<scalar_t, N> &v) {

  return std::sqrt(ROOT::Math::Dot(v, v));
}

/** This method retrieves the pseudo-rapidity from a vector or vector base with
 *rows >= 3
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N, std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST inline scalar_t eta(const storage_type<scalar_t, N> &v) noexcept {

  return std::atanh(v[2] / norm(v));
}

/** This method retrieves the perpenticular magnitude of a vector with rows >= 2
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N, std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST inline scalar_t perp(const storage_type<scalar_t, N> &v) noexcept {

  scalar element0 = v.apply(0);
  scalar element1 = v.apply(1);
  return std::sqrt(element0 * element0 + element1 * element1);
}

}  // namespace algebra::smatrix::math
