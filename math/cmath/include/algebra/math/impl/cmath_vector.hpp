/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/common.hpp"
#include "algebra/math/generic.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::cmath {

/// This method retrieves phi from a vector with rows >= 2
///
/// @param v the input vector
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N,
          std::enable_if_t<(N >= 2), bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t phi(const array_t<scalar_t, N> &v) {
  return algebra::generic::math::phi(v);
}

/// This method retrieves the perpendicular magnitude of a vector with rows >= 2
///
/// @param v the input vector
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N,
          std::enable_if_t<(N >= 2), bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t perp(const array_t<scalar_t, N> &v) {
  return algebra::generic::math::perp(v);
}

/// This method retrieves theta from a vector with rows >= 3
///
/// @param v the input vector
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N,
          std::enable_if_t<(N >= 2), bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t theta(const array_t<scalar_t, N> &v) {
  return algebra::generic::math::theta(v);
}

/// Cross product between two input vectors - 3 Dim
///
/// @tparam vector1_t first vector or column matrix type
/// @tparam vector2_t second vector or column matrix type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector representing the cross product
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N,
          std::enable_if_t<(N >= 3), bool> = true>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, N> cross(
    const array_t<scalar_t, N> &a, const array_t<scalar_t, N> &b) {
  return algebra::generic::math::cross(a, b);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N>
ALGEBRA_HOST_DEVICE inline scalar_t dot(const array_t<scalar_t, N> &a,
                                        const array_t<scalar_t, N> &b) {
  array_t<scalar_t, N> tmp;
  for (size_type i = 0; i < N; i++) {
    tmp[i] = a[i] * b[i];
  }
  scalar_t ret{0.f};
  for (size_type i = 0; i < N; i++) {
    ret += tmp[i];
  }
  return ret;
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input matrix with single column
///
/// @return the scalar dot product value
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N, size_type COLS,
          std::enable_if_t<COLS == 1, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t dot(
    const array_t<scalar_t, N> &a,
    const array_t<array_t<scalar_t, N>, COLS> &b) {
  array_t<scalar_t, N> tmp;
  for (size_type i = 0; i < N; i++) {
    tmp[i] = a[i] * b[0][i];
  }
  scalar_t ret{0.f};
  for (size_type i = 0; i < N; i++) {
    ret += tmp[i];
  }
  return ret;
}

/// Dot product between two input vectors
///
/// @param a the first input matrix with single column
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N, size_type COLS,
          std::enable_if_t<COLS == 1, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t dot(
    const array_t<array_t<scalar_t, N>, COLS> &a,
    const array_t<scalar_t, N> &b) {
  array_t<scalar_t, N> tmp;
  for (size_type i = 0; i < N; i++) {
    tmp[i] = a[0][i] * b[i];
  }
  scalar_t ret{0.f};
  for (size_type i = 0; i < N; i++) {
    ret += tmp[i];
  }
  return ret;
}

/// Dot product between two input vectors
///
/// @param a the first input matrix with single column
/// @param b the second input matrix with single column
///
/// @return the scalar dot product value
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N, size_type COLS,
          std::enable_if_t<COLS == 1, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t dot(
    const array_t<array_t<scalar_t, N>, COLS> &a,
    const array_t<array_t<scalar_t, N>, COLS> &b) {
  array_t<scalar_t, N> tmp;
  for (size_type i = 0; i < N; i++) {
    tmp[i] = a[0][i] * b[0][i];
  }
  scalar_t ret{0.f};
  for (size_type i = 0; i < N; i++) {
    ret += tmp[i];
  }
  return ret;
}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @param v the input vector
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N,
          std::enable_if_t<(N >= 2), bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t norm(const array_t<scalar_t, N> &v) {

  return algebra::math::sqrt(dot(v, v));
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N,
          std::enable_if_t<(N >= 3), bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t eta(
    const array_t<scalar_t, N> &v) noexcept {

  return algebra::math::atanh(v[2] / norm(v));
}

/// Get a normalized version of the input vector
///
/// @param v the input vector
template <concepts::index size_type,
          template <typename, size_type> class array_t,
          concepts::scalar scalar_t, size_type N>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, N> normalize(
    const array_t<scalar_t, N> &v) {

  return (static_cast<scalar_t>(1.) / norm(v)) * v;
}

}  // namespace algebra::cmath
