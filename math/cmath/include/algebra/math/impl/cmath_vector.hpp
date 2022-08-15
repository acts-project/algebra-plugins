/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::cmath::vector {

template <typename size_t, template <typename, size_t> class array_t,
          typename scalar_t>
struct actor {

  using size_type = size_t;
  using scalar_type = scalar_t;
  template <typename arg_t, size_type N>
  using array_type = array_t<arg_t, N>;

  /** This method retrieves phi from a vector, vector base with rows >= 2
   *
   * @param v the input vector
   **/
  template <size_type N, std::enable_if_t<N >= 2, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type phi(
      const array_type<scalar_type, N> &v) noexcept {

    return algebra::math::atan2(v[1], v[0]);
  }

  /** This method retrieves theta from a vector, vector base with rows >= 3
   *
   * @param v the input vector
   **/
  template <size_type N, std::enable_if_t<N >= 3, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type theta(
      const array_type<scalar_type, N> &v) noexcept {

    return algebra::math::atan2(algebra::math::sqrt(v[0] * v[0] + v[1] * v[1]),
                                v[2]);
  }

  /** This method retrieves the perpenticular magnitude of a vector with rows >=
   *2
   *
   * @param v the input vector
   **/
  template <size_type N, std::enable_if_t<N >= 2, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type perp(
      const array_type<scalar_type, N> &v) noexcept {

    return algebra::math::sqrt(v[0] * v[0] + v[1] * v[1]);
  }

  /** This method retrieves the norm of a vector, no dimension restriction
   *
   * @param v the input vector
   **/
  template <size_type N, std::enable_if_t<N == 2, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type norm(
      const array_type<scalar_type, N> &v) {

    return perp(v);
  }

  /** This method retrieves the norm of a vector, no dimension restriction
   *
   * @param v the input vector
   **/
  template <size_type N, std::enable_if_t<N >= 3, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type norm(
      const array_type<scalar_type, N> &v) {

    return algebra::math::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }

  /** This method retrieves the pseudo-rapidity from a vector or vector base
   *with rows >= 3
   *
   * @param v the input vector
   **/
  template <size_type N, std::enable_if_t<N >= 3, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type eta(
      const array_type<scalar_type, N> &v) noexcept {

    return algebra::math::atanh(v[2] / norm(v));
  }

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
  template <
      size_type N,
      std::enable_if_t<N >= 3 && std::is_scalar_v<scalar_type>, bool> = true>
  ALGEBRA_HOST_DEVICE inline array_type<scalar_type, N> cross(
      const array_type<scalar_type, N> &a,
      const array_type<scalar_type, N> &b) {

    return {a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0],
            a[0] * b[1] - b[0] * a[1]};
  }

  /** Cross product between two input vectors - 3 Dim
   *
   * @tparam derived_type_lhs is the first matrix (epresseion) template
   * @tparam derived_type_rhs is the second matrix (epresseion) template
   *
   * @param a the first input vector
   * @param b the second input matrix with single column
   *
   * @return a vector (expression) representing the cross product
   **/
  template <
      size_type N, size_type COLS,
      std::enable_if_t<N >= 3 && COLS == 1 && std::is_scalar_v<scalar_type>,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline array_type<scalar_type, N> cross(
      const array_type<scalar_type, N> &a,
      const array_type<array_type<scalar_type, N>, COLS> &b) {

    return {a[1] * b[0][2] - b[0][1] * a[2], a[2] * b[0][0] - b[0][2] * a[0],
            a[0] * b[0][1] - b[0][0] * a[1]};
  }

  /** Cross product between two input vectors - 3 Dim
   *
   * @tparam derived_type_lhs is the first matrix (epresseion) template
   * @tparam derived_type_rhs is the second matrix (epresseion) template
   *
   * @param a the first input matrix with single column
   * @param b the second input vector
   *
   * @return a vector (expression) representing the cross product
   **/
  template <
      size_type N, size_type COLS,
      std::enable_if_t<N >= 3 && COLS == 1 && std::is_scalar_v<scalar_type>,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline array_type<scalar_type, N> cross(
      const array_type<array_type<scalar_type, N>, COLS> &a,
      const array_type<scalar_type, N> &b) {

    return {a[0][1] * b[2] - b[1] * a[0][2], a[0][2] * b[0] - b[2] * a[0][0],
            a[0][0] * b[1] - b[0] * a[0][1]};
  }

  /** Cross product between two input vectors - 3 Dim
   *
   * @tparam derived_type_lhs is the first matrix (epresseion) template
   * @tparam derived_type_rhs is the second matrix (epresseion) template
   *
   * @param a the first input matrix with single column
   * @param b the second input matrix with single column
   *
   * @return a vector (expression) representing the cross product
   **/
  template <
      size_type N, size_type COLS,
      std::enable_if_t<N >= 3 && COLS == 1 && std::is_scalar_v<scalar_type>,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline array_type<scalar_type, N> cross(
      const array_type<array_type<scalar_type, N>, COLS> &a,
      const array_type<array_type<scalar_type, N>, COLS> &b) {

    return {a[0][1] * b[0][2] - b[0][1] * a[0][2],
            a[0][2] * b[0][0] - b[0][2] * a[0][0],
            a[0][0] * b[0][1] - b[0][0] * a[0][1]};
  }

  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <size_type N,
            std::enable_if_t<std::is_scalar_v<scalar_type>, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type dot(
      const array_type<scalar_type, N> &a,
      const array_type<scalar_type, N> &b) {
    scalar_type ret = 0;
    for (size_type i = 0; i < N; i++) {
      ret += a[i] * b[i];
    }
    return ret;
  }

  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input matrix with single column
   *
   * @return the scalar dot product value
   **/
  template <
      size_type N, size_type COLS,
      std::enable_if_t<COLS == 1 && std::is_scalar_v<scalar_type>, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type dot(
      const array_type<scalar_type, N> &a,
      const array_type<array_type<scalar_type, N>, COLS> &b) {
    scalar_type ret = 0;
    for (size_type i = 0; i < N; i++) {
      ret += a[i] * b[0][i];
    }
    return ret;
  }

  /** Dot product between two input vectors
   *
   * @param a the first input matrix with single column
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <
      size_type N, size_type COLS,
      std::enable_if_t<COLS == 1 && std::is_scalar_v<scalar_type>, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_type dot(
      const array_type<array_type<scalar_type, N>, COLS> &a,
      const array_type<scalar_type, N> &b) {
    scalar_type ret = 0;
    for (size_type i = 0; i < N; i++) {
      ret += a[0][i] * b[i];
    }
    return ret;
  }

  /** Dot product between two input vectors
   *
   * @param a the first input matrix with single column
   * @param b the second input matrix with single column
   *
   * @return the scalar dot product value
   **/
  template <
      size_type N, size_type COLS,
      std::enable_if_t<COLS == 1 && std::is_scalar_v<scalar_t>, bool> = true>
  ALGEBRA_HOST_DEVICE inline scalar_t dot(
      const array_t<array_t<scalar_t, N>, COLS> &a,
      const array_t<array_t<scalar_t, N>, COLS> &b) {
    scalar_t ret = 0;
    for (size_type i = 0; i < N; i++) {
      ret += a[0][i] * b[0][i];
    }
    return ret;
  }

  /** Get a normalized version of the input vector
   *
   * @param v the input vector
   **/
  template <size_type N,
            std::enable_if_t<std::is_scalar_v<scalar_type>, bool> = true>
  ALGEBRA_HOST_DEVICE inline array_type<scalar_type, N> normalize(
      const array_type<scalar_type, N> &v) {

    const scalar_type oon =
        static_cast<scalar_type>(1.) / algebra::math::sqrt(dot(v, v));
    return v * oon;
  }
};

}  // namespace algebra::cmath::vector
