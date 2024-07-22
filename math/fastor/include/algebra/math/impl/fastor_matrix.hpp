/** Algebra plugins, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"
#include "algebra/storage/fastor.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

namespace algebra::fastor::math {

/// Create zero matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t zero() {
  return matrix_t(0);
}

/// Create identity matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t identity() {
  matrix_t identity_matrix;
  identity_matrix.eye();
  return identity_matrix;
}

/// Set input matrix as zero matrix
template <std::size_t ROWS, std::size_t COLS, typename scalar_t>
ALGEBRA_HOST_DEVICE inline void set_zero(
    Fastor::Tensor<scalar_t, ROWS, COLS> &m) {
  m.zeros();
}

/// Set input matrix as identity matrix
template <std::size_t ROWS, std::size_t COLS, typename scalar_t>
ALGEBRA_HOST_DEVICE inline void set_identity(
    Fastor::Tensor<scalar_t, ROWS, COLS> &m) {

  m = identity<ROWS, COLS>();
}

/// Create transpose matrix
template <std::size_t ROWS, std::size_t COLS, typename scalar_t>
ALGEBRA_HOST_DEVICE inline Fastor::Tensor<scalar_t, ROWS, COLS> transpose(
    const Fastor::Tensor<scalar_t, ROWS, COLS> &m) {
  return Fastor::transpose(m);
}

/// @returns the determinant of @param m
template <std::size_t ROWS, std::size_t COLS, typename scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t determinant(
    const Fastor::Tensor<scalar_t, ROWS, COLS> &m) {

  return Fastor::determinant(m);
}

/// @returns the inverse of @param m
template <std::size_t ROWS, std::size_t COLS, typename scalar_t>
ALGEBRA_HOST_DEVICE inline Fastor::Tensor<scalar_t, ROWS, COLS> inverse(
    const Fastor::Tensor<scalar_t, ROWS, COLS> &m) {

  return Fastor::inverse(m);
}

}  // namespace algebra::fastor::math
