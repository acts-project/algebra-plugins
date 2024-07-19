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

// System include(s).
#include <cstddef>  // for the std::size_t type

namespace algebra::fastor::math {

/// Operator getting a block of a const matrix
template <std::size_t ROWS, std::size_t COLS, class input_matrix_type>
ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> block(const input_matrix_type &m,
                                                  std::size_t row,
                                                  std::size_t col) {
  // In `Fastor::seq`, the last element is not included.
  // `Fastor::seq` takes `int`s as input, but `row`, `col`, `ROWS`, and `COLS`
  // have type `std::size_t`.
  return m(Fastor::seq(static_cast<int>(row), static_cast<int>(row + ROWS)),
           Fastor::seq(static_cast<int>(col), static_cast<int>(col + COLS)));
}

/// Operator setting a block with a matrix
template <std::size_t ROWS, std::size_t COLS, class input_matrix_type>
ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                   const matrix_type<ROWS, COLS> &b,
                                   std::size_t row, std::size_t col) {
  // In `Fastor::seq`, the last element is not included.
  // `Fastor::seq` takes `int`s as input, but `ROWS` and `COLS` have type
  // `std::size_t`, which is `std::size_t`.
  m(Fastor::seq(static_cast<int>(row), static_cast<int>(row + ROWS)),
    Fastor::seq(static_cast<int>(col), static_cast<int>(col + COLS))) = b;
}

/// Operator setting a block with a vector
template <std::size_t ROWS, class input_matrix_type>
ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                   const Fastor::Tensor<scalar_t, ROWS> &b,
                                   std::size_t row, std::size_t col) {
  // In `Fastor::seq`, the last element is not included.
  // `Fastor::seq` takes `int`s as input, but `ROWS` and `COLS` have type
  // `std::size_t`, which is `std::size_t`.
  m(Fastor::seq(static_cast<int>(row), static_cast<int>(row + ROWS)),
    static_cast<int>(col)) = b;
}

// Create zero matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t zero() {
  return matrix_t(0);
}

// Create identity matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t identity() {
  matrix_t identity_matrix;
  identity_matrix.eye();
  return identity_matrix;
}

// Set input matrix as zero matrix
template <std::size_t ROWS, std::size_t COLS, typename scalar_t>
ALGEBRA_HOST_DEVICE inline void set_zero(
    Fastor::Tensor<scalar_t, ROWS, COLS> &m) {
  m.zeros();
}

// Set input matrix as identity matrix
template <std::size_t ROWS, std::size_t COLS, typename scalar_t>
ALGEBRA_HOST_DEVICE inline void set_identity(
    Fastor::Tensor<scalar_t, ROWS, COLS> &m) {

  m = identity<ROWS, COLS>();
}

// Create transpose matrix
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
