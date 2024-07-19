/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

namespace algebra::smatrix::math {

/// Operator getting a block of a const matrix
template <unsigned int ROWS, unsigned int COLS, class input_matrix_type>
ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> block(const input_matrix_type &m,
                                                  unsigned int row,
                                                  unsigned int col) {
  return m.template Sub<matrix_type<ROWS, COLS> >(row, col);
}

/// Operator setting a block with a matrix
template <unsigned int ROWS, unsigned int COLS, class input_matrix_type>
ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                   const matrix_type<ROWS, COLS> &b,
                                   unsigned int row, unsigned int col) {
  for (unsigned int i = 0; i < ROWS; ++i) {
    for (unsigned int j = 0; j < COLS; ++j) {
      m(i + row, j + col) = b(i, j);
    }
  }
}

/// Operator setting a block with a vector
template <unsigned int ROWS, class input_matrix_type>
ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                   const vector_type<ROWS> &b, unsigned int row,
                                   unsigned int col) {
  for (unsigned int i = 0; i < ROWS; ++i) {
    m(i + row, col) = b[i];
  }
}

// Create zero matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t zero() {
  return matrix_t();
}

// Create identity matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t identity() {
  return matrix_t(ROOT::Math::SMatrixIdentity());
}

// Set input matrix as zero matrix
template <unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline void set_zero(matrix_type<ROWS, COLS> &m) {

  for (unsigned int i = 0; i < ROWS; ++i) {
    for (unsigned int j = 0; j < COLS; ++j) {
      m(i, j) = 0;
    }
  }
}

// Set input matrix as identity matrix
template <unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline void set_identity(matrix_type<ROWS, COLS> &m) {

  for (unsigned int i = 0; i < ROWS; ++i) {
    for (unsigned int j = 0; j < COLS; ++j) {
      if (i == j) {
        m(i, j) = 1;
      } else {
        m(i, j) = 0;
      }
    }
  }
}

// Create transpose matrix
template <unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline matrix_type<COLS, ROWS> transpose(
    const matrix_type<ROWS, COLS> &m) {
  return ROOT::Math::Transpose(m);
}

/// @returns the determinant of @param m
template <std::size_t ROWS, std::size_t COLS, typename scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline auto determinant(const matrix_type<COLS, ROWS> &m) {

  int ifail = 0;
  return m.Inverse(ifail);
}

/// @returns the inverse of @param m
template <std::size_t ROWS, std::size_t COLS, typename scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> inverse(
    const matrix_type<COLS, ROWS> &m) {

  scalar_t det;
  [[maybe_unused]] bool success = m.Det2(det);

  return det;
}

}  // namespace algebra::smatrix::math
