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

namespace algebra::smatrix::matrix {

/// "Matrix actor", assuming an Eigen matrix
template <typename scalar_t>
struct actor {

  /// 2D matrix type
  template <unsigned int ROWS, unsigned int COLS>
  using matrix_type = ROOT::Math::SMatrix<scalar_t, ROWS, COLS>;

  /// Operator getting a reference to one element of a non-const matrix
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t &element(matrix_type<ROWS, COLS> &m,
                                               unsigned int row,
                                               unsigned int col) const {
    return m(row, col);
  }

  /// Operator getting one value of a const matrix
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t element(const matrix_type<ROWS, COLS> &m,
                                              unsigned int row,
                                              unsigned int col) const {
    return m(row, col);
  }

  /// Operator getting a block of a const matrix
  template <unsigned int ROWS, unsigned int COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> block(const input_matrix_type &m,
                                                    unsigned int row,
                                                    unsigned int col) {
    return m.template Sub<matrix_type<ROWS, COLS> >(row, col);
  }

  // Create zero matrix
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> zero() {
    return matrix_type<ROWS, COLS>();
  }

  // Create identity matrix
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> identity() {
    return matrix_type<ROWS, COLS>(ROOT::Math::SMatrixIdentity());
  }

  // Set input matrix as zero matrix
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline void set_zero(matrix_type<ROWS, COLS> &m) const {

    for (unsigned int i = 0; i < ROWS; ++i) {
      for (unsigned int j = 0; j < COLS; ++j) {
        m(i, j) = 0;
      }
    }
  }

  // Set input matrix as identity matrix
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline void set_identity(
      matrix_type<ROWS, COLS> &m) const {

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

  // Get determinant
  template <unsigned int N>
  ALGEBRA_HOST_DEVICE inline scalar_t determinant(const matrix_type<N, N> &m) {
    scalar_t det;
    bool success = m.Det2(det);

    // suppress unused parameter warning
    (void)success;

    return det;
  }

  // Create inverse matrix
  template <unsigned int N>
  ALGEBRA_HOST_DEVICE inline matrix_type<N, N> inverse(
      const matrix_type<N, N> &m) {
    int ifail = 0;
    return m.Inverse(ifail);
  }
};

}  // namespace algebra::smatrix::matrix
