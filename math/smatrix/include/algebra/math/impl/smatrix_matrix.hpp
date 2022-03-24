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

  // Create zero matrix
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> zero() {
    return matrix_type<ROWS, COLS>();
  }

  // Create identity matrix
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> identity() {
    matrix_type<ROWS, COLS> ret = ROOT::Math::SMatrixIdentity();
    return ret;
  }

  // Create transpose matrix
  template <unsigned int ROWS, unsigned int COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<COLS, ROWS> transpose(
      const matrix_type<ROWS, COLS>& m) {
    matrix_type<COLS, ROWS> ret;

    for (unsigned int i = 0; i < ROWS; ++i) {
      for (unsigned int j = 0; j < COLS; ++j) {
        ret(j, i) = m(i, j);
      }
    }

    return ret;
  }

  // Get determinant
  template <unsigned int N>
  ALGEBRA_HOST_DEVICE inline scalar_t determinant(const matrix_type<N, N>& m) {
    scalar_t det;
    bool success = m.Det2(det);

    // suppress unused parameter warning
    (void)success;

    return det;
  }

  // Create inverse matrix
  template <unsigned int N>
  ALGEBRA_HOST_DEVICE inline matrix_type<N, N> inverse(
      const matrix_type<N, N>& m) {
    int ifail = 0;
    return m.Inverse(ifail);
  }
};

}  // namespace algebra::smatrix::matrix
