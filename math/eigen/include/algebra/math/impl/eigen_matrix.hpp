/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace algebra::eigen::matrix {

/// "Matrix actor", assuming an Eigen matrix
template <typename scalar_t>
struct actor {

  /// 2D matrix type
  template <int ROWS, int COLS>
  using matrix_type = Eigen::Matrix<scalar_t, ROWS, COLS, 0, ROWS, COLS>;

  // Create zero matrix
  template <int ROWS, int COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> zero() {
    return matrix_type<ROWS, COLS>::Zero();
  }

  // Create identity matrix
  template <int ROWS, int COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> identity() {
    return matrix_type<ROWS, COLS>::Identity();
  }

  // Create transpose matrix
  template <int ROWS, int COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<COLS, ROWS> transpose(
      const matrix_type<ROWS, COLS>& m) {
    return m.transpose();
  }

  // Get determinant
  template <int N>
  ALGEBRA_HOST_DEVICE inline scalar_t determinant(const matrix_type<N, N>& m) {
    return m.determinant();
  }

  // Create inverse matrix
  template <int N>
  ALGEBRA_HOST_DEVICE inline matrix_type<N, N> inverse(
      const matrix_type<N, N>& m) {
    return m.inverse();
  }
};

}  // namespace algebra::eigen::matrix
