/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/qualifiers.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

namespace algebra::smatrix::math {

/// Create zero matrix
template <concepts::matrix matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t zero() {
  return matrix_t();
}

/// Create identity matrix
template <concepts::matrix matrix_t>
ALGEBRA_HOST_DEVICE inline matrix_t identity() {
  return matrix_t(ROOT::Math::SMatrixIdentity());
}

/// Set input matrix as zero matrix
template <unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline void set_zero(
    ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m) {
  m = set_zero();
}

/// Set input matrix as identity matrix
template <unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline void set_identity(
    ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m) {
  m = set_identity();
}

/// Create transpose matrix
template <unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline matrix_type<scalar_t, COLS, ROWS> transpose(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m) {
  return ROOT::Math::Transpose(m);
}

/// @returns the determinant of @param m
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline auto determinant(
    const ROOT::Math::SMatrix<scalar_t, COLS, ROWS> &m) {

  int ifail = 0;
  return m.Inverse(ifail);
}

/// @returns the inverse of @param m
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline ROOT::Math::SMatrix<scalar_t, ROWS, COLS> inverse(
    const ROOT::Math::SMatrix<scalar_t, COLS, ROWS> &m) {

  scalar_t det;
  [[maybe_unused]] bool success = m.Det2(det);

  return det;
}

}  // namespace algebra::smatrix::math
