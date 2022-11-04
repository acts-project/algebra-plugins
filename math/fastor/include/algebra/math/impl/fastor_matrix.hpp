/** Algebra plugins, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"
#include "algebra/storage/fastor.hpp"

// Fastor include(s).
#include <Fastor/Fastor.h>

// System include(s).
#include <cstddef>  // for the std::size_t type

namespace algebra::fastor::math {

/// "Matrix actor", assuming an Fastor matrix
template <typename scalar_t>
struct actor {

  /// Size type
  using size_ty = std::size_t;

  /// Scalar type
  using scalar_type = scalar_t;

  /// 2D matrix type
  template <size_ty ROWS, size_ty COLS>
  using matrix_type = algebra::fastor::matrix_type<scalar_t, ROWS, COLS>;

  /// Array type
  template <size_ty N>
  using array_type = storage_type<scalar_type, N>;

  /// 3-element "vector" type
  using vector3 = array_type<3>;

  /// Operator getting a reference to one element of a non-const matrix
  template <size_ty ROWS, size_ty COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t &element(matrix_type<ROWS, COLS> &m,
                                               int row, int col) const {
    return m(row, col);
  }

  /// Operator getting one value of a const matrix
  template <size_ty ROWS, size_ty COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t element(const matrix_type<ROWS, COLS> &m,
                                              int row, int col) const {
    return m(row, col);
  }

  /// Operator getting a block of a const matrix
  template <size_ty ROWS, size_ty COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> block(const input_matrix_type &m,
                                                    size_ty row, size_ty col) {
    // In Fastor::fseq, the last element is not included.
    return m(Fastor::seq(row, row + ROWS), Fastor::seq(col, col + COLS));
  }

  /// Operator setting a block
  template <size_ty ROWS, size_ty COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                     const matrix_type<ROWS, COLS> &b,
                                     size_ty row, size_ty col) {
    m(Fastor::fseq<row, row + ROWS>(), Fastor::fseq<col, col + COLS>()) = b;
  }

  // Create zero matrix
  template <size_ty ROWS, size_ty COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> zero() {
    return matrix_type<ROWS, COLS>().zeros();
  }

  // Create identity matrix
  template <size_ty ROWS, size_ty COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<ROWS, COLS> identity() {
    // There are 2 identity tensor methods in Fastor, eye() and eye2(). The
    // former is for arbitrary order tensors, whereas the latter is specifically
    // for second order tensors. As such, I chose to use eye2() here because it
    // does less and hence would be faster.
    return matrix_type<ROWS, COLS>().eye2();
  }

  // Set input matrix as zero matrix
  template <size_ty ROWS, size_ty COLS>
  ALGEBRA_HOST_DEVICE inline void set_zero(matrix_type<ROWS, COLS> &m) const {
    m.zeros();
  }

  // Set input matrix as identity matrix
  template <size_ty ROWS, size_ty COLS>
  ALGEBRA_HOST_DEVICE inline void set_identity(
      matrix_type<ROWS, COLS> &m) const {
    m.eye2();
  }

  // Create transpose matrix
  template <size_ty ROWS, size_ty COLS>
  ALGEBRA_HOST_DEVICE inline matrix_type<COLS, ROWS> transpose(
      const matrix_type<ROWS, COLS> &m) {
    return Fastor::transpose(m);
  }

  // Get determinant
  template <size_ty N>
  ALGEBRA_HOST_DEVICE inline scalar_t determinant(const matrix_type<N, N> &m) {
    return Fastor::determinant(m);
  }

  // Create inverse matrix
  template <size_ty N>
  ALGEBRA_HOST_DEVICE inline matrix_type<N, N> inverse(
      const matrix_type<N, N> &m) {
    return Fastor::inverse(m);
  }
};

}  // namespace algebra::fastor::math
