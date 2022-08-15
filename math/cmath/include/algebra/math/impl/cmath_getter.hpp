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

// System include(s).
#include <cassert>
#include <cstddef>
#include <type_traits>

namespace algebra::cmath {

/// "Element getter", assuming a simple 2D array access
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t>
struct element_getter {

  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = array_t<array_t<scalar_t, ROWS>, COLS>;

  /// Operator getting a reference to one element of a non-const matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t &operator()(matrix_type<ROWS, COLS> &m,
                                                  std::size_t row,
                                                  std::size_t col) const {

    assert(row < ROWS);
    assert(col < COLS);
    return m[col][row];
  }

  /// Operator getting one value of a const matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const matrix_type<ROWS, COLS> &m, std::size_t row,
      std::size_t col) const {

    assert(row < ROWS);
    assert(col < COLS);
    return m[col][row];
  }
};  // struct element_getter

/// Function extracting an element from a matrix (const)
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t, size_type ROWS, size_type COLS>
ALGEBRA_HOST_DEVICE inline scalar_t element(
    const array_t<array_t<scalar_t, ROWS>, COLS> &m, std::size_t row,
    std::size_t col) {

  return element_getter<size_type, array_t, scalar_t>()(m, row, col);
}

/// Function extracting an element from a matrix (non-const)
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t, size_type ROWS, size_type COLS>
ALGEBRA_HOST_DEVICE inline scalar_t &element(
    array_t<array_t<scalar_t, ROWS>, COLS> &m, std::size_t row,
    std::size_t col) {

  return element_getter<size_type, array_t, scalar_t>()(m, row, col);
}

/// "Vector getter", assuming a simple 2D array access
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t, size_type SIZE,
          typename result_t = array_t<scalar_t, SIZE>>
struct vector_getter {

  /// Result type
  using result_type = result_t;
  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = array_t<array_t<scalar_t, ROWS>, COLS>;

  /// Operator producing a vector out of a const matrix
  template <size_type ROWS, size_type COLS>
  ALGEBRA_HOST_DEVICE inline result_type operator()(
      const matrix_type<ROWS, COLS> &m, std::size_t row, std::size_t col) {

    assert(col < COLS);
    assert(row + SIZE < ROWS);
    result_type subvector{};
    for (std::size_t irow = row; irow < row + SIZE; ++irow) {
      subvector[irow - row] = m[col][irow];
    }
    return subvector;
  }
};  // struct vector_getter

/// "Block getter", assuming a simple 2D array access
template <typename size_type, template <typename, size_type> class array_t,
          typename scalar_t>
struct block_getter {

  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = array_t<array_t<scalar_t, ROWS>, COLS>;

  /// Operator producing a sub-matrix from a const matrix
  template <size_type ROWS, size_type COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> operator()(
      const input_matrix_type &m, std::size_t row, std::size_t col) const {

    matrix_type<ROWS, COLS> submatrix{};
    for (std::size_t icol = col; icol < col + COLS; ++icol) {
      for (std::size_t irow = row; irow < row + ROWS; ++irow) {
        submatrix[icol - col][irow - row] = m[icol][irow];
      }
    }
    return submatrix;
  }
};  // struct block_getter

}  // namespace algebra::cmath
