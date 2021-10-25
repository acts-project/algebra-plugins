/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// System include(s).
#include <cmath>
#include <cstddef>
#include <type_traits>

namespace algebra::cmath {

/** This method retrieves phi from a vector, vector base with rows >= 2
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t, auto N,
          std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t phi(
    const array_t<scalar_t, N> &v) noexcept {

  return std::atan2(v[1], v[0]);
}

/** This method retrieves theta from a vector, vector base with rows >= 3
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t, auto N,
          std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t theta(
    const array_t<scalar_t, N> &v) noexcept {

  return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
}

/** This method retrieves the perpenticular magnitude of a vector with rows >= 2
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t, auto N,
          std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t perp(
    const array_t<scalar_t, N> &v) noexcept {

  return std::sqrt(v[0] * v[0] + v[1] * v[1]);
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t norm(const array_t<scalar_t, 2> &v) {

  return perp<array_t>(v);
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t, auto N,
          std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t norm(const array_t<scalar_t, N> &v) {

  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/** This method retrieves the pseudo-rapidity from a vector or vector base with
 *rows >= 3
 *
 * @param v the input vector
 **/
template <template <typename, auto> class array_t, typename scalar_t, auto N,
          std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST_DEVICE inline scalar_t eta(
    const array_t<scalar_t, N> &v) noexcept {

  return std::atanh(v[2] / norm<array_t>(v));
}

/// "Element getter", assuming a simple 2D array access
template <template <typename, auto> class array_t, typename scalar_t>
struct element_getter {

  /// 2D matrix type
  template <auto ROWS, auto COLS>
  using matrix_type = array_t<array_t<scalar_t, ROWS>, COLS>;

  /// Operator getting a reference to one element of a non-const matrix
  template <auto ROWS, auto COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t &operator()(matrix_type<ROWS, COLS> &m,
                                                  std::size_t row,
                                                  std::size_t col) const {

    return m[col][row];
  }

  /// Operator getting one value of a const matrix
  template <auto ROWS, auto COLS>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const matrix_type<ROWS, COLS> &m, std::size_t row,
      std::size_t col) const {

    return m[col][row];
  }
};  // struct element_getter

/// "Vector getter", assuming a simple 2D array access
template <template <typename, auto> class array_t, typename scalar_t>
struct vector_getter {

  /// 1D vector type
  template <auto SIZE>
  using vector_type = array_t<scalar_t, SIZE>;
  /// 2D matrix type
  template <auto ROWS, auto COLS>
  using matrix_type = array_t<array_t<scalar_t, ROWS>, COLS>;

  /// Operator producing a vector out of a const matrix
  template <auto SIZE, auto ROWS, auto COLS>
  ALGEBRA_HOST_DEVICE inline vector_type<SIZE> operator()(
      const matrix_type<ROWS, COLS> &m, std::size_t row, std::size_t col) {

    vector_type<SIZE> subvector;
    for (std::size_t irow = row; irow < row + SIZE; ++irow) {
      subvector[irow - row] = m[col][irow];
    }
    return subvector;
  }
};  // struct vector_getter

/// "Block getter", assuming a simple 2D array access
template <template <typename, auto> class array_t, typename scalar_t>
struct block_getter {

  /// 2D matrix type
  template <auto ROWS, auto COLS>
  using matrix_type = array_t<array_t<scalar_t, ROWS>, COLS>;

  /// Operator producing a sub-matrix from a const matrix
  template <auto ROWS, auto COLS, class input_matrix_type>
  ALGEBRA_HOST_DEVICE matrix_type<ROWS, COLS> operator()(
      const input_matrix_type &m, std::size_t row, std::size_t col) const {

    matrix_type<ROWS, COLS> submatrix;
    for (std::size_t icol = col; icol < col + COLS; ++icol) {
      for (std::size_t irow = row; irow < row + ROWS; ++irow) {
        submatrix[icol - col][irow - row] = m[icol][irow];
      }
    }
    return submatrix;
  }
};  // struct block_getter

}  // namespace algebra::cmath
