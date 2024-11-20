/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/utils/algorithm_finder.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::cmath {

/// @returns a matrix of dimension @tparam ROW and @tparam COL that contains
/// the submatrix of @param m beginning at row @param row and column
/// @param col
template <std::size_t ROWS, std::size_t COLS, class input_matrix_type>
ALGEBRA_HOST_DEVICE decltype(auto) block(const input_matrix_type &m,
                                         std::size_t row, std::size_t col) {

  return algebra::cmath::block_getter().template operator()<ROWS, COLS>(m, row,
                                                                        col);
}

/// Sets a matrix of dimension @tparam ROW and @tparam COL as submatrix of
/// @param m beginning at row @param row and column @param col
template <std::size_t ROWS, std::size_t COLS, class input_matrix_type,
          typename scalar_t, template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE void set_block(
    input_matrix_type &m, const array_t<array_t<scalar_t, ROWS>, COLS> &b,
    std::size_t row, std::size_t col) {
  for (std::size_t j = 0u; j < COLS; ++j) {
    for (std::size_t i = 0u; i < ROWS; ++i) {
      m[j + col][i + row] = b[j][i];
    }
  }
}

/// Sets a vector of length @tparam ROW as submatrix of
/// @param m beginning at row @param row and column @param col
template <std::size_t ROWS, typename scalar_t,
          template <typename, std::size_t> class vector_t,
          class input_matrix_type>
ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                   const vector_t<scalar_t, ROWS> &b,
                                   std::size_t row, std::size_t col) {
  for (std::size_t i = 0; i < ROWS; ++i) {
    m[col][i + row] = b[i];
  }
}

/// @returns zero matrix of type @tparam matrix_t
template <typename matrix_t>
requires(std::is_scalar_v<typename matrix_t::value_type::value_type>)
    ALGEBRA_HOST_DEVICE inline matrix_t zero() {
  return matrix_t{};
}

/// Create zero matrix - cmath transform3
template <class matrix_t, typename element_getter_t>
ALGEBRA_HOST_DEVICE inline auto zero() {
  matrix_t ret;

  for (std::size_t j = 0; j < algebra::trait::columns<matrix_t>; ++j) {
    for (std::size_t i = 0; i < algebra::trait::rows<matrix_t>; ++i) {
      element_getter_t{}(ret, i, j) = 0;
    }
  }

  return ret;
}

/// @returns identity matrix of type @tparam matrix_t
template <typename matrix_t>
requires(std::is_scalar_v<typename matrix_t::value_type::value_type>)
    ALGEBRA_HOST_DEVICE inline auto identity() {
  auto ret{zero<matrix_t>()};

  for (std::size_t i = 0; i < algebra::trait::rank<matrix_t>; ++i) {
    ret[i][i] = 1;
  }

  return ret;
}

/// Create identity matrix - cmath transform3
template <class matrix_t, typename element_getter_t>
ALGEBRA_HOST_DEVICE inline matrix_t identity() {
  auto ret{zero<matrix_t, element_getter_t>()};

  for (std::size_t i = 0; i < algebra::trait::rank<matrix_t>; ++i) {
    element_getter_t{}(ret, i, i) = 1;
  }

  return ret;
}

/// Set @param m as zero matrix
template <std::size_t ROWS, std::size_t COLS, typename scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline void set_zero(
    array_t<array_t<scalar_t, ROWS>, COLS> &m) {
  m = zero<array_t<array_t<scalar_t, ROWS>, COLS>>();
}

/// Set @param m as identity matrix
template <std::size_t ROWS, std::size_t COLS, typename scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline void set_identity(
    array_t<array_t<scalar_t, ROWS>, COLS> &m) {
  m = identity<array_t<array_t<scalar_t, ROWS>, COLS>>();
}

/// @returns the transpose matrix of @param m
template <std::size_t ROWS, std::size_t COLS, typename scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline array_t<array_t<scalar_t, COLS>, ROWS> transpose(
    const array_t<array_t<scalar_t, ROWS>, COLS> &m) {

  array_t<array_t<scalar_t, COLS>, ROWS> ret;

  for (std::size_t i = 0; i < ROWS; ++i) {
    for (std::size_t j = 0; j < COLS; ++j) {
      ret[i][j] = m[j][i];
    }
  }

  return ret;
}

/// @returns the determinant of @param m
template <std::size_t N, typename scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline scalar_t determinant(
    const array_t<array_t<scalar_t, N>, N> &m) {

  using matrix_t = array_t<array_t<scalar_t, N>, N>;
  using element_getter_t = element_getter<std::size_t, array_t, scalar_t>;

  return determinant_t<N, matrix_t, element_getter_t>{}(m);
}

/// @returns the determinant of @param m
template <std::size_t N, typename scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline array_t<array_t<scalar_t, N>, N> inverse(
    const array_t<array_t<scalar_t, N>, N> &m) {

  using matrix_t = array_t<array_t<scalar_t, N>, N>;
  using element_getter_t = element_getter<std::size_t, array_t, scalar_t>;

  return inversion_t<N, matrix_t, element_getter_t>{}(m);
}

}  // namespace algebra::cmath
