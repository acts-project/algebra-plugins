/** Algebra plugins, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/storage/matrix.hpp"
#include "algebra/storage/vector.hpp"

// System include(s).
#include <cassert>

namespace algebra::storage {

/// Functor used to access elements of a matrix
struct element_getter {

  /// Get const access to a matrix element
  template <template <typename, std::size_t> class array_t,
            concepts::scalar scalar_t, std::size_t ROW, std::size_t COL>
  ALGEBRA_HOST_DEVICE inline decltype(auto) operator()(
      const matrix<array_t, scalar_t, ROW, COL> &m, std::size_t row,
      std::size_t col) const {

    // Make sure that the indices are valid.
    assert(row < ROW);
    assert(col < COL);

    // Return the selected element.
    return m[col][row];
  }

  /// Get non-const access to a matrix element
  template <template <typename, std::size_t> class array_t,
            concepts::scalar scalar_t, std::size_t ROW, std::size_t COL>
  ALGEBRA_HOST_DEVICE inline decltype(auto) operator()(
      matrix<array_t, scalar_t, ROW, COL> &m, std::size_t row,
      std::size_t col) const {
    assert(row < ROW);
    assert(col < COL);

    return m[col][row];
  }

  /// Get const access to a matrix element
  template <template <typename, std::size_t> class array_t,
            concepts::scalar scalar_t, std::size_t ROW>
  ALGEBRA_HOST_DEVICE inline decltype(auto) operator()(
      const matrix<array_t, scalar_t, ROW, 1> &m, std::size_t row) const {

    assert(row < ROW);
    return m[0][row];
  }

  /// Get non-const access to a matrix element
  template <template <typename, std::size_t> class array_t,
            concepts::scalar scalar_t, std::size_t ROW>
  ALGEBRA_HOST_DEVICE inline decltype(auto) operator()(
      matrix<array_t, scalar_t, ROW, 1> &m, std::size_t row) const {

    assert(row < ROW);
    return m[0][row];
  }

  /// Get const access to a vector element
  template <template <typename, std::size_t> class array_t,
            concepts::scalar scalar_t, std::size_t N>
  ALGEBRA_HOST_DEVICE inline decltype(auto) operator()(
      const vector<N, scalar_t, array_t> &v, std::size_t row) const {

    assert(row < N);
    return v[row];
  }

  /// Get non-const access to a vector element
  template <template <typename, std::size_t> class array_t,
            concepts::scalar scalar_t, std::size_t N>
  ALGEBRA_HOST_DEVICE inline decltype(auto) operator()(
      vector<N, scalar_t, array_t> &v, std::size_t row) const {

    assert(row < N);
    return v[row];
  }

};  // struct element_getter

/// Function extracting an element from a matrix (const)
template <std::size_t ROW, std::size_t COL, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline decltype(auto) element(
    const matrix<array_t, scalar_t, ROW, COL> &m, std::size_t row,
    std::size_t col) {
  return element_getter{}(m, row, col);
}

/// Function extracting an element from a matrix (non-const)
template <std::size_t ROW, std::size_t COL, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline decltype(auto) element(
    matrix<array_t, scalar_t, ROW, COL> &m, std::size_t row, std::size_t col) {
  return element_getter{}(m, row, col);
}

/// Function extracting an element from a 1D matrix (const)
template <std::size_t ROW, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline decltype(auto) element(
    const matrix<array_t, scalar_t, ROW, 1> &m, std::size_t row) {
  return element_getter{}(m, row);
}

/// Function extracting an element from a 1D matrix (non-const)
template <std::size_t ROW, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline decltype(auto) element(
    matrix<array_t, scalar_t, ROW, 1> &m, std::size_t row) {
  return element_getter{}(m, row);
}

/// Function extracting an element from a vector (const)
template <std::size_t N, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline decltype(auto) element(
    const vector<N, scalar_t, array_t> &v, std::size_t row) {
  return element_getter{}(v, row);
}

/// Function extracting an element from a vector (non-const)
template <std::size_t N, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline decltype(auto) element(
    vector<N, scalar_t, array_t> &v, std::size_t row) {
  return element_getter{}(v, row);
}

/// Functor used to access a submatrix of a matrix
struct block_getter {

  /// Get a block of a const matrix
  template <std::size_t ROWS, std::size_t COLS, std::size_t mROW,
            std::size_t mCOL, concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  ALGEBRA_HOST_DEVICE constexpr auto operator()(
      const matrix<array_t, scalar_t, mROW, mCOL> &m, const std::size_t row,
      const std::size_t col) noexcept {

    static_assert(ROWS <= mROW);
    static_assert(COLS <= mCOL);
    assert(row + ROWS <= mROW);
    assert(col + COLS <= mCOL);

    using input_matrix_t = matrix<array_t, scalar_t, mROW, mCOL>;
    using matrix_t = matrix<array_t, scalar_t, ROWS, COLS>;

    matrix_t res_m;

    // Don't access single elements in underlying vectors unless necessary
    if constexpr (matrix_t::storage_rows() == input_matrix_t::storage_rows()) {
      if (row == 0u) {
        for (std::size_t j = col; j < col + COLS; ++j) {
          res_m[j - col] = m[j];
        }

        return res_m;
      }
    }

    for (std::size_t j = col; j < col + COLS; ++j) {
      for (std::size_t i = row; i < row + ROWS; ++i) {
        res_m[j - col][i - row] = m[j][i];
      }
    }

    return res_m;
  }

  /// Get a vector of a const matrix
  template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS,
            concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  ALGEBRA_HOST_DEVICE constexpr auto operator()(
      const matrix<array_t, scalar_t, ROWS, COLS> &m, const std::size_t row,
      const std::size_t col) noexcept {

    static_assert(SIZE <= ROWS);
    static_assert(SIZE <= COLS);
    assert(row + SIZE <= ROWS);
    assert(col <= COLS);

    using input_matrix_t = matrix<array_t, scalar_t, ROWS, COLS>;
    using vector_t = vector<SIZE, scalar_t, array_t>;

    vector_t res_v{};

    // Don't access single elements in underlying vectors unless necessary
    if constexpr (SIZE == input_matrix_t::storage_rows()) {
      if (row == 0u) {
        return m[col];
      }
    }
    for (std::size_t i = row; i < row + SIZE; ++i) {
      res_v[i] = m[col][i];
    }

    return res_v;
  }

};  // struct block_getter

/// Get a block of a const matrix
template <std::size_t ROWS, std::size_t COLS, std::size_t mROW,
          std::size_t mCOL, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr auto block(
    const matrix<array_t, scalar_t, mROW, mCOL> &m, const std::size_t row,
    const std::size_t col) noexcept {
  return block_getter{}.template operator()<ROWS, COLS>(m, row, col);
}

/// Get a block of a const matrix
template <std::size_t ROWS, std::size_t COLS, std::size_t mROW,
          std::size_t mCOL, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr void set_block(
    matrix<array_t, scalar_t, mROW, mCOL> &m,
    const matrix<array_t, scalar_t, ROWS, COLS> &b, const std::size_t row,
    const std::size_t col) noexcept {
  static_assert(ROWS <= mROW);
  static_assert(COLS <= mCOL);
  assert(row + ROWS <= mROW);
  assert(col + COLS <= mCOL);

  using input_matrix_t = matrix<array_t, scalar_t, mROW, mCOL>;
  using matrix_t = matrix<array_t, scalar_t, ROWS, COLS>;

  // Don't access single elements in underlying vectors unless necessary
  if constexpr (matrix_t::storage_rows() == input_matrix_t::storage_rows()) {
    if (row == 0u) {
      for (std::size_t j = col; j < mCOL; ++j) {
        m[j] = b[j - col];
      }
      return;
    }
  }
  for (std::size_t j = col; j < col + COLS; ++j) {
    for (std::size_t i = row; i < row + ROWS; ++i) {
      m[j][i] = b[j - col][i - row];
    }
  }
}

/// Operator setting a block with a vector
template <class matrix_t, std::size_t ROW, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE constexpr void set_block(
    matrix_t &m, const storage::vector<ROW, scalar_t, array_t> &b,
    std::size_t row, std::size_t col) noexcept {
  assert(row < ROW);
  assert(row < matrix_t::rows());
  assert(col < matrix_t::columns());

  if constexpr (matrix_t::storage_rows() == ROW) {
    if (row == 0u) {
      m[col] = b;
      return;
    }
  }
  for (std::size_t i = row; i < ROW; ++i) {
    m[col][i] = b[i - row];
  }
}

}  // namespace algebra::storage
