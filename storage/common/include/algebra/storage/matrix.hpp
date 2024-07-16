/** Algebra plugins, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s).
#include "algebra/storage/vector.hpp"

// System include(s).
#include <array>
#include <cassert>

namespace algebra::storage {

/// Generic matrix type that can take vectors as columns
template <template <typename, std::size_t> class array_t, typename value_t,
          std::size_t ROW, std::size_t COL>
struct matrix {

  // The matrix consists of column vectors
  using vector_type = storage::vector<ROW, value_t, array_t>;
  // Value type: Can be simd types
  using value_type = value_t;

  /// Default constructor: Zero matrix
  constexpr matrix() : m_storage{} {}

  /// Construct from given column vectors @param v
  template <typename... vector_t,
            std::enable_if_t<sizeof...(vector_t) == COL, bool> = true>
  explicit matrix(vector_t &&...v) : m_storage{std::forward<vector_t>(v)...} {}

  /// Equality operator between two matrices
  constexpr bool operator==(const matrix &rhs) const {
    return equal(rhs, std::make_index_sequence<COL>());
  }

  /// Subscript operator
  /// @{
  constexpr decltype(auto) operator[](const std::size_t i) const {
    return m_storage[i];
  }
  constexpr decltype(auto) operator[](const std::size_t i) {
    return m_storage[i];
  }
  /// @}

  /// @returns the number of rows
  static constexpr std::size_t rows() { return ROW; }

  /// @returns the number of rows of the underlying vector storage
  static constexpr std::size_t storage_rows() { return vector_type::size(); }

  /// @returns the number of columns
  static constexpr std::size_t columns() { return COL; }

 private:
  /// Sets the trailing uninitialized values to zero.
  /// @{
  // AoS
  template <typename V = value_t, std::size_t... I,
            typename std::enable_if_t<!std::is_scalar_v<V>, bool> = true>
  constexpr bool equal(const matrix &rhs, std::index_sequence<I...>) const {
    return (... && (m_storage[I] == rhs[I]));
  }

  // SoA
  template <typename V = value_t, std::size_t... I,
            typename std::enable_if_t<std::is_scalar_v<V>, bool> = true>
  constexpr bool equal(const matrix &rhs, std::index_sequence<I...>) const {
    return (... && ((m_storage[I].get() == rhs[I].get()).isFull()));
  }
  /// @}

  /// Matrix storage
  std::array<vector_type, COL> m_storage;

};  // struct matrix

/// Functor used to access elements of a matrix
struct element_getter {

  /// Get const access to a matrix element
  template <template <typename, std::size_t> class array_t, typename value_t,
            std::size_t ROW, std::size_t COL>
  ALGEBRA_HOST inline decltype(auto) operator()(
      const matrix<array_t, value_t, ROW, COL> &m, std::size_t row,
      std::size_t col) const {

    // Make sure that the indices are valid.
    assert(row < ROW);
    assert(col < COL);

    // Return the selected element.
    return m[col][row];
  }

  /// Get non-const access to a matrix element
  template <template <typename, std::size_t> class array_t, typename value_t,
            std::size_t ROW, std::size_t COL>
  ALGEBRA_HOST inline decltype(auto) operator()(
      matrix<array_t, value_t, ROW, COL> &m, std::size_t row,
      std::size_t col) const {

    // Make sure that the indices are valid.
    assert(row < ROW);
    assert(col < COL);

    // Return the selected element.
    return m[col][row];
  }

  /// Get const access to a matrix44 element from a flat array
  template <template <typename, std::size_t> class array_t, typename value_t>
  ALGEBRA_HOST inline decltype(auto) operator()(const array_t<value_t, 16> &m,
                                                unsigned int row,
                                                unsigned int col) const {
    // Make sure that the indices are valid.
    assert(row < 4);
    assert(col < 4);
    // Return the selected element.
    return m[row * 4 + col];
  }

};  // struct element_getter

/// Build a zero matrix
template <typename matrix_t>
constexpr matrix_t zero() noexcept {
  return {};
}

/// Build an identity matrix
template <typename matrix_t, std::size_t... I>
constexpr matrix_t identity(std::index_sequence<I...>) noexcept {

  // Zero initialized
  matrix_t m{};

  ((m[I][I] = typename matrix_t::value_type(1)), ...);

  return m;
}

/// Build an identity matrix
template <typename matrix_t>
constexpr matrix_t identity() noexcept {

  return identity<matrix_t>(std::make_index_sequence<std::min(
                                matrix_t::rows(), matrix_t::columns())>());
}

/// Transpose the matrix @param m
template <std::size_t ROW, std::size_t COL, typename value_t,
          template <typename, std::size_t> class array_t>
constexpr auto transpose(const matrix<array_t, value_t, ROW, COL> &m) noexcept {

  matrix<array_t, value_t, COL, ROW> res_m;

  return res_m;
}

/// Get a block of a const matrix
/*template <std::size_t ROWS, std::size_t COLS, typename matrix_t>
ALGEBRA_HOST_DEVICE inline constexpr auto block(const matrix_t &m,
                                                  std::size_t row,
                                                  std::size_t col) noexcept {
  return block_getter().template operator()<ROWS, COLS>(m, row, col);
}

/// Operator setting a block with a vector matrix
template <std::size_t ROWS, std::size_t COLS, class input_matrix_type>
ALGEBRA_HOST_DEVICE inline constexpr void set_block(input_matrix_type &m,
                                    const matrix_type<ROWS, COLS> &b,
                                    std::size_t row, std::size_t col) noexcept {
  for (std::size_t j = 0; j < COLS; ++j) {
    for (std::size_t i = 0; i < ROWS; ++i) {
      element_getter()(m, i + row, j + col) = element_getter()(b, i, j);
    }
  }
}

/// Operator setting a block with a vector
template <size_t ROWS, template <typename, size_t> class vector_t,
          class input_matrix_type>
ALGEBRA_HOST_DEVICE inline constexpr void set_block(input_matrix_type &m,
                                    const vector_t<scalar_t, ROWS> &b,
                                    std::size_t row, std::size_t col) noexcept {
  for (std::size_t i = 0; i < ROWS; ++i) {
    element_getter()(m, i + row, col) = b[i];
  }
}*/

/// Scalar multiplication
template <std::size_t ROW, std::size_t COL, typename value_t, typename scalar_t,
          template <typename, std::size_t> class array_t,
          std::enable_if_t<std::is_scalar_v<scalar_t> ||
                               std::is_same_v<value_t, scalar_t>,
                           bool> = true>
inline constexpr decltype(auto) operator*(
    scalar_t a, const matrix<array_t, value_t, ROW, COL> &rhs) noexcept {

  matrix<array_t, value_t, ROW, COL> res_m;

  for (std::size_t j = 0; j < COL; ++j) {
    res_m[j] = a * rhs[j];
  }
  return res_m;
}
template <std::size_t ROW, std::size_t COL, typename value_t, typename scalar_t,
          template <typename, std::size_t> class array_t,
          std::enable_if_t<std::is_scalar_v<scalar_t> ||
                               std::is_same_v<value_t, scalar_t>,
                           bool> = true>
inline constexpr decltype(auto) operator*(
    const matrix<array_t, value_t, ROW, COL> &lhs, scalar_t a) noexcept {
  return a * lhs;
}

/// Matrix addition
template <std::size_t ROW, std::size_t COL, typename value_t,
          template <typename, std::size_t> class array_t>
inline constexpr decltype(auto) operator+(
    const matrix<array_t, value_t, ROW, COL> &lhs,
    const matrix<array_t, value_t, ROW, COL> &rhs) noexcept {

  matrix<array_t, value_t, ROW, COL> res_m;

  for (std::size_t j = 0; j < COL; ++j) {
    res_m[j] = lhs[j] + rhs[j];
  }

  return res_m;
}

template <std::size_t ROW, std::size_t COL, typename value_t,
          template <typename, std::size_t> class array_t>
inline constexpr decltype(auto) operator-(
    const matrix<array_t, value_t, ROW, COL> &lhs,
    const matrix<array_t, value_t, ROW, COL> &rhs) noexcept {

  matrix<array_t, value_t, ROW, COL> res_m;

  for (std::size_t j = 0; j < COL; ++j) {
    res_m[j] = lhs[j] - rhs[j];
  }

  return res_m;
}

/// Matrix-vector multiplication
template <std::size_t ROW, std::size_t COL, typename value_t,
          template <typename, std::size_t> class array_t>
inline constexpr decltype(auto) operator*(
    const matrix<array_t, value_t, ROW, COL> &lhs,
    const vector<COL, value_t, array_t> &v) noexcept {

  // Init vector
  vector<ROW, value_t, array_t> res_v{v[0] * lhs[0]};

  // Add the rest per column
  for (std::size_t j = 1u; j < COL; ++j) {
    res_v = res_v + v[j] * lhs[j];
  }

  return res_v;
}

/// Matrix-matrix multiplication
template <std::size_t LROW, std::size_t COL, std::size_t RCOL, typename value_t,
          template <typename, std::size_t> class array_t>
inline constexpr decltype(auto) operator*(
    const matrix<array_t, value_t, LROW, COL> &lhs,
    const matrix<array_t, value_t, COL, RCOL> &rhs) noexcept {

  matrix<array_t, value_t, LROW, RCOL> res_m;

  for (std::size_t j = 0u; j < RCOL; ++j) {
    // Init column j
    res_m[j] = rhs[j][0] * lhs[0];

    // Add the rest per column
    for (std::size_t i = 1u; i < COL; ++i) {
      res_m[j] = res_m[j] + rhs[j][i] * lhs[i];
    }
  }

  return res_m;
}

}  // namespace algebra::storage
