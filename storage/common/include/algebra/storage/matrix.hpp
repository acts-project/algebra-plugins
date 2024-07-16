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

/// Build an identity matrix
template <typename matrix_t, std::size_t... I>
constexpr matrix_t identity(std::index_sequence<I...>) {

  // Zero initialized
  matrix_t m{};

  ((m[I][I] = typename matrix_t::value_type(1)), ...);

  return m;
}

/// Build an identity matrix
template <typename matrix_t>
constexpr matrix_t identity() {

  return identity<matrix_t>(std::make_index_sequence<std::min(
                                matrix_t::rows(), matrix_t::columns())>());
}

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

  /// Get const access to a matrix element
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

}  // namespace algebra::storage
