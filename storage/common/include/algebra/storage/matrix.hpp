/** Algebra plugins, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s).
#include "algebra/storage/vector.hpp"
#include "algebra/type_traits.hpp"

// System include(s).
#include <array>
#include <cassert>

namespace algebra::storage {

/// Generic matrix type that can take vectors as columns
template <template <typename, std::size_t> class array_t, typename value_t,
          std::size_t ROW, std::size_t COL>
struct alignas(alignof(storage::vector<ROW, value_t, array_t>)) matrix {

  // The matrix consists of column vectors
  using vector_type = storage::vector<ROW, value_t, array_t>;
  // Value type: Can be simd types
  using value_type = value_t;

  /// Default constructor
  constexpr matrix() = default;

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
  constexpr const vector_type &operator[](const std::size_t i) const {
    assert(i < COL);
    return m_storage[i];
  }
  constexpr vector_type &operator[](const std::size_t i) {
    assert(i < COL);
    return m_storage[i];
  }
  /// @}

  /// @returns the number of rows
  static constexpr std::size_t rows() { return ROW; }

  /// @returns the number of rows of the underlying vector storage
  static constexpr std::size_t storage_rows() {
    return vector_type::simd_size();
  }

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

/// Get a zero-initialized matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE constexpr matrix_t zero() noexcept {

  matrix_t m;

  for (std::size_t j = 0u; j < matrix_t::columns(); ++j) {
    // Fill zero initialized vector
    m[j] = typename matrix_t::vector_type{};
  }

  return m;
}

/// Set a matrix to zero
template <typename matrix_t>
ALGEBRA_HOST_DEVICE constexpr void set_zero(matrix_t &m) noexcept {
  m = zero<matrix_t>();
}

/// Build an identity matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE constexpr matrix_t identity() noexcept {

  // Zero initialized
  matrix_t m{zero<matrix_t>()};

  for (std::size_t i = 0u; i < algebra::trait::rank<matrix_t>; ++i) {
    m[i][i] = typename matrix_t::value_type(1);
  }

  return m;
}

/// Set a matrix to zero
template <typename matrix_t>
ALGEBRA_HOST_DEVICE constexpr void set_identity(matrix_t &m) noexcept {
  m = identity<matrix_t>();
}

/// Transpose the matrix @param m
template <std::size_t ROW, std::size_t COL, typename value_t,
          template <typename, std::size_t> class array_t, std::size_t... I>
ALGEBRA_HOST_DEVICE constexpr auto transpose(
    const matrix<array_t, value_t, ROW, COL> &m,
    std::index_sequence<I...>) noexcept {

  using matrix_T_t = matrix<array_t, value_t, COL, ROW>;
  using column_t = typename matrix_T_t::vector_type;

  matrix_T_t res_m;

  for (std::size_t j = 0u; j < ROW; ++j) {
    res_m[j] = column_t{m[I][j]...};
  }

  return res_m;
}

/// Build an identity matrix
template <typename matrix_t>
ALGEBRA_HOST_DEVICE constexpr auto transpose(const matrix_t &m) noexcept {
  return transpose(m, std::make_index_sequence<matrix_t::columns()>());
}

/// Scalar multiplication
template <typename matrix_t, typename scalar_t, std::size_t... J,
          std::enable_if_t<
              std::is_scalar_v<scalar_t> ||
                  std::is_same_v<typename matrix_t::value_type, scalar_t>,
              bool> = true>
ALGEBRA_HOST_DEVICE inline constexpr matrix_t matrix_scalar_mul(
    scalar_t a, const matrix_t &rhs, std::index_sequence<J...>) noexcept {

  return matrix_t{(a * rhs[J])...};
}

template <std::size_t ROW, std::size_t COL, typename value_t, typename scalar_t,
          template <typename, std::size_t> class array_t,
          std::enable_if_t<std::is_scalar_v<scalar_t> ||
                               std::is_same_v<value_t, scalar_t>,
                           bool> = true>
ALGEBRA_HOST_DEVICE inline constexpr decltype(auto) operator*(
    scalar_t a, const matrix<array_t, value_t, ROW, COL> &rhs) noexcept {

  using matrix_t = matrix<array_t, value_t, ROW, COL>;

  return matrix_scalar_mul(a, rhs,
                           std::make_index_sequence<matrix_t::columns()>());
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
template <typename matrix_t, std::size_t... J>
ALGEBRA_HOST_DEVICE inline constexpr matrix_t matrix_add(
    const matrix_t &lhs, const matrix_t &rhs,
    std::index_sequence<J...>) noexcept {

  return matrix_t{(lhs[J] + rhs[J])...};
}

template <std::size_t ROW, std::size_t COL, typename value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline constexpr decltype(auto) operator+(
    const matrix<array_t, value_t, ROW, COL> &lhs,
    const matrix<array_t, value_t, ROW, COL> &rhs) noexcept {

  using matrix_t = matrix<array_t, value_t, ROW, COL>;

  return matrix_add(lhs, rhs, std::make_index_sequence<matrix_t::columns()>());
}

template <typename matrix_t, std::size_t... J>
ALGEBRA_HOST_DEVICE inline constexpr decltype(auto) matrix_sub(
    const matrix_t &lhs, const matrix_t &rhs,
    std::index_sequence<J...>) noexcept {

  return matrix_t{(lhs[J] - rhs[J])...};
}

template <std::size_t ROW, std::size_t COL, typename value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline constexpr decltype(auto) operator-(
    const matrix<array_t, value_t, ROW, COL> &lhs,
    const matrix<array_t, value_t, ROW, COL> &rhs) noexcept {

  using matrix_t = matrix<array_t, value_t, ROW, COL>;

  return matrix_sub(lhs, rhs, std::make_index_sequence<matrix_t::columns()>());
}

/// Matrix-vector multiplication
template <std::size_t ROW, std::size_t COL, typename value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline constexpr decltype(auto) operator*(
    const matrix<array_t, value_t, ROW, COL> &lhs,
    const vector<COL, value_t, array_t> &v) noexcept {

  // Init vector
  vector<ROW, value_t, array_t> res_v{v[0] * lhs[0]};

  // Add the rest per column
  for (std::size_t j = 1u; j < COL; ++j) {
    // fma
    res_v = res_v + v[j] * lhs[j];
  }

  return res_v;
}

/// Matrix-matrix multiplication
template <std::size_t LROW, std::size_t COL, std::size_t RCOL, typename value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline constexpr decltype(auto) operator*(
    const matrix<array_t, value_t, LROW, COL> &lhs,
    const matrix<array_t, value_t, COL, RCOL> &rhs) noexcept {

  matrix<array_t, value_t, LROW, RCOL> res_m;

  for (std::size_t j = 0u; j < RCOL; ++j) {
    // Init column j
    res_m[j] = rhs[j][0] * lhs[0];

    // Add the rest per column
    for (std::size_t i = 1u; i < COL; ++i) {
      // fma
      res_m[j] = res_m[j] + rhs[j][i] * lhs[i];
    }
  }

  return res_m;
}

}  // namespace algebra::storage
