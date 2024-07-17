/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"
#include "algebra/storage/matrix.hpp"

namespace algebra::vc_aos {

namespace matrix {

/// Explicitly vectorized matrix implementation
template <template <typename, std::size_t> class array_t, typename scalar_t>
struct actor {

  /// Size type
  using size_type = std::size_t;
  using size_ty = size_type;

  /// Scalar type
  using scalar_type = scalar_t;

  /// The type of the matrix elements (scalar for AoS, Vc::Vector for SoA)
  using value_type =
      std::conditional_t<Vc::is_simd_vector<array_t<scalar_t, 4>>::value,
                         scalar_t, Vc::Vector<scalar_t>>;
  /// 2D matrix type
  template <std::size_t ROWS, std::size_t COLS>
  using matrix_type = storage::matrix<array_t, value_type, ROWS, COLS>;

  /// vector type
  // template <std::size_t N>
  // using array_type = array_t<value_type, N>;

  /// Operator getting a reference to one element of a non-const matrix
  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE static constexpr value_type &element(
      matrix_type<ROWS, COLS> &m, std::size_t row, std::size_t col) {
    return m[col][row];
  }

  /// Operator getting one value of a const matrix
  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE static constexpr const value_type &element(
      const matrix_type<ROWS, COLS> &m, std::size_t row, std::size_t col) {
    return m[col][row];
  }

  /// Operator getting a block of a const matrix
  template <std::size_t ROWS, std::size_t COLS, class input_matrix_t>
  ALGEBRA_HOST_DEVICE static constexpr matrix_type<ROWS, COLS> block(
      const input_matrix_t &m, std::size_t row, std::size_t col) {
    return storage::block<ROWS, COLS>(m, row, col);
  }

  /// Operator setting a block with a vector matrix
  template <std::size_t ROWS, std::size_t COLS, class input_matrix_t>
  ALGEBRA_HOST_DEVICE static constexpr void set_block(
      input_matrix_t &m, const matrix_type<ROWS, COLS> &b, std::size_t row,
      std::size_t col) {
    storage::set_block(m, b, row, col);
  }

  /// Operator setting a block with a vector
  template <class input_matrix_t, class vector_t>
  ALGEBRA_HOST_DEVICE static constexpr void set_block(input_matrix_t &m,
                                                      const vector_t &b,
                                                      std::size_t row,
                                                      std::size_t col) {
    storage::set_block(m, b, row, col);
  }

  // Create zero matrix
  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE static constexpr auto zero() {
    return storage::zero<matrix_type<ROWS, COLS>>();
  }

  // Create identity matrix
  template <size_t ROWS, size_t COLS>
  ALGEBRA_HOST_DEVICE static constexpr auto identity() {
    return storage::identity<matrix_type<ROWS, COLS>>();
  }

  // Set input matrix as zero matrix
  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE static constexpr void set_zero(
      matrix_type<ROWS, COLS> &m) {
    m = zero<ROWS, COLS>();
  }

  // Set input matrix as zero matrix
  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE static constexpr void set_identity(
      matrix_type<ROWS, COLS> &m) {
    m = identity<ROWS, COLS>();
  }

  /// Create transpose matrix
  template <std::size_t ROWS, std::size_t COLS>
  ALGEBRA_HOST_DEVICE static constexpr auto transpose(
      const matrix_type<ROWS, COLS> &m) {
    return storage::transpose(m);
  }

  /// Get determinant using a specific algorithm @tparam determinant_t
  /*template <std::size_t N, typename determinant_t = determinant_t>
  ALGEBRA_HOST_DEVICE static constexpr constexpr value_type determinant(
      const matrix_type<N, N> &m) {

    return determinant_t{}(m);
  }

  /// Get inverse using a specific algorithm @tparam inverse_t
  template <std::size_t N, typename inverse_t = inverse_t>
  ALGEBRA_HOST_DEVICE static constexpr constexpr matrix_type<N, N> inverse(
      const matrix_type<N, N> &m) {

    return inverse_t{}(m);
  }*/
};

}  // namespace matrix

namespace math {}  // namespace math
}  // namespace algebra::vc_aos
