/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/matrix.hpp"

namespace algebra::vc_aos::math {

using storage::identity;
using storage::set_identity;
using storage::set_zero;
using storage::transpose;
using storage::zero;

/// General case: Compute the determinant of a square matrix
///
/// @returns the determinant
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline constexpr value_t determinant(
    const storage::matrix<array_t, value_t, N, N> &m) noexcept {

  return 0;
}

template <std::size_t ROW, std::size_t COL, concepts::value value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline constexpr storage::matrix<array_t, value_t, ROW, COL>
inverse(const storage::matrix<array_t, value_t, ROW, COL> &m) noexcept {

  return m;
}

/// Transpose(matrix)-matrix multiplication
/*template <std::size_t LROW, std::size_t COL, std::size_t RCOL, concepts::value
value_t, template <typename, std::size_t> class array_t> ALGEBRA_HOST_DEVICE
inline constexpr decltype(auto) transpose_mul( const matrix<array_t, value_t,
LROW, COL> &lhs, const matrix<array_t, value_t, COL, RCOL> &rhs) noexcept {

  matrix<array_t, value_t, LROW, RCOL> res_m;

  for (std::size_t j = 0u; j < RCOL; ++j) {
    // Add the rest per column
    for (std::size_t i = 1u; i < COL; ++i) {
      res_m[j][i] = vector::sum(lhs[i] * rhs[j]);
    }
  }

  return res_m;
}*/

}  // namespace algebra::vc_aos::math
