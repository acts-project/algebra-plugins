/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/type_traits.hpp"

namespace algebra::cmath {

/// @returns zero matrix of type @tparam matrix_t
template <typename matrix_t>
requires(std::is_scalar_v<typename matrix_t::value_type::value_type>)
    ALGEBRA_HOST_DEVICE inline matrix_t zero() {
  return matrix_t{};
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

/// Set @param m as zero matrix
template <std::size_t ROWS, std::size_t COLS, typename scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline constexpr void set_zero(
    array_t<array_t<scalar_t, ROWS>, COLS> &m) {
  m = zero<array_t<array_t<scalar_t, ROWS>, COLS>>();
}

/// Set @param m as identity matrix
template <std::size_t ROWS, std::size_t COLS, typename scalar_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline constexpr void set_identity(
    array_t<array_t<scalar_t, ROWS>, COLS> &m) {
  m = identity<array_t<array_t<scalar_t, ROWS>, COLS>>();
}

}  // namespace algebra::cmath
