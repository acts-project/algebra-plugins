/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"
#include "algebra/storage/matrix_getter.hpp"

namespace algebra::vc_aos::storage {

using algebra::storage::block;
using algebra::storage::element;
using algebra::storage::set_block;

/// Get a vector of a const matrix
template <std::size_t SIZE, std::size_t ROW, std::size_t COL, typename value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline constexpr auto vector(
    const algebra::storage::matrix<array_t, value_t, ROW, COL> &m,
    const std::size_t row, const std::size_t col) noexcept {
  return algebra::storage::block_getter{}.template operator()<SIZE>(m, row,
                                                                    col);
}

}  // namespace algebra::vc_aos::storage
