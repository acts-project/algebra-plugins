/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/matrix/determinant/hard_coded.hpp"
#include "algebra/math/algorithms/matrix/determinant/partial_pivot_lud.hpp"
#include "algebra/math/algorithms/matrix/inverse/hard_coded.hpp"
#include "algebra/math/algorithms/matrix/inverse/partial_pivot_lud.hpp"

namespace algebra::generic {

/// Get the type of determinant algorithm acording to matrix dimension by
/// partial template specilization
/// @{

// Default algorithm
template <std::size_t N, concepts::square_matrix M>
struct determinant_selector {
  using type = matrix::determinant::partial_pivot_lud<M>;
};

/// Always use hard coded implementation for very small matrices
template <concepts::square_matrix M>
struct determinant_selector<2, M> {
  using type = matrix::determinant::hard_coded<M>;
};
/// @}

/// Get the type of inversion algorithm acording to matrix dimension by
/// partial template specilization
/// @{

// Default algorithm
template <std::size_t N, concepts::square_matrix M>
struct inversion_selector {
  using type = matrix::inverse::partial_pivot_lud<M>;
};

/// Always use hard coded implementation for very small matrices
template <concepts::square_matrix M>
struct inversion_selector<2, M> {
  using type = matrix::inverse::hard_coded<M>;
};
/// @}

}  // namespace algebra::generic
