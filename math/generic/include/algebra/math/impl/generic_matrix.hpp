/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/algorithms/utils/algorithm_finder.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::generic::math {

/// Create zero matrix - generic transform3
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE inline M zero() {

  using index_t = algebra::traits::index_t<M>;
  using element_getter_t = algebra::traits::element_getter_t<M>;

  M ret;

  for (index_t j = 0; j < algebra::traits::columns<M>; ++j) {
    for (index_t i = 0; i < algebra::traits::rows<M>; ++i) {
      element_getter_t{}(ret, i, j) = 0;
    }
  }

  return ret;
}

/// Create identity matrix - generic transform3
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE inline M identity() {

  using index_t = algebra::traits::index_t<M>;
  using element_getter_t = algebra::traits::element_getter_t<M>;

  auto ret{zero<M>()};

  for (index_t i = 0; i < algebra::traits::rank<M>; ++i) {
    element_getter_t{}(ret, i, i) = 1;
  }

  return ret;
}

/// Set @param m as zero matrix
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE inline void set_zero(M &m) {
  m = zero<M>();
}

/// Set @param m as identity matrix
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE inline void set_identity(M &m) {
  m = identity<M>();
}

/// @returns the transpose matrix of @param m
template <concepts::matrix M>
ALGEBRA_HOST_DEVICE inline auto transpose(const M &m) {

  using index_t = algebra::traits::index_t<M>;
  using value_t = algebra::traits::value_t<M>;
  using element_getter_t = algebra::traits::element_getter_t<M>;

  constexpr index_t rows{algebra::traits::rows<M>};
  constexpr index_t columns{algebra::traits::columns<M>};

  algebra::traits::get_matrix_t<M, columns, rows, value_t> ret;

  for (index_t i = 0; i < rows; ++i) {
    for (index_t j = 0; j < columns; ++j) {
      element_getter_t{}(ret, j, i) = element_getter_t{}(m, i, j);
    }
  }

  return ret;
}

/// @returns the determinant of @param m
template <concepts::square_matrix M>
ALGEBRA_HOST_DEVICE inline algebra::traits::scalar_t<M> determinant(
    const M &m) {

  return determinant_t<M>{}(m);
}

/// @returns the determinant of @param m
template <concepts::square_matrix M>
ALGEBRA_HOST_DEVICE inline M inverse(const M &m) {

  return inversion_t<M>{}(m);
}

}  // namespace algebra::generic::math
