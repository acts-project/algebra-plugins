/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/type_traits.hpp"

namespace algebra::cmath::matrix::determinant {

/// "Partial Pivot LU Decomposition", assuming a N X N matrix
template <class matrix_t, class element_getter_t>
struct partial_pivot_lud {

  using scalar_type = algebra::trait::value_t<matrix_t>;
  using size_type = algebra::trait::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  using decomposition_t =
      typename algebra::cmath::matrix::decomposition::partial_pivot_lud<
          matrix_t, element_getter_t>;

  ALGEBRA_HOST_DEVICE inline scalar_type operator()(const matrix_t& m) const {

    constexpr size_type N{algebra::trait::rank<matrix_t>};

    const typename decomposition_t::template lud<algebra::trait::rank<matrix_t>>
        decomp_res = decomposition_t()(m);

    // Get the LU decomposition matrix equal to (L - I) + U
    const auto& lu = decomp_res.lu;
    const auto n_pivot = static_cast<size_type>(decomp_res.n_pivot);

    scalar_type det = element_getter()(lu, 0, 0);

    for (size_type i = 1; i < N; i++) {
      det *= element_getter()(lu, i, i);
    }

    return (n_pivot - N) % 2 == 0 ? det : -det;
  }
};

}  // namespace algebra::cmath::matrix::determinant
