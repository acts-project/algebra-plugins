/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
#include "algebra/qualifiers.hpp"

namespace algebra::cmath::matrix::inverse {

/// "Partial Pivot LU Decomposition", assuming a N X N matrix
template <typename size_type,
          template <typename, size_type, size_type> class matrix_t,
          typename scalar_t, class element_getter_t, size_type... Ds>
struct partial_pivot_lud {

  using _dims = std::integer_sequence<size_type, Ds...>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = matrix_t<scalar_t, ROWS, COLS>;

  using decomposition_t =
      typename algebra::cmath::matrix::decomposition::partial_pivot_lud<
          size_type, matrix_t, scalar_t, element_getter_t>;

  template <size_type N>
  ALGEBRA_HOST_DEVICE inline matrix_type<N, N> operator()(
      const matrix_type<N, N>& m) const {
    // Get the LUDecomposition matrix in the form of (L - I) + U
    const typename decomposition_t::template lud<N> decomp_res =
        decomposition_t()(m);

    const auto& lu = decomp_res.lu;
    const auto& p_mat = decomp_res.p_mat;

    // Inverse matrix
    matrix_type<N, N> inv;

    for (size_type j = 0; j < N; j++) {
      for (size_type i = 0; i < N; i++) {
        if (element_getter_t()(p_mat, i, j) == 1) {
          element_getter_t()(inv, i, j) = 1.0;
        } else {
          element_getter_t()(inv, i, j) = 0.0;
        }

        for (size_type k = 0; k < i; k++) {
          element_getter_t()(inv, i, j) -=
              element_getter_t()(lu, i, k) * element_getter_t()(inv, k, j);
        }
      }

      for (size_type i = N - 1; int(i) >= 0; i--) {
        for (size_type k = i + 1; k < N; k++) {
          element_getter_t()(inv, i, j) -=
              element_getter_t()(lu, i, k) * element_getter_t()(inv, k, j);
        }
        element_getter_t()(inv, i, j) /= element_getter_t()(lu, i, i);
      }
    }

    return inv;
  }
};

}  // namespace algebra::cmath::matrix::inverse