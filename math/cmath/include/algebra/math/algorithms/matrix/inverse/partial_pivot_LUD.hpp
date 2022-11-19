/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

namespace algebra::cmath::matrix::inverse {

/// "Partial Pivot LU Decomposition", assuming a N X N matrix
template <typename size_type,
          template <typename, size_type, size_type> class matrix_t,
          typename scalar_t, class element_getter_t, size_type... Ds>
struct partial_pivot_LUD {

  using _dims = std::integer_sequence<size_type, Ds...>;

  /// Function (object) used for accessing a matrix element
  using element_getter = element_getter_t;

  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = matrix_t<scalar_t, ROWS, COLS>;

  using decomposition_t =
      typename algebra::cmath::matrix::decomposition::partial_pivot_LUD<
          size_type, matrix_t, scalar_t, element_getter_t>;

  template <size_type N>
  ALGEBRA_HOST_DEVICE inline matrix_type<N, N> operator()(
      const matrix_type<N, N>& m) const {
    // Reference: https://en.wikipedia.org/wiki/LU_decomposition

    // Get the LUDecomposition matrix in the form of (L - I) + U
    const LUD<N> decomp_res = decomposition_t()(m);

    const auto& lu = decomp_res.lu;
    const auto& P = decomp_res.P;

    // Inverse matrix
    matrix_type<N, N> Inv;

    for (std::size_t j = 0; j < N; j++) {
      for (std::size_t i = 0; i < N; i++) {
        // IA[i][j] = P[i] == j ? 1.0 : 0.0;
        element_getter_t()(inv, i, j) = P[i] == j ? 1.0 : 0.0;

        for (std::size_t k = 0; k < i; k++) {
          // IA[i][j] -= A[i][k] * IA[k][j];
          element_getter_t()(inv, i, j) -=
              element_getter_t()(lu, j, k) * element_getter_t()(inv, k, j);
        }
      }

      for (std::size_t i = N - 1; i >= 0; i--) {
        for (std::size_t k = i + 1; k < N; k++) {
          // IA[i][j] -= A[i][k] * IA[k][j];
          element_getter_t()(inv, i, j) -=
              element_getter_t()(lu, i, k) * element_getter_t()(inv, k, j);
        }
        // IA[i][j] /= A[i][i];
        element_getter_t()(inv, i, j) /= element_getter_t()(lu, i, i);
      }
    }

    return Inv;
  }

}  // namespace algebra::cmath::matrix::inverse