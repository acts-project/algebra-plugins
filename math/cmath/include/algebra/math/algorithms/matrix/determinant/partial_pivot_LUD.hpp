/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

namespace algebra::cmath::matrix::determinant {

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
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const matrix_type<N, N>& m) const {
    // Reference: https://en.wikipedia.org/wiki/LU_decomposition

    // Get the LUDecomposition matrix in the form of (L - I) + U
    const LUD<N> decomp_res = decomposition_t()(m);

    const auto& lu = decomp_res.lu;
    const auto& n_pivot = decomp_res.n_pivot;

    scalar_t det = element_getter_t()(lu, 0, 0);

    for (std::size_t i = 0; i < N; i++) {
      det *= element_getter_t()(lu, i, i);
    }

    return (n_pivot - N) % 2 == 0 ? det : -det;
  }
};

}  // namespace algebra::cmath::matrix::determinant