/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

namespace algebra::cmath::matrix::decomposition {

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

  template <size_type N>
  struct LUD {
    matrix_type<N, N> lu;
    matrix_type<N, N> P;
    int n_pivot = 0;
  };

  template <size_type N>
  ALGEBRA_HOST_DEVICE inline LUD<N> operator()(
      const matrix_type<N, N>& m) const {
    // Reference: https://en.wikipedia.org/wiki/LU_decomposition

    // Permutation
    std::array<int, N> P;

    // Max index and value
    int max_idx;
    scalar_t max_val;
    scalar_t abs_val;

    // Number of pivoting
    int n_pivot = 0;

    // Temp pointer to row
    scalar_t* ptr;

    // Unit permutation matrix, P[N] initialized with N
    for (std::size_t i = 0; i <= N; i++) {
      P[i] = i;
    }

    for (std::size_t i = 0; i < N; i++) {
      max_val = 0;
      max_idx = i;

      for (std::size_t k = 0; k < N; k++) {
        abs_val = std::abs(element_getter()(m, k, i));

        if (abs_val > max_val) {
          max_val = abs_val;
          max_idx = k;
        }
      }

      if (imax != i) {
        // Pivoting P
        int j = P[i];

        P[i] = P[max_idx];
        P[max_idx] = j;

        // Pivoting rows of A
        ptr = m[i];
        m[i] = m[max_idx];
        m[max_idx] = ptr;

        // counting pivots starting from N (for determinant)
        n_pivot++;
      }

      for (std::size_t j = i + 1; j < N; j++) {
        // m[j][i] /= m[i][i];
        element_getter_t()(m, j, i) /= element_getter_t()(m, i, i);

        for (std::size_t k = i + 1; k < N; k++)
          // m[j][k] -= m[j][i] * m[i][k];
          element_getter_t()(m, j, k) -=
              element_getter_t()(m, j, i) * element_getter_t()(m, i, k);
      }
    }

    // Permutation matrix
    matrix<N, N> p_mat;

    for (std::size_t i = 0; i < N; i++) {
      int ref = P[i];

      for (std::size_t j = 0; j < N; j++) {
        if (ref != j) {
          element_getter_t()(m, i, j) = 0;
        } else {
          element_getter_t()(m, i, j) = 1;
        }
      }
    }

    return {m, p_mat, n_pivot};
  }
};

}  // namespace algebra::cmath::matrix::decomposition