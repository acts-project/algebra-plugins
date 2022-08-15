/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s).
#include <type_traits>

namespace algebra::eigen::math {

/// Functor used to access elements of Eigen matrices
struct element_getter {
  /// Get non-const access to a matrix element
  template <
      typename derived_type,
      std::enable_if_t<std::is_base_of<Eigen::DenseCoeffsBase<
                                           derived_type, Eigen::WriteAccessors>,
                                       Eigen::MatrixBase<derived_type> >::value,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline auto &operator()(
      Eigen::MatrixBase<derived_type> &m, int row, int col) const {

    return m(row, col);
  }
  /// Get const access to a matrix element
  template <typename derived_type>
  ALGEBRA_HOST_DEVICE inline auto operator()(
      const Eigen::MatrixBase<derived_type> &m, int row, int col) const {

    return m(row, col);
  }
};  // struct element_getter

/// Function extracting an element from a matrix (const)
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline auto element(
    const Eigen::MatrixBase<derived_type> &m, int row, int col) {

  return element_getter()(m, row, col);
}

/// Function extracting an element from a matrix (non-const)
template <
    typename derived_type,
    std::enable_if_t<std::is_base_of<Eigen::DenseCoeffsBase<
                                         derived_type, Eigen::WriteAccessors>,
                                     Eigen::MatrixBase<derived_type> >::value,
                     bool> = true>
ALGEBRA_HOST_DEVICE inline auto &element(Eigen::MatrixBase<derived_type> &m,
                                         int row, int col) {

  return element_getter()(m, row, col);
}

/// Functor used to extract a block from Eigen matrices
struct block_getter {
  template <int kROWS, int kCOLS, typename matrix_type>
  ALGEBRA_HOST_DEVICE auto operator()(const matrix_type &m, int row,
                                      int col) const {

    return m.template block<kROWS, kCOLS>(row, col);
  }
};  // struct block_getter

}  // namespace algebra::eigen::math
