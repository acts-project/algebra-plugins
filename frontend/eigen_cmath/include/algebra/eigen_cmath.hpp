/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/eigen.hpp"
#include "algebra/storage/eigen.hpp"

// Eigen include(s).
#include <Eigen/Geometry>

// System include(s).
#include <type_traits>

namespace algebra {
namespace eigen {

/// Functor used to access elements of Eigen matrices
struct element_getter {
  /// Get non-const access to a matrix element
  template <
      typename derived_type,
      std::enable_if_t<std::is_base_of<Eigen::DenseCoeffsBase<
                                           derived_type, Eigen::WriteAccessors>,
                                       Eigen::MatrixBase<derived_type> >::value,
                       bool> = true>
  ALGEBRA_HOST_DEVICE inline auto& operator()(
      Eigen::MatrixBase<derived_type>& m, unsigned int row,
      unsigned int col) const {

    return m(row, col);
  }
  /// Get const access to a matrix element
  template <typename derived_type>
  ALGEBRA_HOST_DEVICE inline auto operator()(
      const Eigen::MatrixBase<derived_type>& m, unsigned int row,
      unsigned int col) const {

    return m(row, col);
  }
};  // struct element_getter

/// Functor used to extract a block from Eigen matrices
struct block_getter {
  template <std::size_t kROWS, std::size_t kCOLS, typename matrix_type>
  ALGEBRA_HOST_DEVICE auto operator()(const matrix_type& m, std::size_t row,
                                      std::size_t col) const {

    return m.template block<kROWS, kCOLS>(row, col);
  }
};  // struct block_getter

/// @name cmath based transforms on @c algebra::eigen::storage_type
/// @{

template <typename T>
using transform3 = cmath::transform3<
    std::size_t, eigen::storage_type, T,
    typename Eigen::Transform<T, 3, Eigen::Affine>::MatrixType,
    algebra::eigen::element_getter, algebra::eigen::block_getter>;
template <typename T>
using cartesian2 = cmath::cartesian2<transform3<T> >;
template <typename T>
using polar2 = cmath::polar2<transform3<T> >;
template <typename T>
using cylindrical2 = cmath::cylindrical2<transform3<T> >;

/// @}

}  // namespace eigen

namespace getter {

/// @name Getter functions on @c algebra::eigen::storage_type
/// @{

using eigen::math::eta;
using eigen::math::norm;
using eigen::math::perp;
using eigen::math::phi;
using eigen::math::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::eigen::transform3
template <unsigned int SIZE, typename derived_type>
ALGEBRA_HOST_DEVICE inline auto vector(const Eigen::MatrixBase<derived_type>& m,
                                       std::size_t row, std::size_t col) {

  return m.template block<SIZE, 1>(row, col);
}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::eigen::storage_type
/// @{

using eigen::math::cross;
using eigen::math::dot;
using eigen::math::normalize;

/// @}

}  // namespace vector
}  // namespace algebra
