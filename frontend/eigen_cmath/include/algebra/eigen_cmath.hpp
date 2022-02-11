/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/eigen.hpp"
#include "algebra/storage/eigen.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Eigen/Geometry>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s).
#include <type_traits>

namespace algebra {
namespace eigen {

/// @name cmath based transforms on @c algebra::eigen::storage_type
/// @{

template <typename T>
using transform3 = cmath::transform3<
    std::size_t, eigen::storage_type, T,
    typename Eigen::Transform<T, 3, Eigen::Affine>::MatrixType,
    math::element_getter, math::block_getter>;
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

/// Function extracting an element from a matrix (const)
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline auto element(
    const Eigen::MatrixBase<derived_type>& m, std::size_t row,
    std::size_t col) {

  return eigen::math::element_getter()(m, row, col);
}

/// Function extracting an element from a matrix (non-const)
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline auto& element(Eigen::MatrixBase<derived_type>& m,
                                         std::size_t row, std::size_t col) {

  return eigen::math::element_getter()(m, row, col);
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
