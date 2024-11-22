/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/eigen.hpp"
#include "algebra/math/generic.hpp"
#include "algebra/print.hpp"
#include "algebra/storage/eigen.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Eigen/Geometry>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

/// Print the linear algebra types of this backend
using algebra::operator<<;

// System include(s).
#include <type_traits>

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::eigen::matrix_type
/// @{

using eigen::storage::block;
using eigen::storage::element;
using eigen::storage::set_block;
using eigen::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::eigen::storage_type
/// @{

using generic::math::cross;
using generic::math::dot;
using generic::math::eta;
using generic::math::norm;
using generic::math::normalize;
using generic::math::perp;
using generic::math::phi;
using generic::math::theta;

/// @}

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::eigen::storage_type
/// @{

using generic::math::determinant;
using generic::math::identity;
using generic::math::inverse;
using generic::math::set_identity;
using generic::math::set_zero;
using generic::math::transpose;
using generic::math::zero;

/// @}

}  // namespace matrix

namespace eigen {

/// @name generic based transforms on @c algebra::eigen
/// @{

template <concepts::scalar T>
using transform3 =
    generic::math::transform3<eigen::size_type, T, eigen::matrix_type,
                              eigen::storage_type>;

/// @}

}  // namespace eigen

namespace plugin {

/// Define the plugin types
/// @{
template <typename V>
struct eigen_generic {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using size_type = algebra::eigen::size_type;
    using transform3D = algebra::eigen::transform3<value_type>;
    using point2D = algebra::eigen::point2<value_type>;
    using point3D = algebra::eigen::point3<value_type>;
    using vector3D = algebra::eigen::vector3<value_type>;

    template <std::size_t ROWS, std::size_t COLS>
    using matrix = algebra::eigen::matrix_type<value_type, ROWS, COLS>;
};
/// @}

} // namespace plugin

}  // namespace algebra
