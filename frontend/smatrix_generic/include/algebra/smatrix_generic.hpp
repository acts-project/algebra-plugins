/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/generic.hpp"
#include "algebra/math/smatrix.hpp"
#include "algebra/print.hpp"
#include "algebra/storage/smatrix.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

/// Print the linear algebra types of this backend
using algebra::operator<<;

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::smatrix::matrix_type
/// @{

using smatrix::storage::block;
using smatrix::storage::element;
using smatrix::storage::set_block;
using smatrix::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::smatrix::storage_type
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

/// @name Matrix functions on @c algebra::smatrix::storage_type
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

namespace smatrix {

/// @name generic based transforms on @c algebra::smatrix
/// @{

template <concepts::scalar T>
using transform3 =
    generic::math::transform3<smatrix::size_type, T, smatrix::matrix_type,
                              smatrix::storage_type>;

/// @}

}  // namespace smatrix

namespace plugin {

/// Define the plugin types
/// @{
template <typename V>
struct smatrix_generic {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using size_type = algebra::smatrix::size_type;
    using transform3D = algebra::smatrix::transform3<value_type>;
    using point2D = algebra::smatrix::point2<value_type>;
    using point3D = algebra::smatrix::point3<value_type>;
    using vector3D = algebra::smatrix::vector3<value_type>;

    template <std::size_t ROWS, std::size_t COLS>
    using matrix = algebra::smatrix::matrix_type<value_type, ROWS, COLS>;
};
/// @}

} // namespace plugin

}  // namespace algebra
