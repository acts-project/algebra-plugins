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
#include "algebra/storage/smatrix.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

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

}  // namespace algebra
