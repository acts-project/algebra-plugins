/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/smatrix.hpp"
#include "algebra/print.hpp"
#include "algebra/storage/smatrix.hpp"

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

using smatrix::math::cross;
using smatrix::math::dot;
using smatrix::math::eta;
using smatrix::math::norm;
using smatrix::math::normalize;
using smatrix::math::perp;
using smatrix::math::phi;
using smatrix::math::theta;

/// @}

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::determinant;
using smatrix::math::identity;
using smatrix::math::inverse;
using smatrix::math::set_identity;
using smatrix::math::set_zero;
using smatrix::math::transpose;
using smatrix::math::zero;

/// @}

}  // namespace matrix

namespace smatrix {

/// @name SMatrix based transforms on @c algebra::smatrix::storage_type
/// @{

template <typename T>
using transform3 = math::transform3<T>;

/// @}

}  // namespace smatrix

}  // namespace algebra
