/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/smatrix.hpp"
#include "algebra/storage/smatrix.hpp"

namespace algebra {
namespace smatrix {

/// @name SMatrix based transforms on @c algebra::smatrix::storage_type
/// @{

template <typename T>
using transform3 = math::transform3<T>;
template <typename T>
using cartesian2 = math::cartesian2<transform3<T> >;
template <typename T>
using polar2 = math::polar2<transform3<T> >;
template <typename T>
using cylindrical2 = math::cylindrical2<transform3<T> >;

/// @}

}  // namespace smatrix

namespace getter {

/// @name Getter functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::eta;
using smatrix::math::norm;
using smatrix::math::perp;
using smatrix::math::phi;
using smatrix::math::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::smatrix::transform3
template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS,
          typename scalar_t>
ALGEBRA_HOST_DEVICE inline auto vector(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS>& m, unsigned int row,
    unsigned int col) {

  return m.template SubCol<smatrix::storage_type<scalar_t, SIZE> >(col, row);
}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::cross;
using smatrix::math::dot;
using smatrix::math::normalize;

/// @}

}  // namespace vector
}  // namespace algebra
