/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/coordinates/coordinates.hpp"
#include "algebra/math/smatrix.hpp"
#include "algebra/storage/smatrix.hpp"

namespace algebra {

namespace getter {

/// Function extracting a slice from the matrix used by
/// @c algebra::smatrix::transform3
template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS,
          typename scalar_t>
ALGEBRA_HOST_DEVICE inline auto vector(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS>& m, unsigned int row,
    unsigned int col) {

  return m.template SubCol<smatrix::storage_type<scalar_t, SIZE>>(col, row);
}

/// @name Getter functions on @c algebra::smatrix::matrix_type
/// @{

using smatrix::math::element;

/// @}

}  // namespace getter

namespace vector {

template <typename scalar_t>
using actor = smatrix::vector::actor<scalar_t>;

}  // namespace vector

namespace matrix {

template <typename scalar_t>
using actor = smatrix::matrix::actor<scalar_t>;

}  // namespace matrix

namespace smatrix {

/// @name SMatrix based transforms on @c algebra::smatrix::storage_type
/// @{

template <typename T>
using transform3 = math::transform3<T, algebra::vector::actor<T>>;
template <typename T>
using cartesian2 = cartesian2<transform3<T>>;
template <typename T>
using polar2 = polar2<transform3<T>>;
template <typename T>
using cylindrical2 = cylindrical2<transform3<T>>;

/// @}

}  // namespace smatrix

}  // namespace algebra
