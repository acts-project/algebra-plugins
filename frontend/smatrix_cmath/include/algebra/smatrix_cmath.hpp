/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/smatrix.hpp"
#include "algebra/storage/smatrix.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

namespace algebra {
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

  return m.template SubCol<smatrix::storage_type<scalar_t, SIZE>>(col, row);
}

/// @name Getter functions on @c algebra::smatrix::matrix_type
/// @{

using smatrix::math::element;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::cross;
using smatrix::math::dot;
using smatrix::math::normalize;

/// @}

}  // namespace vector

namespace matrix {

using smatrix::math::block;
using smatrix::math::determinant;
using smatrix::math::identity;
using smatrix::math::inverse;
using smatrix::math::set_block;
using smatrix::math::set_identity;
using smatrix::math::set_zero;
using smatrix::math::transpose;
using smatrix::math::zero;

}  // namespace matrix

namespace smatrix {

template <typename scalar_t>
using element_getter = smatrix::math::element_getter<scalar_t>;

template <typename scalar_t>
using block_getter = smatrix::math::block_getter<scalar_t>;

/// @name cmath based transforms on @c algebra::smatrix
/// @{

template <typename T>
using transform3 =
    cmath::transform3<smatrix::size_type, T, smatrix::matrix_type,
                      smatrix::storage_type, element_getter<T>,
                      block_getter<T>>;

/// @}

}  // namespace smatrix

}  // namespace algebra
