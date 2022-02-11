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
namespace smatrix {

/// @name cmath based transforms on @c algebra::smatrix::storage_type
/// @{

template <typename T>
using transform3 =
    cmath::transform3<unsigned int, smatrix::storage_type, T,
                      ROOT::Math::SMatrix<T, 4, 4>, math::element_getter<T>,
                      math::block_getter<T> >;
template <typename T>
using cartesian2 = cmath::cartesian2<transform3<T> >;
template <typename T>
using polar2 = cmath::polar2<transform3<T> >;
template <typename T>
using cylindrical2 = cmath::cylindrical2<transform3<T> >;

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

/// Function extracting an element from a matrix (const)
template <typename scalar_t, unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline scalar_t element(
    const smatrix::matrix_type<scalar_t, ROWS, COLS>& m, std::size_t row,
    std::size_t col) {

  return smatrix::math::element_getter<scalar_t>()(m, row, col);
}

/// Function extracting an element from a matrix (non-const)
template <typename scalar_t, unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline scalar_t& element(
    smatrix::matrix_type<scalar_t, ROWS, COLS>& m, std::size_t row,
    std::size_t col) {

  return smatrix::math::element_getter<scalar_t>()(m, row, col);
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
