/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/vc_aos.hpp"
#include "algebra/storage/vc_aos.hpp"

/// @name Operators on @c algebra::vc_aos types
/// @{

using algebra::cmath::operator*;
using algebra::cmath::operator-;
using algebra::cmath::operator+;

/// @}

namespace algebra {
namespace getter {

/// @name Getter functions on @c algebra::vc_aos types
/// @{

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

using vc_aos::math::eta;
using vc_aos::math::norm;
using vc_aos::math::perp;
using vc_aos::math::phi;
using vc_aos::math::theta;

/// @|

/// Function extracting a slice from the matrix used by
/// @c algebra::vc_aos::transform3
template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS,
          typename scalar_t>
ALGEBRA_HOST_DEVICE inline Vc::array<scalar_t, SIZE> vector(
    const vc_aos::matrix_type<scalar_t, ROWS, COLS>& m, std::size_t row,
    std::size_t col) {

  return cmath::vector_getter<std::size_t, Vc::array, scalar_t, SIZE,
                              Vc::array<scalar_t, SIZE>>()(m, row, col);
}

/// @name Getter functions on @c algebra::vc_aos::matrix_type
/// @{

using cmath::element;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc_aos::storage_type
/// @{

using cmath::cross;
using cmath::dot;
using cmath::normalize;

using vc_aos::math::cross;
using vc_aos::math::dot;
using vc_aos::math::normalize;

/// @}

}  // namespace vector

namespace matrix {

using cmath::block;
using cmath::determinant;
using cmath::identity;
using cmath::inverse;
using cmath::set_block;
using cmath::set_identity;
using cmath::set_zero;
using cmath::transpose;
using cmath::zero;

}  // namespace matrix

namespace vc_aos {

template <typename scalar_t>
using element_getter = cmath::element_getter<size_type, Vc::array, scalar_t>;

/// @name cmath based transforms on @c algebra::vc

template <typename T>
using transform3 =
    cmath::transform3<vc_aos::size_type, T, vc_aos::matrix_type, Vc::array,
                      element_getter<T>, cmath::block_getter>;

/// @}

}  // namespace vc_aos

}  // namespace algebra
