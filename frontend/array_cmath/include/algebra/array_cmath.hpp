/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/generic.hpp"
#include "algebra/print.hpp"
#include "algebra/storage/array.hpp"

/// @name Operators on @c algebra::array::storage_type
/// @{

using algebra::cmath::operator*;
using algebra::cmath::operator-;
using algebra::cmath::operator+;

/// Print the linear algebra types of this backend
using algebra::operator<<;

/// @}

namespace algebra {

namespace getter {

/// @name Getter functions on @c algebra::array::storage_type
/// @{

using cmath::storage::block;
using cmath::storage::element;
using cmath::storage::set_block;
using cmath::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::array::storage_type
/// @{

// array specific implementations
using cmath::dot;
using cmath::normalize;

// generic implementations
using cmath::cross;
using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

/// @}

}  // namespace vector

// Use special algorithms for 4 dimensional matrices
namespace generic {

// Determinant algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct determinant_selector<4, array::matrix_type<T, ROWS, COLS>> {
  using type =
      matrix::determinant::hard_coded<array::matrix_type<T, ROWS, COLS>,
                                      array::element_getter>;
};

// Inversion algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct inversion_selector<4, array::matrix_type<T, ROWS, COLS>> {
  using type = matrix::inverse::hard_coded<array::matrix_type<T, ROWS, COLS>,
                                           array::element_getter>;
};

}  // namespace generic

namespace matrix {

/// @name Matrix functions on @c algebra::array::storage_type
/// @{

using cmath::identity;
using cmath::set_identity;
using cmath::set_zero;
using cmath::zero;

// Uses generic implementation in the background
using cmath::determinant;
using cmath::inverse;
using cmath::transpose;

/// @}

}  // namespace matrix

namespace array {

/// @name cmath based transforms on @c algebra::array
/// @{

template <concepts::scalar T>
using transform3 =
    generic::math::transform3<array::size_type, T, array::matrix_type,
                              array::storage_type>;

/// @}

}  // namespace array

}  // namespace algebra
