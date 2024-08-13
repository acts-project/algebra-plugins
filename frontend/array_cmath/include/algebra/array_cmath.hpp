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
#include "algebra/storage/array.hpp"

/// @name Operators on @c algebra::array::storage_type
/// @{

using algebra::cmath::operator*;
using algebra::cmath::operator-;
using algebra::cmath::operator+;

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
using generic::math::cross;
using generic::math::eta;
using generic::math::norm;
using generic::math::perp;
using generic::math::phi;
using generic::math::theta;

/// @}

}  // namespace vector

// Use special algorithms for 4 dimensional matrices by partial template
// specilization
// @see "algebra/math/algorithms/utils/algorithm_finder.hpp"
namespace generic {

// Determinant algorithms
template <typename T, auto ROWS, auto COLS>
struct determinant_selector<4, array::matrix_type<T, ROWS, COLS>> {
  using type =
      matrix::determinant::hard_coded<array::matrix_type<T, ROWS, COLS>>;
};

// Inversion algorithms
template <typename T, auto ROWS, auto COLS>
struct inversion_selector<4, array::matrix_type<T, ROWS, COLS>> {
  using type = matrix::inverse::hard_coded<array::matrix_type<T, ROWS, COLS>>;
};

}  // namespace generic

namespace matrix {

/// @name Matrix functions on @c algebra::array::storage_type
/// @{

using cmath::identity;
using cmath::set_identity;
using cmath::set_zero;
using cmath::zero;

using generic::math::determinant;
using generic::math::inverse;
using generic::math::transpose;

/// @}

}  // namespace matrix

namespace array {

/// @name cmath based transforms on @c algebra::array
/// @{

template <typename T>
using transform3 =
    generic::math::transform3<array::size_type, T, array::matrix_type,
                              array::storage_type>;

/// @}

}  // namespace array

}  // namespace algebra
