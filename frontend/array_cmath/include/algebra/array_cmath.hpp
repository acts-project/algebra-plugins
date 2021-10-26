/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"
#include "algebra/math/cmath.hpp"
#include "algebra/storage/array.hpp"

namespace algebra {

/// @name Operators on @c algebra::array::storage_type
/// @{

using cmath::operator*;
using cmath::operator-;
using cmath::operator+;

/// @}

namespace array {

/// @name cmath based transforms on @c algebra::array::storage_type
/// @{

using transform3 = cmath::transform3<array::storage_type, scalar>;
using cartesian2 = cmath::cartesian2<transform3>;
using polar2 = cmath::polar2<transform3>;
using cylindrical2 = cmath::cylindrical2<transform3>;

/// @}

}  // namespace array

namespace getter {

/// @name Getter functions on @c algebra::array::storage_type
/// @{

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::array::transform3
template <auto SIZE, auto ROWS, auto COLS>
ALGEBRA_HOST_DEVICE inline array::storage_type<scalar, SIZE> vector(
    const array::storage_type<array::storage_type<scalar, ROWS>, COLS>& m,
    std::size_t row, std::size_t col) {

  return cmath::vector_getter<array::storage_type, scalar, SIZE>()(m, row, col);
}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::array::storage_type
/// @{

using cmath::cross;
using cmath::dot;
using cmath::normalize;

/// @}

}  // namespace vector
}  // namespace algebra
