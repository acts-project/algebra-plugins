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
#include "algebra/storage/vecmem.hpp"

/// @name Operators on @c algebra::vecmem::storage_type
/// @{

using algebra::cmath::operator*;
using algebra::cmath::operator-;
using algebra::cmath::operator+;

/// @}

namespace algebra {
namespace vecmem {

/// @name cmath based transforms on @c algebra::vecmem::storage_type
/// @{

using transform3 = cmath::transform3<std::size_t, vecmem::storage_type, scalar>;
using cartesian2 = cmath::cartesian2<transform3>;
using polar2 = cmath::polar2<transform3>;
using cylindrical2 = cmath::cylindrical2<transform3>;

/// @}

}  // namespace vecmem

namespace getter {

/// @name Getter functions on @c algebra::vecmem::storage_type
/// @{

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::vecmem::transform3
template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS>
ALGEBRA_HOST_DEVICE inline vecmem::storage_type<scalar, SIZE> vector(
    const vecmem::storage_type<vecmem::storage_type<scalar, ROWS>, COLS>& m,
    std::size_t row, std::size_t col) {

  return cmath::vector_getter<std::size_t, vecmem::storage_type, scalar,
                              SIZE>()(m, row, col);
}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vecmem::storage_type
/// @{

using cmath::cross;
using cmath::dot;
using cmath::normalize;

/// @}

}  // namespace vector
}  // namespace algebra
