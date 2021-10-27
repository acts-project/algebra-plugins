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
#include "algebra/math/vc.hpp"
#include "algebra/storage/vc.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

namespace algebra {

/// @name Operators on @c algebra::vc types
/// @{

using cmath::operator*;
using cmath::operator-;
using cmath::operator+;

/// @}

namespace vc {

/// @name Vc based transforms on @c algebra::vc::storage_type
/// @{

// Pull in the definitions needed by the cmath transforms, into this namespace.
using math::perp;
using math::phi;

using transform3 = math::transform3<storage_type, scalar, vector3, point2>;
using cartesian2 = cmath::cartesian2<transform3>;
using polar2 = cmath::polar2<transform3>;
using cylindrical2 = cmath::cylindrical2<transform3>;

/// @}

}  // namespace vc

namespace getter {

/// @name Getter functions on @c algebra::vc types
/// @{

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

using vc::math::eta;
using vc::math::norm;
using vc::math::perp;
using vc::math::phi;
using vc::math::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::array::transform3
template <auto SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline auto vector(const vc::transform3::matrix44& m,
                                       std::size_t row, std::size_t col) {

  assert(row == 0);
  assert(col < 4);
  switch (col) {
    case 0:
      return m.x;
    case 1:
      return m.y;
    case 2:
      return m.z;
    case 3:
      return m.t;
    default:
      return m.x;
  }
}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc types
/// @{

using cmath::dot;
using cmath::normalize;
using vc::math::cross;
using vc::math::dot;
using vc::math::normalize;

/// @}

}  // namespace vector
}  // namespace algebra
