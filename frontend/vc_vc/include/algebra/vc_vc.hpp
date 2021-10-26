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
#include <type_traits>

namespace algebra {

/// @name Operators on @c algebra::vc::storage_type
/// @{

using cmath::operator*;
using cmath::operator-;
using cmath::operator+;

/// @}

namespace vc {

/// @name Vc based transforms on @c algebra::vc::storage_type
/// @{

using transform3 = math::transform3<storage_type, scalar, vector3, point2>;
using cartesian2 = math::cartesian2<transform3>;
using polar2 = math::polar2<storage_type, transform3>;
using cylindrical2 = math::cylindrical2<storage_type, transform3>;

/// @}

}  // namespace vc

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
template <auto SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline vc::storage_type<scalar, SIZE> vector(
    const vc::transform3::matrix44& m, std::size_t row, std::size_t col) {

  vc::storage_type<scalar, SIZE> result;
  const Vc::SimdArray<scalar, 4>* a = nullptr;
  switch (col) {
    case 0:
      a = &(m.x);
      break;
    case 1:
      a = &(m.y);
      break;
    case 2:
      a = &(m.z);
      break;
    case 3:
      a = &(m.t);
      break;
    default:
      return result;
  }
  for (std::size_t irow = row; irow < row + SIZE; ++irow) {
    result[irow - row] = (*a)[irow];
  }
  return result;
}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc::storage_type
/// @{

using cmath::dot;
using cmath::normalize;
using vc::math::cross;
using vc::math::dot;
using vc::math::normalize;

/// @}

}  // namespace vector
}  // namespace algebra
