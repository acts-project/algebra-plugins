/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/vc_aos.hpp"
#include "algebra/storage/vc_aos.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

namespace algebra {
namespace vc_aos {

/// @name Vc based transforms on @c algebra::vc_aos::storage_type
/// @{

template <typename T>
using transform3 = math::transform3<storage_type, T, vector3<T>, point2<T>>;

/// @}

}  // namespace vc_aos

namespace getter {

/// @name Getter functions on @c algebra::vc_aos types
/// @{

using vc_aos::math::eta;
using vc_aos::math::norm;
using vc_aos::math::perp;
using vc_aos::math::phi;
using vc_aos::math::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::vc_aos::transform3<float>
template <std::size_t SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline vc_aos::vector3<float> vector(
    const vc_aos::transform3<float>::matrix44& m,
    std::size_t
#ifndef NDEBUG
        row
#endif  // not NDEBUG
    ,
    std::size_t col) {

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

/// Function extracting a slice from the matrix used by
/// @c algebra::vc_aos::transform3<double>
template <std::size_t SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline vc_aos::vector3<double> vector(
    const vc_aos::transform3<double>::matrix44& m,
    std::size_t
#ifndef NDEBUG
        row
#endif  // not NDEBUG
    ,
    std::size_t col) {

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

/// @name Getter functions on @c algebra::vc_aos::matrix_type
/// @{

using cmath::element;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc_aos types
/// @{

using vc_aos::math::cross;
using vc_aos::math::dot;
using vc_aos::math::normalize;

/// @}

}  // namespace vector

}  // namespace algebra
