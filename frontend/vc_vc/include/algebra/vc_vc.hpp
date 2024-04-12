/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/vc.hpp"
#include "algebra/storage/vc.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

namespace algebra {
namespace vc {

/// @name Vc based transforms on @c algebra::vc::storage_type
/// @{

template <typename T>
using transform3 = math::transform3<storage_type, T, vector3<T>, point2<T>>;

/// @}

}  // namespace vc

namespace getter {

/// @name Getter functions on @c algebra::vc types
/// @{

using vc::math::eta;
using vc::math::norm;
using vc::math::perp;
using vc::math::phi;
using vc::math::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::vc::transform3<float>
template <std::size_t SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline vc::vector3<float> vector(
    const vc::transform3<float>::matrix44& m,
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
      return {m.x[0], m.x[1], m.x[2]};
    case 1:
      return {m.y[0], m.y[1], m.y[2]};
    case 2:
      return {m.z[0], m.z[1], m.z[2]};
    case 3:
      return {m.t[0], m.t[1], m.t[2]};
    default:
      return {m.x[0], m.x[1], m.x[2]};
  }
}

/// Function extracting a slice from the matrix used by
/// @c algebra::vc::transform3<double>
template <std::size_t SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline vc::vector3<double> vector(
    const vc::transform3<double>::matrix44& m,
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
      return {m.x[0], m.x[1], m.x[2]};
    case 1:
      return {m.y[0], m.y[1], m.y[2]};
    case 2:
      return {m.z[0], m.z[1], m.z[2]};
    case 3:
      return {m.t[0], m.t[1], m.t[2]};
    default:
      return {m.x[0], m.x[1], m.x[2]};
  }
}

/// @name Getter functions on @c algebra::vc::matrix_type
/// @{

using cmath::element;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc types
/// @{

using vc::math::cross;
using vc::math::dot;
using vc::math::normalize;

/// @}

}  // namespace vector

namespace matrix {

using size_type = vc::size_type;

template <typename T, size_type N>
using array_type = vc::storage_type<T, N>;

template <typename T, size_type ROWS, size_type COLS>
using matrix_type = vc::matrix_type<T, ROWS, COLS>;

}  // namespace matrix
}  // namespace algebra
