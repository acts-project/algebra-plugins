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
using transform3 = math::transform3<storage_type, T>;

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

/// @name Getter functions on @c algebra::vc_aos::matrix_type
/// @{

using cmath::element;

/// Function extracting a slice from matrix44 - const
template <std::size_t SIZE, template <typename, std::size_t> class array_t,
          typename value_t, std::size_t N,
          std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline const auto& vector(
    const storage::matrix44<array_t, value_t, N>& m,
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
#ifndef _MSC_VER
      __builtin_unreachable();
#else
      return m.x;
#endif
  }
}

/// Function extracting a slice from matrix44 - non-const
template <std::size_t SIZE, template <typename, std::size_t> class array_t,
          typename value_t, std::size_t N,
          std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline auto& vector(
    storage::matrix44<array_t, value_t, N>& m,
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
#ifndef _MSC_VER
      __builtin_unreachable();
#else
      return m.x;
#endif
  }
}

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
