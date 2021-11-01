/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

namespace algebra::cmath {

/// @name Operators on 2-element arrays
/// @{

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator*(
    const array_t<scalar_t, 2> &a, float s) {

  return {a[0] * static_cast<scalar_t>(s), a[1] * static_cast<scalar_t>(s)};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator*(
    const array_t<scalar_t, 2> &a, double s) {

  return {a[0] * static_cast<scalar_t>(s), a[1] * static_cast<scalar_t>(s)};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator*(
    float s, const array_t<scalar_t, 2> &a) {

  return {static_cast<scalar_t>(s) * a[0], static_cast<scalar_t>(s) * a[1]};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator*(
    double s, const array_t<scalar_t, 2> &a) {

  return {static_cast<scalar_t>(s) * a[0], static_cast<scalar_t>(s) * a[1]};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator-(
    const array_t<scalar_t, 2> &a, const array_t<scalar_t, 2> &b) {

  return {a[0] - b[0], a[1] - b[1]};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator+(
    const array_t<scalar_t, 2> &a, const array_t<scalar_t, 2> &b) {

  return {a[0] + b[0], a[1] + b[1]};
}

/// @}

/// @name Operators on 3-element arrays
/// @{

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator*(
    const array_t<scalar_t, 3> &a, float s) {

  return {a[0] * static_cast<scalar_t>(s), a[1] * static_cast<scalar_t>(s),
          a[2] * static_cast<scalar_t>(s)};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator*(
    const array_t<scalar_t, 3> &a, double s) {

  return {a[0] * static_cast<scalar_t>(s), a[1] * static_cast<scalar_t>(s),
          a[2] * static_cast<scalar_t>(s)};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator*(
    float s, const array_t<scalar_t, 3> &a) {

  return {static_cast<scalar_t>(s) * a[0], static_cast<scalar_t>(s) * a[1],
          static_cast<scalar_t>(s) * a[2]};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator*(
    double s, const array_t<scalar_t, 3> &a) {

  return {static_cast<scalar_t>(s) * a[0], static_cast<scalar_t>(s) * a[1],
          static_cast<scalar_t>(s) * a[2]};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator-(
    const array_t<scalar_t, 3> &a, const array_t<scalar_t, 3> &b) {

  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

template <template <typename, auto> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator+(
    const array_t<scalar_t, 3> &a, const array_t<scalar_t, 3> &b) {

  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

/// @}

}  // namespace algebra::cmath
