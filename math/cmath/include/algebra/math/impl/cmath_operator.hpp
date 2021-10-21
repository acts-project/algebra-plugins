/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// System include(s).
#include <cstddef>

namespace algebra::cmath {

/// @name Operators on 2-element arrays
/// @{

template <template <typename, std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator*(
    const array_t<scalar_t, 2> &a, scalar_t s) {

  return {a[0] * s, a[1] * s};
}

template <template <typename, std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator*(
    scalar_t s, const array_t<scalar_t, 2> &a) {

  return {s * a[0], s * a[1]};
}

template <template <typename, std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator-(
    const array_t<scalar_t, 2> &a, const array_t<scalar_t, 2> &b) {

  return {a[0] - b[0], a[1] - b[1]};
}

template <template <typename, std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 2> operator+(
    const array_t<scalar_t, 2> &a, const array_t<scalar_t, 2> &b) {

  return {a[0] + b[0], a[1] + b[1]};
}

/// @}

/// @name Operators on 3-element arrays
/// @{

template <template <typename, std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator*(
    const array_t<scalar_t, 3> &a, scalar_t s) {

  return {a[0] * s, a[1] * s, a[2] * s};
}

template <template <typename, std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator*(
    scalar_t s, const array_t<scalar_t, 3> &a) {

  return {s * a[0], s * a[1], s * a[2]};
}

template <template <typename, std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator-(
    const array_t<scalar_t, 3> &a, const array_t<scalar_t, 3> &b) {

  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

template <template <typename, std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE inline array_t<scalar_t, 3> operator+(
    const array_t<scalar_t, 3> &a, const array_t<scalar_t, 3> &b) {

  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

/// @}

}  // namespace algebra::cmath
