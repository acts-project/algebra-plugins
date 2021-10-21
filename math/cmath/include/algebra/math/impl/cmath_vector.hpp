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
#include <cmath>

namespace algebra {
namespace vector {

/** Cross product between two input vectors - 3 Dim
 *
 * @tparam derived_type_lhs is the first matrix (epresseion) template
 * @tparam derived_type_rhs is the second matrix (epresseion) template
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector (expression) representing the cross product
 **/
template<template <typename,std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE
inline array_t<scalar_t, 3> cross(const array_t<scalar_t, 3> &a,
                                  const array_t<scalar_t, 3> &b) {

  return {a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0],
          a[0] * b[1] - b[0] * a[1]};
}

/** Dot product between two input vectors - 2 Dim
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template<template <typename,std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE
inline scalar_t dot(const array_t<scalar_t, 2> &a,
                    const array_t<scalar_t, 2> &b) {

  return a[0] * b[0] + a[1] * b[1];
}

/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
template<template <typename,std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE
inline array_t<scalar_t, 2> normalize(const array_t<scalar_t, 2> &v) {

  scalar_t oon = 1. / std::sqrt(dot(v, v));
  return {v[0] * oon, v[1] * oon};
}

/** Dot product between two input vectors - 3 Dim
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template<template <typename,std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE
inline scalar_t dot(const array_t<scalar_t, 3> &a,
                    const array_t<scalar_t, 3> &b) {

  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
template<template <typename,std::size_t> class array_t, typename scalar_t>
ALGEBRA_HOST_DEVICE
inline array_t<scalar_t, 3> normalize(const array_t<scalar_t, 3> &v) {

  scalar oon = 1. / std::sqrt(dot(v, v));
  return {v[0] * oon, v[1] * oon, v[2] * oon};
}

}  // namespace vector
}  // namespace algebra
