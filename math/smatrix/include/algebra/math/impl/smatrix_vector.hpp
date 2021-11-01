/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// ROOT/Smatrix include(s).
#include <Math/Expression.h>
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <TMath.h>

namespace algebra::smatrix::math {

/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST inline ROOT::Math::SVector<scalar_t, N> normalize(
    const ROOT::Math::SVector<scalar_t, N> &v) {

  return ROOT::Math::Unit(v);
}
/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
template <typename scalar_t, class A, auto N>
ALGEBRA_HOST inline ROOT::Math::SVector<scalar_t, N> normalize(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) {

  return ROOT::Math::Unit(v);
}

/** Dot product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST inline scalar_t dot(const ROOT::Math::SVector<scalar_t, N> &a,
                                 const ROOT::Math::SVector<scalar_t, N> &b) {

  return ROOT::Math::Dot(a, b);
}
/** Dot product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template <typename scalar_t, class A, auto N>
ALGEBRA_HOST inline scalar_t dot(const ROOT::Math::SVector<scalar_t, N> &a,
                                 const ROOT::Math::VecExpr<A, scalar_t, N> &b) {

  return ROOT::Math::Dot(a, b);
}
/** Dot product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template <typename scalar_t, class A, auto N>
ALGEBRA_HOST inline scalar_t dot(const ROOT::Math::VecExpr<A, scalar_t, N> &a,
                                 const ROOT::Math::SVector<scalar_t, N> &b) {

  return ROOT::Math::Dot(a, b);
}
/** Dot product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template <typename scalar_t, class A, auto N>
ALGEBRA_HOST inline scalar_t dot(const ROOT::Math::VecExpr<A, scalar_t, N> &a,
                                 const ROOT::Math::VecExpr<A, scalar_t, N> &b) {

  return ROOT::Math::Dot(a, b);
}

/** Cross product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector (expression) representing the cross product
 **/
template <typename scalar_t>
ALGEBRA_HOST inline ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::SVector<scalar_t, 3> &a,
    const ROOT::Math::SVector<scalar_t, 3> &b) {

  return ROOT::Math::Cross(a, b);
}
/** Cross product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector (expression) representing the cross product
 **/
template <typename scalar_t, class A>
ALGEBRA_HOST inline ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::SVector<scalar_t, 3> &a,
    const ROOT::Math::VecExpr<A, scalar_t, 3> &b) {

  return ROOT::Math::Cross(a, b);
}
/** Cross product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector (expression) representing the cross product
 **/
template <typename scalar_t, class A>
ALGEBRA_HOST inline ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::VecExpr<A, scalar_t, 3> &a,
    const ROOT::Math::SVector<scalar_t, 3> &b) {

  return ROOT::Math::Cross(a, b);
}
/** Cross product between two input vectors
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector (expression) representing the cross product
 **/
template <typename scalar_t, class A>
ALGEBRA_HOST inline ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::VecExpr<A, scalar_t, 3> &a,
    const ROOT::Math::VecExpr<A, scalar_t, 3> &b) {

  return ROOT::Math::Cross(a, b);
}

}  // namespace algebra::smatrix::math
