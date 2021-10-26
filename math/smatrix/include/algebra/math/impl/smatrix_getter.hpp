/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
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

/** This method retrieves phi from a vector, vector base with rows >= 2
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N, std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST inline scalar_t phi(
    const ROOT::Math::SVector<scalar_t, N> &v) noexcept {

  return TMath::ATan2(v[1], v[0]);
}

template <typename scalar_t, class A, auto N,
          std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST inline scalar_t phi(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) noexcept {

  return TMath::ATan2(v.apply(1), v.apply(0));
}

/** This method retrieves theta from a vector, vector base with rows >= 3
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N, std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST inline scalar_t theta(
    const ROOT::Math::SVector<scalar_t, N> &v) noexcept {

  return TMath::ATan2(TMath::Sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
}

template <typename scalar_t, class A, auto N,
          std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST inline scalar_t theta(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) noexcept {

  return TMath::ATan2(
      TMath::Sqrt(v.apply(0) * v.apply(0) + v.apply(1) * v.apply(1)),
      v.apply(2));
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST inline scalar_t norm(const ROOT::Math::SVector<scalar_t, N> &v) {

  return TMath::Sqrt(ROOT::Math::Dot(v, v));
}

template <typename scalar_t, class A, auto N>
ALGEBRA_HOST inline scalar_t norm(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) {

  return TMath::Sqrt(ROOT::Math::Dot(v, v));
}

/** This method retrieves the pseudo-rapidity from a vector or vector base with
 *rows >= 3
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N, std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST inline scalar_t eta(
    const ROOT::Math::SVector<scalar_t, N> &v) noexcept {

  return TMath::ATanH(v[2] / norm(v));
}

template <typename scalar_t, class A, auto N,
          std::enable_if_t<N >= 3, bool> = true>
ALGEBRA_HOST inline scalar_t eta(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) noexcept {

  return TMath::ATanH(v.apply(2) / norm(v));
}

/** This method retrieves the perpenticular magnitude of a vector with rows >= 2
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N, std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST inline scalar_t perp(
    const ROOT::Math::SVector<scalar_t, N> &v) noexcept {

  return TMath::Sqrt(v[0] * v[0] + v[1] * v[1]);
}

template <typename scalar_t, class A, auto N,
          std::enable_if_t<N >= 2, bool> = true>
ALGEBRA_HOST inline scalar_t perp(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) noexcept {

  return TMath::Sqrt(v.apply(0) * v.apply(0) + v.apply(1) * v.apply(1));
}

}  // namespace algebra::smatrix::math
