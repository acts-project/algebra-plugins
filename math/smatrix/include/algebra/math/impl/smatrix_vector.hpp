/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// ROOT/Smatrix include(s).
#include <Math/Expression.h>
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <TMath.h>

namespace algebra::smatrix::vector {

template <typename scalar_t>
struct actor {

  using scalar_type = scalar_t;

  /** This method retrieves phi from a vector, vector base with rows >= 2
   *
   * @param v the input vector
   **/
  template <auto N, std::enable_if_t<N >= 2, bool> = true>
  ALGEBRA_HOST inline scalar_type phi(
      const ROOT::Math::SVector<scalar_type, N> &v) noexcept {

    return TMath::ATan2(v[1], v[0]);
  }

  template <class A, auto N, std::enable_if_t<N >= 2, bool> = true>
  ALGEBRA_HOST inline scalar_type phi(
      const ROOT::Math::VecExpr<A, scalar_type, N> &v) noexcept {

    return TMath::ATan2(v.apply(1), v.apply(0));
  }

  /** This method retrieves theta from a vector, vector base with rows >= 3
   *
   * @param v the input vector
   **/
  template <auto N, std::enable_if_t<N >= 3, bool> = true>
  ALGEBRA_HOST inline scalar_type theta(
      const ROOT::Math::SVector<scalar_type, N> &v) noexcept {

    return TMath::ATan2(TMath::Sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
  }

  template <class A, auto N, std::enable_if_t<N >= 3, bool> = true>
  ALGEBRA_HOST inline scalar_type theta(
      const ROOT::Math::VecExpr<A, scalar_type, N> &v) noexcept {

    return TMath::ATan2(
        TMath::Sqrt(v.apply(0) * v.apply(0) + v.apply(1) * v.apply(1)),
        v.apply(2));
  }

  /** This method retrieves the norm of a vector, no dimension restriction
   *
   * @param v the input vector
   **/
  template <auto N>
  ALGEBRA_HOST inline scalar_type norm(
      const ROOT::Math::SVector<scalar_type, N> &v) {

    return TMath::Sqrt(ROOT::Math::Dot(v, v));
  }

  template <class A, auto N>
  ALGEBRA_HOST inline scalar_type norm(
      const ROOT::Math::VecExpr<A, scalar_type, N> &v) {

    return TMath::Sqrt(ROOT::Math::Dot(v, v));
  }

  /** This method retrieves the pseudo-rapidity from a vector or vector base
   *with rows >= 3
   *
   * @param v the input vector
   **/
  template <auto N, std::enable_if_t<N >= 3, bool> = true>
  ALGEBRA_HOST inline scalar_type eta(
      const ROOT::Math::SVector<scalar_type, N> &v) noexcept {

    return TMath::ATanH(v[2] / norm(v));
  }

  template <class A, auto N, std::enable_if_t<N >= 3, bool> = true>
  ALGEBRA_HOST inline scalar_type eta(
      const ROOT::Math::VecExpr<A, scalar_type, N> &v) noexcept {

    return TMath::ATanH(v.apply(2) / norm(v));
  }

  /** This method retrieves the perpenticular magnitude of a vector with rows >=
   *2
   *
   * @param v the input vector
   **/
  template <auto N, std::enable_if_t<N >= 2, bool> = true>
  ALGEBRA_HOST inline scalar_type perp(
      const ROOT::Math::SVector<scalar_type, N> &v) noexcept {

    return TMath::Sqrt(v[0] * v[0] + v[1] * v[1]);
  }

  template <class A, auto N, std::enable_if_t<N >= 2, bool> = true>
  ALGEBRA_HOST inline scalar_type perp(
      const ROOT::Math::VecExpr<A, scalar_type, N> &v) noexcept {

    return TMath::Sqrt(v.apply(0) * v.apply(0) + v.apply(1) * v.apply(1));
  }

  /** Get a normalized version of the input vector
   *
   * @param v the input vector
   **/
  template <auto N>
  ALGEBRA_HOST inline ROOT::Math::SVector<scalar_type, N> normalize(
      const ROOT::Math::SVector<scalar_type, N> &v) {

    return ROOT::Math::Unit(v);
  }
  /** Get a normalized version of the input vector
   *
   * @param v the input vector
   **/
  template <class A, auto N>
  ALGEBRA_HOST inline ROOT::Math::SVector<scalar_type, N> normalize(
      const ROOT::Math::VecExpr<A, scalar_type, N> &v) {

    return ROOT::Math::Unit(v);
  }

  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <auto N>
  ALGEBRA_HOST inline scalar_type dot(
      const ROOT::Math::SVector<scalar_type, N> &a,
      const ROOT::Math::SVector<scalar_type, N> &b) {

    return ROOT::Math::Dot(a, b);
  }
  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <class A, auto N>
  ALGEBRA_HOST inline scalar_type dot(
      const ROOT::Math::SVector<scalar_type, N> &a,
      const ROOT::Math::VecExpr<A, scalar_type, N> &b) {

    return ROOT::Math::Dot(a, b);
  }
  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <class A, auto N>
  ALGEBRA_HOST inline scalar_type dot(
      const ROOT::Math::VecExpr<A, scalar_type, N> &a,
      const ROOT::Math::SVector<scalar_type, N> &b) {

    return ROOT::Math::Dot(a, b);
  }

  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <class A, auto N>
  ALGEBRA_HOST inline scalar_type dot(
      const ROOT::Math::VecExpr<A, scalar_type, N> &a,
      const ROOT::Math::VecExpr<A, scalar_type, N> &b) {

    return ROOT::Math::Dot(a, b);
  }

  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <class A, auto N>
  ALGEBRA_HOST inline scalar_type dot(
      const ROOT::Math::SMatrix<scalar_type, N, 1> &a,
      const ROOT::Math::VecExpr<A, scalar_type, N> &b) {

    return ROOT::Math::Dot(a.Col(0), b);
  }

  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <class A, auto N>
  ALGEBRA_HOST inline scalar_type dot(
      const ROOT::Math::VecExpr<A, scalar_type, N> &a,
      const ROOT::Math::SMatrix<scalar_type, N, 1> &b) {
    return dot(b, a);
  }

  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <auto N>
  ALGEBRA_HOST inline scalar_type dot(
      const ROOT::Math::SMatrix<scalar_type, N, 1> &a,
      const ROOT::Math::SVector<scalar_type, N> &b) {

    return ROOT::Math::Dot(a.Col(0), b);
  }

  /** Dot product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return the scalar dot product value
   **/
  template <auto N>
  ALGEBRA_HOST inline scalar_type dot(
      const ROOT::Math::SVector<scalar_type, N> &a,
      const ROOT::Math::SMatrix<scalar_type, N, 1> &b) {
    return dot(b, a);
  }

  /** Cross product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return a vector (expression) representing the cross product
   **/
  ALGEBRA_HOST inline ROOT::Math::SVector<scalar_type, 3> cross(
      const ROOT::Math::SVector<scalar_type, 3> &a,
      const ROOT::Math::SVector<scalar_type, 3> &b) {

    return ROOT::Math::Cross(a, b);
  }
  /** Cross product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return a vector (expression) representing the cross product
   **/
  template <class A>
  ALGEBRA_HOST inline ROOT::Math::SVector<scalar_type, 3> cross(
      const ROOT::Math::SVector<scalar_type, 3> &a,
      const ROOT::Math::VecExpr<A, scalar_type, 3> &b) {

    return ROOT::Math::Cross(a, b);
  }
  /** Cross product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return a vector (expression) representing the cross product
   **/
  template <class A>
  ALGEBRA_HOST inline ROOT::Math::SVector<scalar_type, 3> cross(
      const ROOT::Math::VecExpr<A, scalar_type, 3> &a,
      const ROOT::Math::SVector<scalar_type, 3> &b) {

    return ROOT::Math::Cross(a, b);
  }
  /** Cross product between two input vectors
   *
   * @param a the first input vector
   * @param b the second input vector
   *
   * @return a vector (expression) representing the cross product
   **/
  template <class A>
  ALGEBRA_HOST inline ROOT::Math::SVector<scalar_type, 3> cross(
      const ROOT::Math::VecExpr<A, scalar_type, 3> &a,
      const ROOT::Math::VecExpr<A, scalar_type, 3> &b) {

    return ROOT::Math::Cross(a, b);
  }

  /** Cross product between vector3 and matrix<3,1>
   *
   * @param a the first input vector
   * @param b the second input matrix<3,1>
   *
   * @return a vector (expression) representing the cross product
   **/
  ALGEBRA_HOST inline ROOT::Math::SVector<scalar_type, 3> cross(
      const ROOT::Math::SVector<scalar_type, 3> &a,
      const ROOT::Math::SMatrix<scalar_type, 3, 1> &b) {

    return ROOT::Math::Cross(a, b.Col(0));
  }

  /** Cross product between matrix<3,1> and vector3
   *
   * @param a the second input matrix<3,1>
   * @param b the first input vector
   *
   * @return a vector (expression) representing the cross product
   **/
  ALGEBRA_HOST inline ROOT::Math::SVector<scalar_type, 3> cross(
      const ROOT::Math::SMatrix<scalar_type, 3, 1> &a,
      const ROOT::Math::SVector<scalar_type, 3> &b) {

    return ROOT::Math::Cross(a.Col(0), b);
  }

  /** Cross product between two matrix<3,1>
   *
   * @param a the second input matrix<3,1>
   * @param b the first input matrix<3,1>
   *
   * @return a vector (expression) representing the cross product
   **/
  ALGEBRA_HOST inline ROOT::Math::SVector<scalar_type, 3> cross(
      const ROOT::Math::SMatrix<scalar_type, 3, 1> &a,
      const ROOT::Math::SMatrix<scalar_type, 3, 1> &b) {

    return ROOT::Math::Cross(a.Col(0), b.Col(0));
  }
};

}  // namespace algebra::smatrix::vector
