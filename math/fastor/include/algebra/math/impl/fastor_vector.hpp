/** Algebra plugins, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

namespace algebra::fastor::math {

/// This method retrieves phi from a vector @param v
template <concepts::scalar scalar_t, auto N>
requires(N >= 2) ALGEBRA_HOST_DEVICE
    constexpr auto phi(const Fastor::Tensor<scalar_t, N> &v) {
  return algebra::math::atan2(v[1], v[0]);
}

/// This method retrieves theta from a vector, vector base with rows >= 3
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
requires(N >= 3) ALGEBRA_HOST constexpr scalar_t
    theta(const Fastor::Tensor<scalar_t, N> &v) noexcept {

  return algebra::math::atan2(Fastor::norm(v(Fastor::fseq<0, 2>())), v[2]);
}

/// This method retrieves the perpendicular magnitude of a vector with rows >= 2
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
requires(N >= 2) ALGEBRA_HOST constexpr scalar_t
    perp(const Fastor::Tensor<scalar_t, N> &v) noexcept {

  return algebra::math::sqrt(
      Fastor::inner(v(Fastor::fseq<0, 2>()), v(Fastor::fseq<0, 2>())));
}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
ALGEBRA_HOST constexpr scalar_t norm(const Fastor::Tensor<scalar_t, N> &v) {

  return Fastor::norm(v);
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
requires(N >= 3) ALGEBRA_HOST constexpr scalar_t
    eta(const Fastor::Tensor<scalar_t, N> &v) noexcept {

  return algebra::math::atanh(v[2] / Fastor::norm(v));
}

/// Get a normalized version of the input vector
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
ALGEBRA_HOST constexpr Fastor::Tensor<scalar_t, N> normalize(
    const Fastor::Tensor<scalar_t, N> &v) {

  return (static_cast<scalar_t>(1.0) / Fastor::norm(v)) * v;
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
ALGEBRA_HOST_DEVICE constexpr scalar_t dot(
    const Fastor::Tensor<scalar_t, N> &a,
    const Fastor::Tensor<scalar_t, N> &b) {
  return Fastor::inner(a, b);
}

/// Dot product between Tensor<scalar_t, N> and Tensor<scalar_t, N, 1>
///
/// @param a the first input vector
/// @param b the second input Tensor<scalar_t, N, 1>
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
ALGEBRA_HOST constexpr scalar_t dot(const Fastor::Tensor<scalar_t, N> &a,
                                    const Fastor::Tensor<scalar_t, N, 1> &b) {

  // We need to specify the type of the Tensor slice because Fastor by default
  // is lazy, so it returns an intermediate type which does not play well with
  // the Fastor::inner function.
  return Fastor::inner(a,
                       Fastor::Tensor<scalar_t, N>(b(Fastor::fseq<0, N>(), 0)));
}

/// Dot product between Tensor<scalar_t, N> and Tensor<scalar_t, N, 1>
///
/// @param a the second input Tensor<scalar_t, N, 1>
/// @param b the first input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
ALGEBRA_HOST constexpr scalar_t dot(const Fastor::Tensor<scalar_t, N, 1> &a,
                                    const Fastor::Tensor<scalar_t, N> &b) {

  return Fastor::inner(Fastor::Tensor<scalar_t, N>(a(Fastor::fseq<0, N>(), 0)),
                       b);
}

/// Dot product between two Tensor<scalar_t, 3, 1>
///
/// @param a the second input Tensor<scalar_t, 3, 1>
/// @param b the first input Tensor<scalar_t, 3, 1>
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
ALGEBRA_HOST constexpr scalar_t dot(const Fastor::Tensor<scalar_t, N, 1> &a,
                                    const Fastor::Tensor<scalar_t, N, 1> &b) {

  return Fastor::inner(Fastor::Tensor<scalar_t, N>(a(Fastor::fseq<0, 3>(), 0)),
                       Fastor::Tensor<scalar_t, N>(b(Fastor::fseq<0, 3>(), 0)));
}

/// Cross product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE constexpr Fastor::Tensor<scalar_t, 3> cross(
    const Fastor::Tensor<scalar_t, 3> &a,
    const Fastor::Tensor<scalar_t, 3> &b) {
  return Fastor::cross(a, b);
}

/// Cross product between Tensor<scalar_t, 3> and Tensor<scalar_t, 3, 1>
///
/// @param a the first input vector
/// @param b the second input Tensor<scalar_t, 3, 1>
///
/// @return a vector representing the cross product
template <concepts::scalar scalar_t>
ALGEBRA_HOST constexpr Fastor::Tensor<scalar_t, 3> cross(
    const Fastor::Tensor<scalar_t, 3> &a,
    const Fastor::Tensor<scalar_t, 3, 1> &b) {

  // We need to specify the type of the Tensor slice because Fastor by default
  // is lazy, so it returns an intermediate type which does not play well with
  // the Fastor::cross function.
  return Fastor::cross(a,
                       Fastor::Tensor<scalar_t, 3>(b(Fastor::fseq<0, 3>(), 0)));
}

/// Cross product between Tensor<scalar_t, 3> and Tensor<scalar_t, 3, 1>
///
/// @param a the second input Tensor<scalar_t, 3, 1>
/// @param b the first input vector
///
/// @return a vector representing the cross product
template <concepts::scalar scalar_t>
ALGEBRA_HOST constexpr Fastor::Tensor<scalar_t, 3> cross(
    const Fastor::Tensor<scalar_t, 3, 1> &a,
    const Fastor::Tensor<scalar_t, 3> &b) {

  return Fastor::cross(Fastor::Tensor<scalar_t, 3>(a(Fastor::fseq<0, 3>(), 0)),
                       b);
}

/// Cross product between two Tensor<scalar_t, 3, 1>
///
/// @param a the second input Tensor<scalar_t, 3, 1>
/// @param b the first input Tensor<scalar_t, 3, 1>
///
/// @return a vector representing the cross product
template <concepts::scalar scalar_t>
ALGEBRA_HOST constexpr Fastor::Tensor<scalar_t, 3> cross(
    const Fastor::Tensor<scalar_t, 3, 1> &a,
    const Fastor::Tensor<scalar_t, 3, 1> &b) {

  return Fastor::cross(Fastor::Tensor<scalar_t, 3>(a(Fastor::fseq<0, 3>(), 0)),
                       Fastor::Tensor<scalar_t, 3>(b(Fastor::fseq<0, 3>(), 0)));
}

}  // namespace algebra::fastor::math
