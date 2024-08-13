/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/vector.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s).
#include <type_traits>
#include <utility>

namespace algebra::vc_aos::math {

/// This method retrieves phi from a vector @param v
template <typename vector_t,
          std::enable_if_t<(Vc::is_simd_vector<vector_t>::value ||
                            algebra::detail::is_storage_vector_v<vector_t>),
                           bool> = true>
ALGEBRA_HOST_DEVICE inline auto phi(const vector_t &v) {
  return algebra::math::atan2(v[1], v[0]);
}

/// This method retrieves the perpendicular magnitude of a vector @param v
template <typename vector_t,
          std::enable_if_t<(Vc::is_simd_vector<vector_t>::value ||
                            algebra::detail::is_storage_vector_v<vector_t>),
                           bool> = true>
ALGEBRA_HOST_DEVICE inline auto perp(const vector_t &v) {
  return algebra::math::sqrt(algebra::math::fma(v[0], v[0], v[1] * v[1]));
}

/// This method retrieves theta from a vector @param v
template <typename vector_t,
          std::enable_if_t<(Vc::is_simd_vector<vector_t>::value ||
                            algebra::detail::is_storage_vector_v<vector_t>),
                           bool> = true>
ALGEBRA_HOST_DEVICE inline auto theta(const vector_t &v) {
  return algebra::math::atan2(perp(v), v[2]);
}

/// Dot product between two input vectors
///
/// @tparam vector_type generic input vector type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <typename vector_type1, typename vector_type2,
          std::enable_if_t<
              ((Vc::is_simd_vector<vector_type1>::value ||
                algebra::detail::is_storage_vector_v<
                    vector_type1>)&&(Vc::is_simd_vector<vector_type2>::value ||
                                     algebra::detail::is_storage_vector_v<
                                         vector_type2>)),
              bool> = true>
ALGEBRA_HOST_DEVICE inline auto dot(const vector_type1 &a,
                                    const vector_type2 &b) {

  return (a * b).sum();
}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @param v the input vector
template <typename vector_t,
          std::enable_if_t<(Vc::is_simd_vector<vector_t>::value ||
                            algebra::detail::is_storage_vector_v<vector_t>),
                           bool> = true>
ALGEBRA_HOST_DEVICE inline auto norm(const vector_t &v) {

  return algebra::math::sqrt(dot(v, v));
}

/// Get a normalized version of the input vector
///
/// @tparam vector_type generic input vector type
///
/// @param v the input vector
template <typename vector_type,
          std::enable_if_t<(Vc::is_simd_vector<vector_type>::value ||
                            algebra::detail::is_storage_vector_v<vector_type>),
                           bool> = true>
ALGEBRA_HOST_DEVICE inline vector_type normalize(const vector_type &v) {

  return v / norm(v);
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <typename vector_t,
          std::enable_if_t<(Vc::is_simd_vector<vector_t>::value ||
                            algebra::detail::is_storage_vector_v<vector_t>),
                           bool> = true>
ALGEBRA_HOST_DEVICE inline auto eta(const vector_t &v) noexcept {

  return algebra::math::atanh(v[2] / norm(v));
}

/// Cross product between two input vectors - 3 Dim
///
/// @tparam vector_type generic input vector type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector representing the cross product
template <typename vector_type1, typename vector_type2,
          std::enable_if_t<
              ((Vc::is_simd_vector<vector_type1>::value ||
                algebra::detail::is_storage_vector_v<
                    vector_type1>)&&(Vc::is_simd_vector<vector_type2>::value ||
                                     algebra::detail::is_storage_vector_v<
                                         vector_type2>)),
              bool> = true>
ALGEBRA_HOST_DEVICE inline vector_type1 cross(const vector_type1 &a,
                                              const vector_type2 &b) {

  return {algebra::math::fma(a[1], b[2], -b[1] * a[2]),
          algebra::math::fma(a[2], b[0], -b[2] * a[0]),
          algebra::math::fma(a[0], b[1], -b[0] * a[1]), 0.f};
}

/// Elementwise sum
///
/// @tparam vector_type generic input vector type
///
/// @param v the vector whose elements should be summed
///
/// @return the sum of the elements
template <typename vector_type,
          std::enable_if_t<(Vc::is_simd_vector<vector_type>::value ||
                            algebra::detail::is_storage_vector_v<vector_type>),
                           bool> = true>
ALGEBRA_HOST_DEVICE inline auto sum(const vector_type &v) {
  return v.get().sum();
}

}  // namespace algebra::vc_aos::math
