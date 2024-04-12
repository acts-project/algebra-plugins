/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
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

namespace algebra::vc::math {

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

/// Get a normalized version of the input vector
///
/// @tparam vector_type generic input vector type
///
/// @param v the input vector
template <typename vector_type,
          std::enable_if_t<(Vc::is_simd_vector<vector_type>::value ||
                            algebra::detail::is_storage_vector_v<vector_type>),
                           bool> = true>
ALGEBRA_HOST_DEVICE inline auto normalize(const vector_type &v) {

  return v / algebra::math::sqrt(dot(v, v));
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
ALGEBRA_HOST_DEVICE inline auto cross(const vector_type1 &a,
                                      const vector_type2 &b)
    -> decltype(a * b - a * b) {

  return {a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0],
          a[0] * b[1] - b[0] * a[1], 0};
}

}  // namespace algebra::vc::math
