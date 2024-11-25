/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/qualifiers.hpp"

// SYCL include(s).
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

// System include(s).
#include <algorithm>
#include <cmath>

namespace algebra::math {

/// Namespace to pick up math functions from
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
namespace math_ns = ::sycl;
#else
namespace math_ns = std;
#endif  // SYCL

/// Absolute value of arg
template <concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline auto fabs(scalar_t arg) {
  return math_ns::fabs(arg);
}

/// Fused multiply add
template <concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline auto fma(scalar_t x, scalar_t y, scalar_t z) {
  return math_ns::fma(x, y, z);
}

/// Arc tangent of y/x
template <concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t atan2(scalar_t y, scalar_t x) {
  return math_ns::atan2(y, x);
}

/// Square root of arg
template <concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t sqrt(scalar_t arg) {
  return math_ns::sqrt(arg);
}

/// Inverse hyperbolic tangent of arg
template <concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t atanh(scalar_t arg) {
  return math_ns::atanh(arg);
}

/// Minimum of two values
template <concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t min(scalar_t a, scalar_t b) {
  return math_ns::min(a, b);
}

/// Maximum of two values
template <concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline scalar_t max(scalar_t a, scalar_t b) {
  return math_ns::max(a, b);
}

}  // namespace algebra::math
