/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/boolean.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace algebra::boolean {

/// Boolean utilities on single values
/// @{
using algebra::boolean::all_of;
using algebra::boolean::any_of;
using algebra::boolean::none_of;
/// @}

/// Vc overloads of boolean utilities
/// @{
template <typename T>
requires Vc::Traits::is_simd_mask<T>::value inline bool any_of(T &&mask) {
  return Vc::any_of(std::forward<T>(mask));
}

template <typename T>
requires Vc::Traits::is_simd_mask<T>::value inline bool all_of(T &&mask) {
  return Vc::all_of(std::forward<T>(mask));
}

template <typename T>
requires Vc::Traits::is_simd_mask<T>::value inline bool none_of(T &&mask) {
  return Vc::none_of(std::forward<T>(mask));
}
/// @}

}  // namespace algebra::boolean
