/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace algebra::math {

/// Vc overloads of common math functions
/// @{
template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) abs(T &&vec) {
  return Vc::abs(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) fabs(T &&vec) {
  return Vc::abs(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) sqrt(T &&vec) {
  return Vc::sqrt(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) exp(T &&vec) {
  return Vc::exp(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) log(T &&vec) {
  return Vc::log(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) sin(T &&vec) {
  return Vc::sin(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) asin(T &&vec) {
  return Vc::asin(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) cos(T &&vec) {
  return Vc::cos(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) tan(T &&vec) {
  // It seems there is no dedicated @c Vc::tan function ?
  return Vc::sin(std::forward<T>(vec)) / Vc::cos(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) atan(T &&vec) {
  return Vc::atan(std::forward<T>(vec));
}

template <typename T, typename S>
  requires Vc::Traits::is_simd_vector<T>::value &&
           Vc::Traits::is_simd_vector<S>::value
inline decltype(auto) copysign(T &&mag, S &&sgn) {
  return Vc::copysign(std::forward<T>(mag), std::forward<S>(sgn));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) min(T &&vec) {
  return Vc::min(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) max(T &&vec) {
  return Vc::max(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) signbit(T &&vec) {
  return Vc::isnegative(std::forward<T>(vec));
}

template <typename T>
  requires Vc::Traits::is_simd_vector<T>::value
inline decltype(auto) fma(T &&x, T &&y, T &&z) {
  return Vc::fma(std::forward<T>(x), std::forward<T>(y), std::forward<T>(z));
}
/// @}

}  // namespace algebra::math
