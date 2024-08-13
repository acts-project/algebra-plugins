/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/type_traits.hpp"

// System include(s).
#include <concepts>

namespace algebra::concepts {

// Value concept: Single entry
template <typename T>
concept value = std::is_arithmetic_v<std::decay_t<T>>;

/// Scalar concept: Elements of vectors/matrices (can be simd vectors)
template <typename T>
concept scalar = !algebra::traits::is_matrix<T> &&
                 !algebra::traits::is_vector<T> && requires(T a, T b) {
                   { a + b } -> std::convertible_to<T>;
                   { a - b } -> std::convertible_to<T>;
                   { a* b } -> std::convertible_to<T>;
                   { a / b } -> std::convertible_to<T>;
                 };

/// Check if a scalar is simd
template <typename T>
concept simd_scalar = scalar<T> && !std::is_scalar_v<T>;

/// Index concept to access vector/matrix elements
template <typename T>
concept index = std::is_integral_v<T>;

/// Vector concepts
/// @{
template <typename V>
concept vector = algebra::traits::is_vector<V>;

template <typename V>
concept vector2D = vector<V> && (algebra::traits::size<V> == 2);

template <typename V>
concept vector3D = vector<V> && (algebra::traits::size<V> == 3);
/// @}

/// Point concepts
/// @{
template <typename V>
concept point = vector<V>;

template <typename V>
concept point2D = point<V> && (algebra::traits::size<V> == 2);

template <typename V>
concept point3D = point<V> && (algebra::traits::size<V> == 3);
/// @}

/// Matrix concepts
/// @{
template <typename M>
concept matrix = algebra::traits::is_matrix<M>;

template <typename M>
concept square_matrix = matrix<M> && algebra::traits::is_square<M>;
/// @}

/// Transform concept
template <typename T>
concept transform3D = requires(T trf) {
  // Local type definitions
  requires scalar<typename T::scalar_type>;
  requires vector3D<typename T::vector3>;
  requires point2D<typename T::point2>;
  requires point3D<typename T::point3>;

  // Methods
  trf.rotation();
  trf.translation();
  trf.point_to_global(typename T::vector3());
  trf.point_to_local(typename T::vector3());
  trf.vector_to_global(typename T::vector3());
  trf.vector_to_local(typename T::vector3());
};

}  // namespace algebra::concepts
