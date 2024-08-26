/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/matrix44.hpp"
#include "algebra/storage/vector.hpp"

// System include(s).
#include <array>
#include <cstddef>

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace algebra::vc_aos {

/// size type for Vc storage model
using size_type = std::size_t;
/// Array type used to store Vc::Vectors
template <typename T, size_type N>
using storage_type = Vc::SimdArray<T, N>;
/// value type in a linear algebra vector: SoA layout
template <typename T>
using value_type = T;
/// Vector type used in the Vc SoA storage model
template <typename T, std::size_t N>
using vector_type = storage::vector<N, T, storage_type>;
/// Matrix type used in the Vc SoA storage model
template <typename T, size_type ROWS, size_type COLS>
using matrix_type = storage::vector<ROWS * COLS, T, storage_type>;

/// 2-element "vector" type, using @c Vc::SimdArray
template <typename T>
using vector2 = vector_type<T, 2>;
/// Point in 2D space, using @c Vc::SimdArray
template <typename T>
using point2 = vector2<T>;
/// 3-element "vector" type, using @c Vc::SimdArray
template <typename T>
using vector3 = vector_type<T, 4>;
/// Point in 3D space, using @c Vc::SimdArray
template <typename T>
using point3 = vector3<T>;
/// 6-element "vector" type, using @c Vc::SimdArray
template <typename T>
using vector6 = vector_type<T, 6>;
/// 8-element "vector" type, using @c Vc::SimdArray
template <typename T>
using vector8 = vector_type<T, 8>;

}  // namespace algebra::vc_aos
