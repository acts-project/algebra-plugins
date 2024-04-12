/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/matrix.hpp"
#include "algebra/storage/vector.hpp"
#include "algebra/type_traits.hpp"

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

namespace algebra {

namespace vc_soa {

/// size type for Vc storage model
using size_type = std::size_t;
/// Array type used to store Vc::Vectors
template <typename T, size_type N>
using storage_type = std::array<T, N>;
/// value type in a linear algebra vector: SoA layout
template <typename T>
using value_type = Vc::Vector<T>;
/// Vector type used in the Vc SoA storage model
template <typename T, std::size_t N>
using vector_type = storage::vector<N, value_type<T>, storage_type>;
/// Matrix type used in the Vc SoA storage model
template <typename T, size_type ROWS, size_type COLS>
using matrix_type = storage::matrix<storage_type, value_type<T>, ROWS, COLS>;

/// 2-element "vector" type, using @c Vc::Vector in every element
template <typename T>
using vector2 = vector_type<T, 2>;
/// Point in 2D space, using @c Vc::Vector in every element
template <typename T>
using point2 = vector2<T>;
/// 3-element "vector" type, using @c Vc::Vector in every element
template <typename T>
using vector3 = vector_type<T, 3>;
/// Point in 3D space, using @c Vc::Vector in every element
template <typename T>
using point3 = vector3<T>;
/// 6-element "vector" type, using @c Vc::Vector in every element
template <typename T>
using vector6 = vector_type<T, 6>;
/// 8-element "vector" type, using @c Vc::Vector in every element
template <typename T>
using vector8 = vector_type<T, 8>;

}  // namespace vc_soa

namespace trait {

/// Type trait specializations
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<storage::matrix<std::array, Vc::Vector<T>, ROWS, COLS>> {
  using type = algebra::vc_soa::size_type;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<storage::matrix<std::array, Vc::Vector<T>, ROWS, COLS>> {

  using size_type =
      index_t<storage::matrix<std::array, Vc::Vector<T>, ROWS, COLS>>;

  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<storage::matrix<std::array, Vc::Vector<T>, ROWS, COLS>> {
  using type = T;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct scalar<storage::matrix<std::array, Vc::Vector<T>, ROWS, COLS>> {
  using type = Vc::Vector<T>;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<storage::matrix<std::array, Vc::Vector<T>, ROWS, COLS>> {
  using type = storage::vector<ROWS, Vc::Vector<T>, std::array>;
};
/// @}

}  // namespace trait

}  // namespace algebra
