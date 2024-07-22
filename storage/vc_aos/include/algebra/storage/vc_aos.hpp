/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/impl/vc_aos_getter.hpp"
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

namespace vc_aos {

/// size type for Vc storage model
using size_type = std::size_t;
/// Array type used to store Vc::Vectors
template <typename T, size_type N>
using storage_type = Vc::SimdArray<T, N>;
/// value type in a linear algebra vector: AoS layout
template <typename T>
using value_type = T;
/// Vector type used in the Vc AoS storage model
template <typename T, std::size_t N>
using vector_type = algebra::storage::vector<N, T, storage_type>;
/// Matrix type used in the Vc AoS storage model
template <typename T, size_type ROWS, size_type COLS>
using matrix_type = algebra::storage::matrix<storage_type, T, ROWS, COLS>;

/// 2-element "vector" type, using @c Vc::SimdArray
template <typename T>
using vector2 = vector_type<T, 2>;
/// Point in 2D space, using @c Vc::SimdArray
template <typename T>
using point2 = vector2<T>;
/// 3-element "vector" type, using @c Vc::SimdArray
template <typename T>
using vector3 = vector_type<T, 3>;
/// Point in 3D space, using @c Vc::SimdArray
template <typename T>
using point3 = vector3<T>;
/// 6-element "vector" type, using @c Vc::SimdArray
template <typename T>
using vector6 = vector_type<T, 6>;
/// 8-element "vector" type, using @c Vc::SimdArray
template <typename T>
using vector8 = vector_type<T, 8>;

}  // namespace vc_aos

namespace trait {

/// Type trait specializations
/// @{

/// Index
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct index<algebra::storage::matrix<vc_aos::storage_type, T, ROWS, COLS>> {
  using type = algebra::vc_aos::size_type;
};

template <typename T, std::size_t N>
struct index<algebra::storage::vector<N, T, vc_aos::storage_type>> {
  using type = algebra::vc_aos::size_type;
};
/// @}

/// Dimension
/// @{
template <typename T, std::size_t N>
struct dimensions<algebra::storage::vector<N, T, vc_aos::storage_type>> {

  using size_type =
      index_t<algebra::storage::vector<N, T, vc_aos::storage_type>>;

  static constexpr size_type dim{1};
  static constexpr size_type rows{N};
  static constexpr size_type columns{1};
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct dimensions<
    algebra::storage::matrix<vc_aos::storage_type, T, ROWS, COLS>> {

  using size_type =
      index_t<algebra::storage::matrix<vc_aos::storage_type, T, ROWS, COLS>>;

  static constexpr size_type dim{2};
  static constexpr size_type rows{ROWS};
  static constexpr size_type columns{COLS};
};
/// @}

/// Value
/// @{
template <typename T, std::size_t N>
struct value<algebra::storage::vector<N, T, vc_aos::storage_type>> {
  using type = T;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct value<algebra::storage::matrix<vc_aos::storage_type, T, ROWS, COLS>> {
  using type = T;
};
/// @}

/// Matrix
/// @{
template <typename T, std::size_t N>
struct vector<algebra::storage::vector<N, T, vc_aos::storage_type>> {

  template <typename other_T, std::size_t other_N>
  using other_type =
      algebra::storage::vector<other_N, other_T, vc_aos::storage_type>;

  using type = other_type<T, N>;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct vector<algebra::storage::matrix<vc_aos::storage_type, T, ROWS, COLS>> {

  template <typename other_T, std::size_t other_N>
  using other_type =
      algebra::storage::vector<other_N, other_T, vc_aos::storage_type>;

  using type = other_type<T, ROWS>;
};
/// @}

/// Matrix
/// @{
template <typename T, std::size_t ROWS, std::size_t COLS>
struct matrix<algebra::storage::matrix<vc_aos::storage_type, T, ROWS, COLS>> {
  template <typename other_T, std::size_t other_ROWS, std::size_t other_COLS>
  using other_type = algebra::storage::matrix<vc_aos::storage_type, other_T,
                                              other_ROWS, other_COLS>;

  using type = algebra::storage::matrix<vc_aos::storage_type, T, ROWS, COLS>;
};

template <typename T, int N>
struct matrix<algebra::storage::vector<N, T, vc_aos::storage_type>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type = algebra::storage::matrix<vc_aos::storage_type, other_T,
                                              other_ROWS, other_COLS>;

  using type = other_type<T, N, 1>;
};
/// @}

/// Elemet/Block Getter
/// @{
template <typename T, std::size_t N>
struct element_getter<algebra::storage::vector<N, T, vc_aos::storage_type>> {
  using type = algebra::storage::element_getter;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct element_getter<
    algebra::storage::matrix<vc_aos::storage_type, T, ROWS, COLS>> {
  using type = algebra::storage::element_getter;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
struct block_getter<
    algebra::storage::matrix<vc_aos::storage_type, T, ROWS, COLS>> {
  using type = algebra::storage::block_getter;
};
/// @}
}  // namespace trait

}  // namespace algebra
