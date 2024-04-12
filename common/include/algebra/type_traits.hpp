/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <algorithm>
#include <cmath>
#include <type_traits>

namespace algebra::trait {

/// Matrix traits
/// @{

/// Type of the matrix indices
/// @{
template <class M>
struct index {};

template <class M>
using index_t = typename index<M>::type;
/// @}

/// Value type that is used with the matrix (e.g. float or double precision)
/// @{
template <class M>
struct value {};

template <class M>
using value_t = typename value<M>::type;
/// @}

/// Scalar type that is used with the matrix (can be multiple values in SoA)
/// @{
template <class M>
struct scalar {
  using type = value_t<M>;
};

template <class M>
using scalar_t = typename scalar<M>::type;
/// @}

/// Vector type that is compatible with the matrix
/// @{
template <class M>
struct vector {};

template <class M>
using vector_t = typename vector<M>::type;

template <class M, index_t<M> N, typename T>
using get_vector_t = typename vector<M>::template other_type<T, N>;

template <class M>
struct matrix {};

template <class M, index_t<M> ROWS, index_t<M> COLS, typename T>
using other_matrix_t = typename matrix<M>::template other_type<T, ROWS, COLS>;
/// @}

/// Matrix dimensions
/// @{
template <typename M>
struct dimensions {};

/// Specilization for scalar types
template <typename M>
  requires std::is_fundamental_v<M>
struct dimensions<M> {

  using size_type = std::size_t;

  static constexpr size_type dim{0};
  static constexpr size_type rows{1};
  static constexpr size_type columns{1};
};

template <class M>
inline constexpr index_t<M> dim{dimensions<M>::sim};

template <class M>
inline constexpr index_t<M> rows{dimensions<M>::rows};

template <class M>
inline constexpr index_t<M> columns{dimensions<M>::columns};

template <class M>
inline constexpr index_t<M> rank{std::min(rows<M>, columns<M>)};

template <class M>
inline constexpr index_t<M> size{rows<M> * columns<M>};

template <class M>
inline constexpr bool is_matrix{dimensions<M>::dim == 2};

template <class M>
inline constexpr bool is_vector{dimensions<M>::dim == 1};

template <class M>
inline constexpr bool is_square{(rows<M> == columns<M>)};
/// @}

/// Getter types
/// @{
template <class M>
struct element_getter {};

template <class M>
using element_getter_t = typename element_getter<M>::type;

template <class M>
struct block_getter {};

template <class M>
using block_getter_t = typename block_getter<M>::type;
/// @}

/// @}

}  // namespace algebra::trait
