/** Algebra plugins, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// @TODO: Remove this when Vc fixes their false positives.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
#endif

// Project include(s)
#include "algebra/concepts.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/array_operators.hpp"

// System include(s).
#include <array>
#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <utility>

namespace algebra {

namespace storage {

namespace detail {
/// Make sure the vector/matrix dimension aligns with simd sizes
/// @see
/// https://gitlab.in2p3.fr/CodeursIntensifs/Fast5x5/-/blob/master/fast5x5.hpp?ref_type=heads
constexpr std::size_t nearest_power_of_two(std::size_t min_value,
                                           std::size_t current_value) {
  // Computes the nearest power of two relative to `min_value` starting from the
  // power of two `current_value`
  return min_value <= current_value
             ? current_value
             : nearest_power_of_two(min_value, current_value * 2u);
}

}  // namespace detail

/// Vector wrapper for AoS vs interleaved SoA data. @c scalar_t can e.g. be a
/// SIMD vector.
template <std::size_t N, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
class alignas(
    alignof(array_t<scalar_t, detail::nearest_power_of_two(N, 2u)>)) vector {

 public:
  /// @returns the required size of the underlying array storage
  static constexpr std::size_t simd_size() {
    return concepts::value<scalar_t> ? detail::nearest_power_of_two(N, 2u) : N;
  }

  // Value type is a simd vector in SoA and a scalar in AoS
  using scalar_type = scalar_t;
  /// Underlying data array type
  using array_type = array_t<scalar_t, simd_size()>;

  /// Default contructor sets all entries to zero.
  constexpr vector() {
    if constexpr (!concepts::simd_scalar<scalar_type>) {
      zero_fill(std::make_index_sequence<simd_size()>{});
    }
  }

  /// Construct from other vector type
  template <typename other_vector_t, typename S = scalar_type,
            std::enable_if_t<std::is_scalar_v<S>, bool> = true>
  constexpr vector(const other_vector_t &v) {
    for (std::size_t i = 0; i < N; ++i) {
      m_data[i] = v[i];
    }
  }

  /// Construct from element values @param vals .
  ///
  /// @{

  /// Construct vector in SoA layout from simd scalars
  template <typename... Scalars,
            std::enable_if_t<sizeof...(Scalars) == simd_size(), bool> = true>
    requires((concepts::simd_scalar<Scalars> ||
              std::convertible_to<Scalars, scalar_type>) &&
             ...)
  constexpr vector(Scalars &&...scals)
      : m_data{std::forward<Scalars>(scals)...} {}

  /// In order to avoid uninitialized values, which deteriorate the performance
  /// in explicitely vectorized code, the underlying data array is filled with
  /// zeroes if too few arguments are given.
  template <typename... Values,
            std::enable_if_t<(sizeof...(Values) < simd_size()), bool> = true>
    requires(std::convertible_to<Values, scalar_type> && ...)
  constexpr vector(Values &&...vals) {

    // Fill up last entries, if necessary (explicitly for now)
    if constexpr ((m_data.size() - sizeof...(Values)) == 1) {
      // 3D vectors
      m_data = {scalar_type{std::forward<Values>(vals)}..., scalar_type{0}};
    } else if constexpr ((m_data.size() - sizeof...(Values)) == 2) {
      // 6D vectors
      m_data = {scalar_type{std::forward<Values>(vals)}..., scalar_type{0}, scalar_type{0}};
    } else {
      // @TODO: Does this actually work, yet?
      zero_fill(std::make_index_sequence<simd_size() - sizeof...(Values)>{});
    }
  }

  /// @}

  /// Construct from existing array storage @param vals .
  constexpr vector(array_type &&vals) : m_data{std::move(vals)} {}
  constexpr vector(const array_type &vals) : m_data{vals} {}

  /// Conversion operator from wrapper to underlying data array.
  /// @{
  constexpr operator array_type &() { return m_data; }
  constexpr operator const array_type &() const { return m_data; }
  /// @}

  constexpr const auto &get() const { return m_data; }

  /// Subscript operator[]
  /// @{
  constexpr decltype(auto) operator[](std::size_t i) { return m_data[i]; }
  constexpr decltype(auto) operator[](std::size_t i) const { return m_data[i]; }
  /// @}

  /// @returns the size of the underlying data storage
  static constexpr std::size_t storage_size() { return simd_size(); }

  /// Operator*=
  ///
  /// @return Vector expression/return type according to the operation.
  constexpr decltype(auto) operator*=(scalar_type factor) noexcept {
    return m_data *= factor;
  }

  /// Equality operators
  /// @{
  template <std::size_t M, concepts::scalar o_scalar_t,
            template <typename, std::size_t> class o_array_t,
            template <typename, std::size_t> class p_array_t>
  friend constexpr bool operator==(
      const vector<M, o_scalar_t, o_array_t> &,
      const vector<M, o_scalar_t, p_array_t> &) noexcept;

  template <std::size_t M, concepts::scalar o_scalar_t,
            template <typename, std::size_t> class o_array_t,
            template <typename, std::size_t> class p_array_t, bool>
  friend constexpr bool operator==(const vector<M, o_scalar_t, o_array_t> &,
                                   const p_array_t<o_scalar_t, M> &) noexcept;
  /// @}

  /// Inequality operator
  template <typename other_type>
  constexpr bool operator!=(const other_type &rhs) const noexcept {
    return ((*this == rhs) == false);
  }

  /// Elementwise comparison. Can result in a vector-of-masks for SoA vectors
  template <typename other_type>
  constexpr auto compare(const other_type &rhs) const noexcept {
    using result_t = decltype(m_data[0] == rhs[0]);

    std::array<result_t, N> comp;

    for (unsigned int i{0u}; i < N; ++i) {
      comp[i] = (m_data[i] == rhs[i]);
    }

    return comp;
  }

  /// Holds the data value for every vector element
  array_t<scalar_t, simd_size()> m_data;

 private:
  /// Sets the trailing uninitialized values to zero.
  template <std::size_t... Is>
  constexpr void zero_fill(std::index_sequence<Is...>) noexcept {
    ((m_data[simd_size() - sizeof...(Is) + Is] = scalar_t(0)), ...);
  }
};

/// Friend operators
/// @{

template <std::size_t N, std::size_t M, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t,
          template <typename, std::size_t> class o_array_t,
          std::enable_if_t<std::is_scalar_v<scalar_t>, bool> = true>
constexpr bool operator==(const vector<N, scalar_t, array_t> &lhs,
                          const o_array_t<scalar_t, M> &rhs) noexcept {

  const auto comp = lhs.compare(rhs);
  bool is_full = false;

  for (unsigned int i{0u}; i < N; ++i) {
    is_full |= comp[i];
  }

  return is_full;
}

template <std::size_t N, std::size_t M, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t,
          template <typename, std::size_t> class o_array_t,
          std::enable_if_t<!std::is_scalar_v<scalar_t>, bool> = true>
constexpr bool operator==(const vector<N, scalar_t, array_t> &lhs,
                          const o_array_t<scalar_t, M> &rhs) noexcept {

  const auto comp = lhs.compare(rhs);
  bool is_full = false;

  for (unsigned int i{0u}; i < N; ++i) {
    // Ducktyping the Vc::Vector::MaskType
    is_full |= comp[i].isFull();
  }

  return is_full;
}

template <std::size_t N, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t,
          template <typename, std::size_t> class o_array_t>
constexpr bool operator==(const vector<N, scalar_t, array_t> &lhs,
                          const vector<N, scalar_t, o_array_t> &rhs) noexcept {
  return (rhs.m_data == lhs);
}

/// @}

/// Macro declaring all instances of a specific arithmetic operator
#define DECLARE_vector_OPERATORS(OP)                                           \
  template <std::size_t N, concepts::scalar scalar_t, concepts::value value_t, \
            template <typename, std::size_t> class array_t>                    \
  inline constexpr decltype(auto) operator OP(                                 \
      const vector<N, scalar_t, array_t> &lhs, value_t rhs) noexcept {         \
    return lhs.m_data OP static_cast<scalar_t>(rhs);                           \
  }                                                                            \
  template <std::size_t N, concepts::scalar scalar_t, concepts::value value_t, \
            template <typename, std::size_t> class array_t>                    \
  inline decltype(auto) operator OP(                                           \
      value_t lhs, const vector<N, scalar_t, array_t> &rhs) noexcept {         \
    return static_cast<scalar_t>(lhs) OP rhs.m_data;                           \
  }                                                                            \
  template <std::size_t N, concepts::scalar scalar_t,                          \
            template <typename, std::size_t> class array_t>                    \
  inline constexpr decltype(auto) operator OP(                                 \
      const vector<N, scalar_t, array_t> &lhs,                                 \
      const vector<N, scalar_t, array_t> &rhs) noexcept {                      \
    return lhs.m_data OP rhs.m_data;                                           \
  }                                                                            \
  template <std::size_t N, concepts::scalar scalar_t,                          \
            template <typename, std::size_t> class array_t,                    \
            typename other_type>                                               \
    requires(!concepts::value<other_type>)                                     \
  inline constexpr decltype(auto) operator OP(                                 \
      const vector<N, scalar_t, array_t> &lhs,                                 \
      const other_type &rhs) noexcept {                                        \
    return lhs.m_data OP rhs;                                                  \
  }                                                                            \
  template <std::size_t N, concepts::scalar scalar_t,                          \
            template <typename, std::size_t> class array_t,                    \
            typename other_type>                                               \
    requires(!concepts::value<other_type>)                                     \
  inline constexpr decltype(auto) operator OP(                                 \
      const other_type &lhs,                                                   \
      const vector<N, scalar_t, array_t> &rhs) noexcept {                      \
    return lhs OP rhs.m_data;                                                  \
  }

// Implement all arithmetic operations on top of @c vector.
// clang-format off
DECLARE_vector_OPERATORS(+)
DECLARE_vector_OPERATORS(-)
DECLARE_vector_OPERATORS(*)
DECLARE_vector_OPERATORS(/)
// clang-format on

// Clean up.
#undef DECLARE_vector_OPERATORS

}  // namespace storage

namespace detail {

template <typename T>
struct is_storage_vector : public std::false_type {};

template <std::size_t N, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
struct is_storage_vector<storage::vector<N, scalar_t, array_t>>
    : public std::true_type {};

template <typename T>
inline constexpr bool is_storage_vector_v = is_storage_vector<T>::value;

}  // namespace detail

}  // namespace algebra
