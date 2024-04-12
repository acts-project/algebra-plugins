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
#include "algebra/qualifiers.hpp"
#include "algebra/storage/array_operators.hpp"

// System include(s).
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
ALGEBRA_HOST_DEVICE
consteval std::size_t nearest_power_of_two(std::size_t min_value,
                                           std::size_t current_value) {
  // Computes the nearest power of two relative to `min_value` starting from the
  // power of two `current_value`
  return min_value <= current_value
             ? current_value
             : nearest_power_of_two(min_value, current_value * 2u);
}

}  // namespace detail

/// Vector wrapper for AoS vs interleaved SoA data. @c value_t can e.g. be a
/// SIMD vector.
template <std::size_t N, typename value_t,
          template <typename, std::size_t> class array_t>
class alignas(
    alignof(array_t<value_t, detail::nearest_power_of_two(N, 2u)>)) vector {

 public:
  /// @returns the required size of the underlying array storage
  ALGEBRA_HOST_DEVICE
  static consteval std::size_t simd_size() {
    return std::is_scalar_v<value_t> ? detail::nearest_power_of_two(N, 2u) : N;
  }

  // Value type is a simd vector in SoA and a scalar in AoS
  using value_type = value_t;
  /// Underlying data array type
  using array_type = array_t<value_t, simd_size()>;

  /// Default contructor sets all entries to zero.
  ALGEBRA_HOST_DEVICE
  constexpr vector() {
    if constexpr (std::is_scalar_v<value_type>) {
      zero_fill(std::make_index_sequence<simd_size()>{});
    }
  }

  /// Construct from element values @param vals .
  ///
  /// In order to avoid uninitialized values, which deteriorate the performance
  /// in explicitely vectorized code, the underlying data array is filled with
  /// zeroes if too few arguments are given.
  template <typename... Values>
    requires(std::conjunction_v<std::is_convertible<Values, value_type>...> &&
             sizeof...(Values) <= N && !(N == 3 && simd_size() == 4) &&
             !(N == 6 && simd_size() == 8))
  ALGEBRA_HOST_DEVICE constexpr vector(Values &&...vals)
      : m_data{std::forward<Values>(vals)...} {
    // Fill the uninitialized part of the vector register with zero
    if constexpr ((sizeof...(Values) < simd_size()) &&
                  std::is_scalar_v<value_type> &&
                  (!std::conjunction_v<std::is_same<array_type, Values>...>)) {
      zero_fill(std::make_index_sequence<simd_size() - sizeof...(Values)>{});
    }
  }

  template <typename... Values>
    requires(std::conjunction_v<std::is_convertible<Values, value_type>...> &&
             N == 3 && simd_size() == 4)
  ALGEBRA_HOST_DEVICE constexpr vector(Values &&...vals)
      : m_data{std::forward<Values>(vals)..., 0.f} {}

  template <typename... Values>
    requires(std::conjunction_v<std::is_convertible<Values, value_type>...> &&
             N == 6 && simd_size() == 8)
  ALGEBRA_HOST_DEVICE constexpr vector(Values &&...vals)
      : m_data{std::forward<Values>(vals)..., 0.f, 0.f} {}

  /// Construct from existing array storage @param vals
  ALGEBRA_HOST_DEVICE
  constexpr vector(const array_type &vals) : m_data{vals} {}

  /// Assignment operator from wrapped data.
  ///
  /// @param lhs wrap a copy of this data.
  ALGEBRA_HOST_DEVICE
  constexpr const vector &operator=(const array_type &lhs) {
    m_data = lhs;
    return *this;
  }

  /// Assignment operator from @c std::initializer_list .
  ///
  /// @param list wrap an array of this data.
  ALGEBRA_HOST_DEVICE
  constexpr vector &operator=(std::initializer_list<value_type> &list) {
    m_data = array_type(list);
    return *this;
  }

  /// Conversion operator from wrapper to underlying data array.
  /// @{
  ALGEBRA_HOST_DEVICE
  constexpr operator array_type &() { return m_data; }
  ALGEBRA_HOST_DEVICE
  constexpr operator const array_type &() const { return m_data; }
  /// @}

  ALGEBRA_HOST_DEVICE
  constexpr const auto &get() const { return m_data; }

  /// Subscript operator[]
  /// @{
  ALGEBRA_HOST_DEVICE
  constexpr decltype(auto) operator[](std::size_t i) { return m_data[i]; }
  ALGEBRA_HOST_DEVICE
  constexpr decltype(auto) operator[](std::size_t i) const { return m_data[i]; }
  /// @}

  /// Operator*=
  ///
  /// @return Vector expression/return type according to the operation.
  ALGEBRA_HOST_DEVICE
  constexpr decltype(auto) operator*=(value_type factor) noexcept {
    return m_data *= factor;
  }

  /// Equality operators
  /// @{
  template <std::size_t M, typename o_value_t,
            template <typename, std::size_t> class o_array_t,
            template <typename, std::size_t> class p_array_t>
  ALGEBRA_HOST_DEVICE friend constexpr bool operator==(
      const vector<M, o_value_t, o_array_t> &,
      const vector<M, o_value_t, p_array_t> &) noexcept;

  template <std::size_t M, typename o_value_t,
            template <typename, std::size_t> class o_array_t,
            template <typename, std::size_t> class p_array_t, bool>
  ALGEBRA_HOST_DEVICE friend constexpr bool operator==(
      const vector<M, o_value_t, o_array_t> &,
      const p_array_t<o_value_t, M> &) noexcept;
  /// @}

  /// Inequality operator
  template <typename other_type>
  ALGEBRA_HOST_DEVICE constexpr bool operator!=(
      const other_type &rhs) const noexcept {
    return ((*this == rhs) == false);
  }

  /// Elementwise comparison. Can result in a vector-of-masks for SoA vectors
  template <typename other_type>
  ALGEBRA_HOST_DEVICE constexpr auto compare(
      const other_type &rhs) const noexcept {
    using result_t = decltype(m_data[0] == rhs[0]);

    array_t<result_t, N> comp;

    for (unsigned int i{0u}; i < N; ++i) {
      comp[i] = (m_data[i] == rhs[i]);
    }

    return comp;
  }

  /// Holds the data value for every vector element
  array_t<value_t, simd_size()> m_data;

 private:
  /// Sets the trailing uninitialized values to zero.
  template <std::size_t... Is>
  ALGEBRA_HOST_DEVICE constexpr void zero_fill(
      std::index_sequence<Is...>) noexcept {
    ((m_data[simd_size() - sizeof...(Is) + Is] = value_t(0)), ...);
  }
};

/// Friend operators
/// @{

template <std::size_t N, typename value_t,
          template <typename, std::size_t> class array_t,
          template <typename, std::size_t> class o_array_t>
  requires(std::is_scalar_v<value_t>)
ALGEBRA_HOST_DEVICE constexpr bool operator==(
    const vector<N, value_t, array_t> &lhs,
    const o_array_t<value_t, N> &rhs) noexcept {

  const auto comp = lhs.compare(rhs);
  bool is_full = false;

  for (unsigned int i{0u}; i < N; ++i) {
    is_full |= comp[i];
  }

  return is_full;
}

template <std::size_t N, typename value_t,
          template <typename, std::size_t> class array_t,
          template <typename, std::size_t> class o_array_t>
  requires(!std::is_scalar_v<value_t>)
ALGEBRA_HOST_DEVICE constexpr bool operator==(
    const vector<N, value_t, array_t> &lhs,
    const o_array_t<value_t, N> &rhs) noexcept {

  const auto comp = lhs.compare(rhs);
  bool is_full = false;

  for (unsigned int i{0u}; i < N; ++i) {
    // Ducktyping the Vc::Vector::MaskType
    is_full |= comp[i].isFull();
  }

  return is_full;
}
template <std::size_t N, typename value_t,
          template <typename, std::size_t> class array_t,
          template <typename, std::size_t> class o_array_t>
ALGEBRA_HOST_DEVICE constexpr bool operator==(
    const vector<N, value_t, array_t> &lhs,
    const vector<N, value_t, o_array_t> &rhs) noexcept {
  return (lhs == rhs.m_data);
}

/// @}

/// Macro declaring all instances of a specific arithmetic operator
#define DECLARE_vector_OPERATORS(OP)                                          \
  template <std::size_t N, typename value_t, typename scalar_t,               \
            template <typename, std::size_t> class array_t>                   \
    requires(std::is_scalar_v<scalar_t>)                                      \
  ALGEBRA_HOST_DEVICE constexpr decltype(auto) operator OP(                   \
      const vector<N, value_t, array_t> &lhs, scalar_t rhs) noexcept {        \
    return lhs.m_data OP static_cast<value_t>(rhs);                           \
  }                                                                           \
  template <std::size_t N, typename value_t, typename scalar_t,               \
            template <typename, std::size_t> class array_t>                   \
    requires(std::is_scalar_v<scalar_t>)                                      \
  ALGEBRA_HOST_DEVICE inline decltype(auto) operator OP(                      \
      scalar_t lhs, const vector<N, value_t, array_t> &rhs) noexcept {        \
    return static_cast<value_t>(lhs) OP rhs.m_data;                           \
  }                                                                           \
  template <std::size_t N, typename value_t,                                  \
            template <typename, std::size_t> class array_t>                   \
  ALGEBRA_HOST_DEVICE constexpr decltype(auto) operator OP(                   \
      const vector<N, value_t, array_t> &lhs,                                 \
      const vector<N, value_t, array_t> &rhs) noexcept {                      \
    return lhs.m_data OP rhs.m_data;                                          \
  }                                                                           \
  template <std::size_t N, typename value_t,                                  \
            template <typename, std::size_t> class array_t,                   \
            typename other_type>                                              \
    requires(                                                                 \
        std::is_object<decltype(std::declval<typename vector<                 \
                                    N, value_t, array_t>::array_type>()       \
                                    OP std::declval<other_type>())>::value && \
        !std::is_scalar_v<other_type>)                                        \
  ALGEBRA_HOST_DEVICE constexpr decltype(auto) operator OP(                   \
      const vector<N, value_t, array_t> &lhs,                                 \
      const other_type &rhs) noexcept {                                       \
    return lhs.m_data OP rhs;                                                 \
  }                                                                           \
  template <std::size_t N, typename value_t,                                  \
            template <typename, std::size_t> class array_t,                   \
            typename other_type>                                              \
    requires(                                                                 \
        std::is_object<decltype(std::declval<typename vector<                 \
                                    N, value_t, array_t>::array_type>()       \
                                    OP std::declval<other_type>())>::value && \
        !std::is_scalar_v<other_type>)                                        \
  ALGEBRA_HOST_DEVICE constexpr decltype(auto) operator OP(                   \
      const other_type &lhs,                                                  \
      const vector<N, value_t, array_t> &rhs) noexcept {                      \
    return lhs OP rhs.m_data;                                                 \
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

template <std::size_t N, typename value_t,
          template <typename, std::size_t> class array_t>
struct is_storage_vector<storage::vector<N, value_t, array_t>>
    : public std::true_type {};

template <typename T>
inline constexpr bool is_storage_vector_v = is_storage_vector<T>::value;

}  // namespace detail

}  // namespace algebra
