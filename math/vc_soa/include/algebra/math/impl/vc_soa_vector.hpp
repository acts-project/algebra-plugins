/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/vector.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace algebra::vc_soa::math {

/// This method retrieves phi from a vector, vector base with rows >= 2
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
requires(N >= 2) ALGEBRA_HOST_DEVICE
    inline auto phi(const storage::vector<N, Vc::Vector<value_t>, array_t> &v) {

  return Vc::atan2(v[1], v[0]);
}

/// This method retrieves the perpenticular magnitude of a vector with rows >= 2
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
requires(N >= 2) ALGEBRA_HOST_DEVICE inline auto perp(
    const storage::vector<N, Vc::Vector<value_t>, array_t> &v) {

  return Vc::sqrt(Vc::fma(v[0], v[0], v[1] * v[1]));
}

/// This method retrieves theta from a vector, vector base with rows >= 3
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
requires(N >= 3) ALGEBRA_HOST_DEVICE inline auto theta(
    const storage::vector<N, Vc::Vector<value_t>, array_t> &v) {

  return Vc::atan2(perp(v), v[2]);
}

/// Cross product between two input vectors - 3 Dim
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam array_t array type that holds the vector elements
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
requires(N == 3) ALGEBRA_HOST_DEVICE
    inline storage::vector<N, Vc::Vector<value_t>, array_t> cross(
        const storage::vector<N, Vc::Vector<value_t>, array_t> &a,
        const storage::vector<N, Vc::Vector<value_t>, array_t> &b) {

  return {Vc::fma(a[1], b[2], -b[1] * a[2]), Vc::fma(a[2], b[0], -b[2] * a[0]),
          Vc::fma(a[0], b[1], -b[0] * a[1])};
}

/// Dot product between two input vectors
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam array_t array type that holds the vector elements
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline Vc::Vector<value_t> dot(
    const storage::vector<N, Vc::Vector<value_t>, array_t> &a,
    const storage::vector<N, Vc::Vector<value_t>, array_t> &b) {

  auto ret = a[0] * b[0];

  for (unsigned int i{1u}; i < N; i++) {
    ret = Vc::fma(a[i], b[i], ret);
  }

  return ret;
}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline auto norm(
    const storage::vector<N, Vc::Vector<value_t>, array_t> &v) {

  return Vc::sqrt(dot(v, v));
}

/// Get a normalized version of the input vector
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline storage::vector<N, Vc::Vector<value_t>, array_t>
normalize(const storage::vector<N, Vc::Vector<value_t>, array_t> &v) {

  return (Vc::Vector<value_t>::One() / norm(v)) * v;

  // Less accurate, but faster
  // return Vc::reciprocal(norm(v)) * v;

  // Even less accurate, but even faster
  // return Vc::rsqrt(dot(v, v)) * v;
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
requires(N >= 3) ALGEBRA_HOST_DEVICE
    inline auto eta(const storage::vector<N, Vc::Vector<value_t>, array_t> &v) {

  // atanh does not exist in Vc
  auto atanh_func = [](value_t e) { return std::atanh(e); };

  return (v[2] / norm(v)).apply(atanh_func);

  // Faster, but less accurate
  // return (Vc::reciprocal(norm(v)) * v[2]).apply(atanh_func);

  // Even faster, but even less accurate
  // return (Vc::rsqrt(dot(v, v)) * v[2]).apply(atanh_func);
}

/// Elementwise sum
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam array_t array type that holds the vector elements
///
/// @param v the vector whose elements should be summed
///
/// @return the sum of the elements
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
ALGEBRA_HOST_DEVICE inline Vc::Vector<value_t> sum(
    const storage::vector<N, Vc::Vector<value_t>, array_t> &v) {

  Vc::Vector<value_t> res{v[0]};

  for (std::size_t i = 1u; i < Vc::Vector<value_t>::size(); ++i) {
    res = res + v[i];
  }

  return res;
}

}  // namespace algebra::vc_soa::math
