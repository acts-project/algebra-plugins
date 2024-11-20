/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/algorithms/matrix/determinant/hard_coded.hpp"
#include "algebra/math/algorithms/matrix/determinant/partial_pivot_lud.hpp"
#include "algebra/math/algorithms/matrix/inverse/hard_coded.hpp"
#include "algebra/math/algorithms/matrix/inverse/partial_pivot_lud.hpp"

namespace algebra::generic {

/// Get the type of determinant algorithm acording to matrix dimension
/// @{
template <std::size_t N, typename... Args>
struct determinant_selector {
  using type = matrix::determinant::partial_pivot_lud<Args...>;
};

template <typename... Args>
struct determinant_selector<2, Args...> {
  using type = matrix::determinant::hard_coded<Args...>;
};

template <typename... Args>
struct determinant_selector<4, Args...> {
  using type = matrix::determinant::hard_coded<Args...>;
};

/// @tparam M matrix type
template <class M>
using determinant_t =
    typename determinant_selector<algebra::trait::rank<M>, M,
                                  algebra::trait::element_getter_t<M>>::type;
/// @}

/// Get the type of inversion algorithm acording to matrix dimension
/// @{
template <std::size_t N, typename... Args>
struct inversion_selector {
  using type = matrix::inverse::partial_pivot_lud<Args...>;
};

template <typename... Args>
struct inversion_selector<2, Args...> {
  using type = matrix::inverse::hard_coded<Args...>;
};

template <typename... Args>
struct inversion_selector<4, Args...> {
  using type = matrix::inverse::hard_coded<Args...>;
};

/// @tparam M matrix type
template <class M>
using inversion_t =
    typename inversion_selector<algebra::trait::rank<M>, M,
                                algebra::trait::element_getter_t<M>>::type;
/// @}generic

}  // namespace algebra::generic
