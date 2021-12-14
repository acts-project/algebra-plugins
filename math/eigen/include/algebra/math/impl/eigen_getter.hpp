/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// Eigen include(s).
#include <Eigen/Core>

// System include(s).
#include <cmath>
#include <type_traits>

namespace algebra::eigen::math {

/** This method retrieves phi from a vector, vector base with rows >= 2
 *
 * @param v the input vector
 **/
template <
    typename derived_type,
    std::enable_if_t<Eigen::MatrixBase<derived_type>::RowsAtCompileTime >= 2,
                     bool> = true>
ALGEBRA_HOST_DEVICE inline auto phi(
    const Eigen::MatrixBase<derived_type> &v) noexcept {

  return std::atan2(v[1], v[0]);
}

/** This method retrieves theta from a vector, vector base with rows >= 3
 *
 * @param v the input vector
 **/
template <
    typename derived_type,
    std::enable_if_t<Eigen::MatrixBase<derived_type>::RowsAtCompileTime >= 3,
                     bool> = true>
ALGEBRA_HOST_DEVICE inline auto theta(
    const Eigen::MatrixBase<derived_type> &v) noexcept {

  return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
}

/** This method retrieves the perpenticular magnitude of a vector with rows >= 2
 *
 * @param v the input vector
 **/
template <
    typename derived_type,
    std::enable_if_t<Eigen::MatrixBase<derived_type>::RowsAtCompileTime >= 2,
                     bool> = true>
ALGEBRA_HOST_DEVICE inline auto perp(
    const Eigen::MatrixBase<derived_type> &v) noexcept {

  return std::sqrt(v[0] * v[0] + v[1] * v[1]);
}

/** This method retrieves the norm of a vector, no dimension restriction
 *
 * @param v the input vector
 **/
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline auto norm(const Eigen::MatrixBase<derived_type> &v) {

  return v.norm();
}

/** This method retrieves the pseudo-rapidity from a vector or vector base with
 *rows >= 3
 *
 * @param v the input vector
 **/
template <
    typename derived_type,
    std::enable_if_t<Eigen::MatrixBase<derived_type>::RowsAtCompileTime >= 3,
                     bool> = true>
ALGEBRA_HOST_DEVICE inline auto eta(
    const Eigen::MatrixBase<derived_type> &v) noexcept {

  return std::atanh(v[2] / v.norm());
}

}  // namespace algebra::eigen::math
