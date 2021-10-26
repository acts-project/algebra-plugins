/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"

// Eigen include(s).
#include <Eigen/Core>

namespace algebra::eigen::math {

/** Get a normalized version of the input vector
 *
 * @tparam derived_type is the matrix template
 *
 * @param v the input vector
 **/
template <typename derived_type>
ALGEBRA_HOST_DEVICE inline auto normalize(
    const Eigen::MatrixBase<derived_type> &v) {
  return v.normalized();
}

/** Dot product between two input vectors
 *
 * @tparam derived_type_lhs is the first matrix (epresseion) template
 * @tparam derived_type_rhs is the second matrix (epresseion) template
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return the scalar dot product value
 **/
template <typename derived_type_lhs, typename derived_type_rhs>
ALGEBRA_HOST_DEVICE inline auto dot(
    const Eigen::MatrixBase<derived_type_lhs> &a,
    const Eigen::MatrixBase<derived_type_rhs> &b) {
  return a.dot(b);
}

/** Cross product between two input vectors
 *
 * @tparam derived_type_lhs is the first matrix (epresseion) template
 * @tparam derived_type_rhs is the second matrix (epresseion) template
 *
 * @param a the first input vector
 * @param b the second input vector
 *
 * @return a vector (expression) representing the cross product
 **/
template <typename derived_type_lhs, typename derived_type_rhs>
ALGEBRA_HOST_DEVICE inline auto cross(
    const Eigen::MatrixBase<derived_type_lhs> &a,
    const Eigen::MatrixBase<derived_type_rhs> &b) {
  return a.cross(b);
}

}  // namespace algebra::eigen::math
