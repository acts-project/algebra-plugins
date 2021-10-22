/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"
#include "algebra/storage/smatrix.hpp"

// ROOT/Smatrix include(s).
#include <Math/Functions.h>

namespace algebra::smatrix::math {

/** Get a normalized version of the input vector
 *
 * @tparam derived_type is the matrix template
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST inline storage_type<scalar_t, N> normalize(
    const storage_type<scalar_t, N> &v) {

  return ROOT::Math::Unit(v);
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
template <typename scalar_t, auto N>
ALGEBRA_HOST inline scalar_t dot(const storage_type<scalar_t, N> &a,
                                 const storage_type<scalar_t, N> &b) {

  return ROOT::Math::Dot(a, b);
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
template <typename scalar_t, auto N>
ALGEBRA_HOST inline storage_type<scalar_t, N> cross(
    const storage_type<scalar_t, N> &a, const storage_type<scalar_t, N> &b) {

  return ROOT::Math::Cross(a, b);
}

}  // namespace algebra::smatrix::math
