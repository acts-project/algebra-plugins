/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

// System include(s).

namespace algebra::fastor::math {

/** Get a normalized version of the input vector
 *
 * @param v the input vector
 **/
template <typename scalar_t, auto N>
ALGEBRA_HOST inline Fastor::Tensor<scalar_t, N> normalize(
    const Fastor::Tensor<scalar_t, N> &v) {

  return Fastor::norm(v);
}

}  // namespace algebra::fastor::math
