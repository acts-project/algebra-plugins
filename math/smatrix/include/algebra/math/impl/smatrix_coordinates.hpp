/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/coordinates.hpp"

namespace algebra::smatrix::coordinate {

template <typename transform3_t>
using cartesian2 = algebra::common::cartesian2<transform3_t>;

template <typename transform3_t>
using polar2 = algebra::common::polar2<transform3_t>;

template <typename transform3_t>
using cylindrical2 = algebra::common::cylindrical2<transform3_t>;

}  // namespace algebra::smatrix::coordinate