/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"
#include "algebra/math/smatrix.hpp"
#include "algebra/storage/smatrix.hpp"

namespace algebra {
namespace smatrix {

using transform3 = math::transform3<scalar>;
using cartesian2 = math::cartesian2<transform3>;
using polar2 = math::polar2<transform3>;
using cylindrical2 = math::cylindrical2<transform3>;

}  // namespace smatrix

namespace getter {

using smatrix::math::eta;
using smatrix::math::norm;
using smatrix::math::perp;
using smatrix::math::phi;
using smatrix::math::theta;

template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline auto vector(
    const ROOT::Math::SMatrix<scalar, ROWS, COLS>& m, unsigned int row,
    unsigned int col) {

  return m.template SubCol<smatrix::storage_type<scalar, SIZE> >(col, row);
}

}  // namespace getter

namespace vector {

using smatrix::math::cross;
using smatrix::math::dot;
using smatrix::math::normalize;

}  // namespace vector
}  // namespace algebra
