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

auto phi = [](const auto& a) { return smatrix::math::phi(a); };
auto theta = [](const auto& a) { return smatrix::math::theta(a); };
auto perp = [](const auto& a) { return smatrix::math::perp(a); };
auto norm = [](const auto& a) { return smatrix::math::norm(a); };
auto eta = [](const auto& a) { return smatrix::math::eta(a); };

template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS>
ALGEBRA_HOST_DEVICE inline auto vector(
    const ROOT::Math::SMatrix<scalar, ROWS, COLS>& m, unsigned int row,
    unsigned int col) {

  return m.template SubCol<smatrix::storage_type<scalar, SIZE> >(col, row);
}

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) {
  return smatrix::math::cross(a, b);
};
auto dot = [](const auto& a, const auto& b) {
  return smatrix::math::dot(a, b);
};
auto normalize = [](const auto& a) { return smatrix::math::normalize(a); };

}  // namespace vector
}  // namespace algebra
