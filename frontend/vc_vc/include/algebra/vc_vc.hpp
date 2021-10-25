/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"
#include "algebra/math/vc.hpp"
#include "algebra/storage/vc.hpp"

namespace algebra {
namespace vc {

using transform3 = math::transform3<storage_type, scalar>;
using cartesian2 = math::cartesian2<transform3>;
using polar2 = math::polar2<storage_type, transform3>;
using cylindrical2 = math::cylindrical2<storage_type, transform3>;

}  // namespace vc

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<vc::storage_type>(a); };
auto theta = [](const auto& a) { return cmath::theta<vc::storage_type>(a); };
auto perp = [](const auto& a) { return cmath::perp<vc::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<vc::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<vc::storage_type>(a); };

template <auto SIZE, typename input_matrix_type>
ALGEBRA_HOST_DEVICE inline vc::storage_type<scalar, 4> vector(
    const input_matrix_type& m, std::size_t row, std::size_t col) {

  switch (row) {
    case 0:
      return m.x;
    case 1:
      return m.y;
    case 2:
      return m.z;
    case 3:
      return m.t;
    default:
      return {};
  }
}

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) {
  return vc::math::cross<vc::storage_type>(a, b);
};
auto dot = [](const auto& a, const auto& b) {
  return vc::math::dot<vc::storage_type>(a, b);
};
auto normalize = [](const auto& a) {
  return vc::math::normalize<vc::storage_type>(a);
};

}  // namespace vector
}  // namespace algebra
