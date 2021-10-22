/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/vc.hpp"
#include "algebra/storage/vecmem.hpp"

namespace algebra {

using cmath::operator*;
using cmath::operator-;
using cmath::operator+;

namespace vecmem {

using transform3 = vc::math::transform3<vecmem::storage_type, scalar>;
using cartesian2 = vc::math::cartesian2<transform3>;
using polar2 = vc::math::polar2<vecmem::storage_type, transform3>;
using cylindrical2 = vc::math::cylindrical2<vecmem::storage_type, transform3>;

}  // namespace vecmem

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<vecmem::storage_type>(a); };
auto theta = [](const auto& a) {
  return cmath::theta<vecmem::storage_type>(a);
};
auto perp = [](const auto& a) { return cmath::perp<vecmem::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<vecmem::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<vecmem::storage_type>(a); };

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) {
  return cmath::cross<vecmem::storage_type>(a, b);
};
auto dot = [](const auto& a, const auto& b) {
  return cmath::dot<vecmem::storage_type>(a, b);
};
auto normalize = [](const auto& a) {
  return cmath::normalize<vecmem::storage_type>(a);
};

}  // namespace vector
}  // namespace algebra
