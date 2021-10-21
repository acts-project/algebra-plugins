/** Algebra plugin library, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/storage/array.hpp"
#include "algebra/math/cmath.hpp"

namespace algebra {
namespace array {

using transform3 = cmath::transform3<storage_type, scalar>;
using cartesian2 = cmath::cartesian2<storage_type, scalar>;
using polar2 = cmath::polar2<storage_type, scalar>;
using cylindrical2 = cmath::cylindrical2<storage_type, scalar>;

}  // namespace array

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<array::storage_type>(a); };
auto theta = [](const auto& a) { return cmath::theta<array::storage_type>(a); };
auto perp = [](const auto& a) { return cmath::perp<array::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<array::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<array::storage_type>(a); };

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) { return cmath::cross<array::storage_type>(a, b); };
auto dot = [](const auto& a, const auto& b) { return cmath::dot<array::storage_type>(a, b); };
auto normalize = [](const auto& a) { return cmath::normalize<array::storage_type>(a); };

}  // namespace vector
}  // namespace algebra

// Test include(s).
#define __plugin algebra::array
#define __plugin_without_matrix_element_accessor 1
#include "test_plugin.inl"
