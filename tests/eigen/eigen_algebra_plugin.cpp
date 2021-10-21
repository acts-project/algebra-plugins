/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/storage/eigen.hpp"
#include "algebra/math/cmath.hpp"

// Eigen include(s).
#include <Eigen/Geometry>

namespace algebra {
namespace eigen {

using transform3 = cmath::transform3<storage_type, scalar, Eigen::Transform<scalar, 3, Eigen::Affine>::MatrixType, algebra::eigen::element_getter, algebra::eigen::block_getter >;
using cartesian2 = cmath::cartesian2<storage_type, scalar, transform3>;
using polar2 = cmath::polar2<storage_type, scalar, transform3>;
using cylindrical2 = cmath::cylindrical2<storage_type, scalar, transform3>;

}  // namespace vc

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<eigen::storage_type>(a); };
auto theta = [](const auto& a) { return cmath::theta<eigen::storage_type>(a); };
auto perp = [](const auto& a) { return cmath::perp<eigen::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<eigen::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<eigen::storage_type>(a); };

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) { return cmath::cross<eigen::storage_type>(a, b); };
auto dot = [](const auto& a, const auto& b) { return cmath::dot<eigen::storage_type>(a, b); };
auto normalize = [](const auto& a) { return cmath::normalize<eigen::storage_type>(a); };

}  // namespace vector
}  // namespace algebra

// Test include(s).
#define __plugin algebra::eigen
#include "test_plugin.inl"
