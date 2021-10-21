/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"

// Eigen include(s).
#include <Eigen/Core>

// System include(s).
#include <cstddef>

namespace algebra {
namespace eigen {

/// Array type used in this storage model
template <typename T, std::size_t N>
using storage_type = Eigen::Matrix<T, N, 1>;

/// 3-element "vector" type, using @c Eigen::Matrix
using vector3 = storage_type<scalar, 3>;
/// Point in 3D space, using @c Eigen::Matrix
using point3 = vector3;
/// Point in 2D space, using @c Eigen::Matrix
using point2 = storage_type<scalar, 2>;

}  // namespace eigen
}  // namespace algebra
