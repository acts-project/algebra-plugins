/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/impl/eigen_array.hpp"

// System include(s).
#include <cstddef>

namespace algebra::eigen {

/// Array type used in the Eigen storage model
template <typename T, int N>
using storage_type = array<T, N>;
/// Matrix type used in the Eigen storage model
template <typename T, int ROWS, int COLS>
using matrix_type = Eigen::Matrix<T, ROWS, COLS, 0, ROWS, COLS>;

/// 3-element "vector" type, using @c algebra::eigen::array
template <typename T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c algebra::eigen::array
template <typename T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c algebra::eigen::array
template <typename T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c algebra::eigen::array
template <typename T>
using point2 = vector2<T>;

}  // namespace algebra::eigen
