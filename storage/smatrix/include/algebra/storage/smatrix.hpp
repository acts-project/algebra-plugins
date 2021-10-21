/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/scalar.hpp"

// ROOT/Smatrix include(s).
#include <Math/SVector.h>

// System include(s).
#include <cstddef>

namespace algebra::smatrix {

/// Array type used in the Smatrix storage model
template <typename T, std::size_t N>
using storage_type = ROOT::Math::SVector<T, N>;

/// 3-element "vector" type, using @c ROOT::Math::SVector
using vector3 = storage_type<scalar, 3>;
/// Point in 3D space, using @c ROOT::Math::SVector
using point3 = vector3;
/// 2-element "vector" type, using @c ROOT::Math::SVector
using vector2 = storage_type<scalar, 2>;
/// Point in 2D space, using @c ROOT::Math::SVector
using point2 = vector2;

}  // namespace algebra::smatrix
