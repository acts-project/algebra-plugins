/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"
#include "algebra/common/scalar.hpp"

// Vc include(s).
#include <Vc/Vc>

// System include(s).
#include <vector>

namespace algebra {
namespace vc {

/// Array type used in this storage model
template <typename T, std::size_t N>
using storage_type = Vc::SimdArray<T, N>;

/// 3-element "vector" type, using @c Vc::SimdArray
using vector3 = storage_type<scalar, 3>;
/// Point in 3D space, using @c Vc::SimdArray
using point3 = vector3;
/// Point in 2D space, using @c Vc::SimdArray
using point2 = storage_type<scalar, 2>;

}  // namespace vc
}  // namespace algebra
