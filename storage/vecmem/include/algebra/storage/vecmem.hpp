/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "algebra/concepts.hpp"
#include "algebra/storage/impl/cmath_getter.hpp"
#include "algebra/type_traits.hpp"

// VecMem include(s).
#include <vecmem/containers/static_array.hpp>

// System include(s).
#include <cstddef>

namespace algebra {

namespace vecmem {

/// size type for VecMem storage model
using size_type = std::size_t;
/// Value type for VecMem storage model
template <concepts::value T>
using value_type = T;
/// Scalar type for VecMem storage model
template <concepts::value T>
using scalar_type = T;
/// Array type used in the VecMem storage model
template <typename T, std::size_t N>
using storage_type = ::vecmem::static_array<T, N>;
/// Vector type used in the VecMem storage model
template <concepts::scalar T, std::size_t N>
using vector_type = storage_type<T, N>;
/// Matrix type used in the VecMem storage model
template <concepts::scalar T, std::size_t ROWS, std::size_t COLS>
using matrix_type = storage_type<storage_type<T, ROWS>, COLS>;

/// 3-element "vector" type, using @c vecmem::static_array
template <concepts::scalar T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c vecmem::static_array
template <concepts::scalar T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c vecmem::static_array
template <concepts::scalar T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c vecmem::static_array
template <concepts::scalar T>
using point2 = vector2<T>;

/// Element Getter
using element_getter = cmath::storage::element_getter;
/// Block Getter
using block_getter = cmath::storage::block_getter;

}  // namespace vecmem

ALGEBRA_PLUGINS_DEFINE_TYPE_TRAITS(vecmem)

}  // namespace algebra
