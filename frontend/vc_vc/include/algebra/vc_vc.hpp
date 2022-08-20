/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Operators include(s).
#include "algebra/math/cmath_operators.hpp"

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/vc.hpp"
#include "algebra/storage/vc.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

namespace algebra {

using size_type = vc::size_type;
template <typename T, size_type N>
using array_type = vc::storage_type<T, N>;

namespace vector {

template <typename scalar_t>
using actor = cmath::vector::actor<size_type, Vc::array, scalar_t>;

}  // namespace vector

namespace matrix {

template <typename T, size_type ROWS, size_type COLS>
using matrix_type = vc::matrix_type<T, ROWS, COLS>;

template <typename scalar_t>
using element_getter = cmath::element_getter<size_type, Vc::array, scalar_t>;

template <typename scalar_t>
using block_getter = cmath::block_getter<size_type, Vc::array, scalar_t>;

// matrix actor
template <typename scalar_t, typename determinant_actor_t,
          typename inverse_actor_t>
using actor =
    cmath::matrix::actor<size_type, Vc::array, matrix_type, scalar_t,
                         determinant_actor_t, inverse_actor_t,
                         element_getter<scalar_t>, block_getter<scalar_t>>;

namespace determinant {

// determinant aggregation
template <typename scalar_t, class... As>
using actor =
    cmath::matrix::determinant::actor<size_type, matrix_type, scalar_t, As...>;

// determinant::cofactor
template <typename scalar_t, size_type... Ds>
using cofactor =
    cmath::matrix::determinant::cofactor<size_type, matrix_type, scalar_t,
                                         element_getter<scalar_t>, Ds...>;

// determinant::hard_coded
template <typename scalar_t, size_type... Ds>
using hard_coded =
    cmath::matrix::determinant::hard_coded<size_type, matrix_type, scalar_t,
                                           element_getter<scalar_t>, Ds...>;

// preset(s) as standard option(s) for user's convenience
template <typename scalar_t>
using preset0 = actor<scalar_t, cofactor<scalar_t>, hard_coded<scalar_t, 2, 4>>;

}  // namespace determinant

namespace inverse {

// inverion aggregation
template <typename scalar_t, class... As>
using actor =
    cmath::matrix::inverse::actor<size_type, matrix_type, scalar_t, As...>;

// inverse::cofactor
template <typename scalar_t, size_type... Ds>
using cofactor =
    cmath::matrix::inverse::cofactor<size_type, matrix_type, scalar_t,
                                     element_getter<scalar_t>, Ds...>;

// inverse::hard_coded
template <typename scalar_t, size_type... Ds>
using hard_coded =
    cmath::matrix::inverse::hard_coded<size_type, matrix_type, scalar_t,
                                       element_getter<scalar_t>, Ds...>;

// preset(s) as standard option(s) for user's convenience
template <typename scalar_t>
using preset0 = actor<scalar_t, cofactor<scalar_t>, hard_coded<scalar_t, 2, 4>>;

}  // namespace inverse

}  // namespace matrix

namespace vc {

/// @name Vc based transforms on @c algebra::vc::storage_type
/// @{

template <typename T>
using matrix_actor = algebra::matrix::actor<T, matrix::determinant::preset0<T>,
                                            matrix::inverse::preset0<T>>;
template <typename T>
using vector_actor = algebra::vector::actor<T>;
template <typename T>
using transform3 = math::transform3<storage_type, T, matrix_actor<T>,
                                    vector_actor<T>, vector3<T>, point2<T>>;

/// @}

/// @name cmath based track indices

using track_indices = cmath::index::track_indices;

/// @}

/// @name cmath based common algebras
/// @{

template <typename T>
using column_wise_operator =
    common::column_wise_operator<matrix_actor<T>, vector_actor<T>>;

/// @}

}  // namespace vc

namespace getter {

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::vc::transform3<float>
template <std::size_t SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline Vc::array<float, 3> vector(
    const vc::transform3<float>::matrix44& m,
    std::size_t
#ifndef NDEBUG
        row
#endif  // not NDEBUG
    ,
    std::size_t col) {

  assert(row == 0);
  assert(col < 4);
  switch (col) {
    case 0:
      return {m.x[0], m.x[1], m.x[2]};
    case 1:
      return {m.y[0], m.y[1], m.y[2]};
    case 2:
      return {m.z[0], m.z[1], m.z[2]};
    case 3:
      return {m.t[0], m.t[1], m.t[2]};
    default:
      return {m.x[0], m.x[1], m.x[2]};
  }
}

/// Function extracting a slice from the matrix used by
/// @c algebra::vc::transform3<double>
template <std::size_t SIZE, std::enable_if_t<SIZE <= 4, bool> = true>
ALGEBRA_HOST_DEVICE inline Vc::array<double, 3> vector(
    const vc::transform3<double>::matrix44& m,
    std::size_t
#ifndef NDEBUG
        row
#endif  // not NDEBUG
    ,
    std::size_t col) {

  assert(row == 0);
  assert(col < 4);
  switch (col) {
    case 0:
      return {m.x[0], m.x[1], m.x[2]};
    case 1:
      return {m.y[0], m.y[1], m.y[2]};
    case 2:
      return {m.z[0], m.z[1], m.z[2]};
    case 3:
      return {m.t[0], m.t[1], m.t[2]};
    default:
      return {m.x[0], m.x[1], m.x[2]};
  }
}

/// @name Getter functions on @c algebra::vc::matrix_type
/// @{

using cmath::element;

/// @}

}  // namespace getter

}  // namespace algebra
