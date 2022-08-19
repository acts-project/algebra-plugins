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
#include "algebra/storage/array.hpp"

namespace algebra {

namespace getter {

/// Function extracting a slice from a matrix
template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS,
          typename scalar_t>
ALGEBRA_HOST_DEVICE inline array::storage_type<scalar_t, SIZE> vector(
    const array::matrix_type<scalar_t, ROWS, COLS>& m, std::size_t row,
    std::size_t col) {

  return cmath::vector_getter<std::size_t, array::storage_type, scalar_t,
                              SIZE>()(m, row, col);
}

/// @name Getter functions on @c algebra::array::matrix_type
/// @{

using cmath::element;

/// @}

}  // namespace getter

using size_type = array::size_type;

template <typename T, size_type N>
using array_type = array::storage_type<T, N>;

namespace vector {

template <typename scalar_t>
using actor = cmath::vector::actor<size_type, array_type, scalar_t>;

}  // namespace vector

namespace matrix {

template <typename T, size_type ROWS, size_type COLS>
using matrix_type = array::matrix_type<T, ROWS, COLS>;

template <typename scalar_t>
using element_getter = cmath::element_getter<size_type, array_type, scalar_t>;

template <typename scalar_t>
using block_getter = cmath::block_getter<size_type, array_type, scalar_t>;

// matrix actor
template <typename scalar_t, typename determinant_actor_t,
          typename inverse_actor_t>
using actor =
    cmath::matrix::actor<size_type, array_type, matrix_type, scalar_t,
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

namespace array {

/// @name cmath based transforms on @c algebra::matrix::actor
/// @{

template <typename T>
using matrix_actor = matrix::actor<T, matrix::determinant::preset0<T>,
                                   matrix::inverse::preset0<T>>;
template <typename T>
using vector_actor = vector::actor<T>;
template <typename T>
using transform3 = cmath::transform3<matrix_actor<T>, vector_actor<T>>;

/// @}

/// @name cmath based track indices

using track_indices = cmath::index::track_indices;

/// @}

/// @name cmath based common algebras
/// @{

template <typename T>
using cartesian2 = common::cartesian2<transform3<T>, track_indices>;
template <typename T>
using cartesian3 = common::cartesian3<transform3<T>, track_indices>;
template <typename T>
using polar2 = common::polar2<transform3<T>, track_indices>;
template <typename T>
using cylindrical2 = common::cylindrical2<transform3<T>, track_indices>;
template <typename T>
using cylindrical3 = common::cylindrical3<transform3<T>, track_indices>;
template <typename T>
using line2 = common::line2<transform3<T>, track_indices>;

template <typename T>
using bound_track_parameters =
    common::bound_track_parameters<matrix_actor<T>, vector_actor<T>,
                                   track_indices>;
template <typename T>
using free_track_parameters =
    common::free_track_parameters<matrix_actor<T>, vector_actor<T>,
                                  track_indices>;

template <typename T>
using column_wise_operator =
    common::column_wise_operator<matrix_actor<T>, vector_actor<T>>;

/// @}

}  // namespace array

}  // namespace algebra
