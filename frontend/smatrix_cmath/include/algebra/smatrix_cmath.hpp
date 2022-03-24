/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/math/smatrix.hpp"
#include "algebra/storage/smatrix.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

namespace algebra {
namespace smatrix {

/// @name cmath based transforms on @c algebra::smatrix::storage_type
/// @{

template <typename T>
using transform3 =
    cmath::transform3<unsigned int, smatrix::storage_type, T,
                      ROOT::Math::SMatrix<T, 4, 4>, math::element_getter<T>,
                      math::block_getter<T>>;
template <typename T>
using cartesian2 = cmath::cartesian2<transform3<T>>;
template <typename T>
using polar2 = cmath::polar2<transform3<T>>;
template <typename T>
using cylindrical2 = cmath::cylindrical2<transform3<T>>;

/// @}

}  // namespace smatrix

namespace getter {

/// @name Getter functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::eta;
using smatrix::math::norm;
using smatrix::math::perp;
using smatrix::math::phi;
using smatrix::math::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::smatrix::transform3
template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS,
          typename scalar_t>
ALGEBRA_HOST_DEVICE inline auto vector(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS>& m, unsigned int row,
    unsigned int col) {

  return m.template SubCol<smatrix::storage_type<scalar_t, SIZE>>(col, row);
}

/// @name Getter functions on @c algebra::smatrix::matrix_type
/// @{

using smatrix::math::element;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::smatrix::storage_type
/// @{

using smatrix::math::cross;
using smatrix::math::dot;
using smatrix::math::normalize;

/// @}

}  // namespace vector

namespace matrix {

template <typename T, unsigned int ROWS, unsigned int COLS>
using matrix_type = smatrix::matrix_type<T, ROWS, COLS>;
template <typename scalar_t>
using element_getter_type = smatrix::math::element_getter<scalar_t>;

// matrix actor
template <typename size_type, typename scalar_t, typename determinant_actor_t,
          typename inverse_actor_t>
using actor =
    cmath::matrix::actor<size_type, matrix_type, scalar_t, determinant_actor_t,
                         inverse_actor_t, element_getter_type<scalar_t>>;

namespace determinant {

// determinant aggregation
template <typename size_type, typename scalar_t, class... As>
using actor =
    cmath::matrix::determinant::actor<size_type, matrix_type, scalar_t, As...>;

// determinant::cofactor
template <typename size_type, typename scalar_t, size_type... Ds>
using cofactor =
    cmath::matrix::determinant::cofactor<size_type, matrix_type, scalar_t,
                                         element_getter_type<scalar_t>, Ds...>;

// determinant::hard_coded
template <typename size_type, typename scalar_t, size_type... Ds>
using hard_coded = cmath::matrix::determinant::hard_coded<
    size_type, matrix_type, scalar_t, element_getter_type<scalar_t>, Ds...>;

// preset(s) as standard option(s) for user's convenience
template <typename size_type, typename scalar_t>
using preset0 = actor<size_type, scalar_t, cofactor<size_type, scalar_t>,
                      hard_coded<size_type, scalar_t, 2>>;

}  // namespace determinant

namespace inverse {

// inverion aggregation
template <typename size_type, typename scalar_t, class... As>
using actor =
    cmath::matrix::inverse::actor<size_type, matrix_type, scalar_t, As...>;

// inverse::cofactor
template <typename size_type, typename scalar_t, size_type... Ds>
using cofactor =
    cmath::matrix::inverse::cofactor<size_type, matrix_type, scalar_t,
                                     element_getter_type<scalar_t>, Ds...>;

// inverse::hard_coded
template <typename size_type, typename scalar_t, size_type... Ds>
using hard_coded =
    cmath::matrix::inverse::hard_coded<size_type, matrix_type, scalar_t,
                                       element_getter_type<scalar_t>, Ds...>;

// preset(s) as standard option(s) for user's convenience
template <typename size_type, typename scalar_t>
using preset0 = actor<size_type, scalar_t, cofactor<size_type, scalar_t>,
                      hard_coded<size_type, scalar_t, 2>>;

}  // namespace inverse

}  // namespace matrix

}  // namespace algebra
