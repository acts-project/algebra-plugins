/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/math/cmath.hpp"
#include "algebra/storage/vecmem.hpp"

/// @name Operators on @c algebra::vecmem::storage_type
/// @{

using algebra::cmath::operator*;
using algebra::cmath::operator-;
using algebra::cmath::operator+;

/// @}

namespace algebra {
namespace vecmem {

/// @name cmath based transforms on @c algebra::vecmem::storage_type
/// @{

template <typename T>
using transform3 = cmath::transform3<std::size_t, vecmem::storage_type, T>;
template <typename T>
using cartesian2 = cmath::cartesian2<transform3<T> >;
template <typename T>
using polar2 = cmath::polar2<transform3<T> >;
template <typename T>
using cylindrical2 = cmath::cylindrical2<transform3<T> >;

/// @}

}  // namespace vecmem

namespace getter {

/// @name Getter functions on @c algebra::vecmem::storage_type
/// @{

using cmath::eta;
using cmath::norm;
using cmath::perp;
using cmath::phi;
using cmath::theta;

/// @}

/// Function extracting a slice from the matrix used by
/// @c algebra::vecmem::transform3
template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS,
          typename scalar_t>
ALGEBRA_HOST_DEVICE inline vecmem::storage_type<scalar_t, SIZE> vector(
    const vecmem::storage_type<vecmem::storage_type<scalar_t, ROWS>, COLS>& m,
    std::size_t row, std::size_t col) {

  return cmath::vector_getter<std::size_t, vecmem::storage_type, scalar_t,
                              SIZE>()(m, row, col);
}

/// Function extracting an element from a matrix (const)
template <typename scalar_t, std::size_t ROWS, std::size_t COLS>
ALGEBRA_HOST_DEVICE inline scalar_t element(
    const vecmem::matrix_type<scalar_t, ROWS, COLS>& m, std::size_t row,
    std::size_t col) {

  return cmath::element_getter<std::size_t, vecmem::storage_type, scalar_t>()(
      m, row, col);
}

/// Function extracting an element from a matrix (non-const)
template <typename scalar_t, std::size_t ROWS, std::size_t COLS>
ALGEBRA_HOST_DEVICE inline scalar_t& element(
    vecmem::matrix_type<scalar_t, ROWS, COLS>& m, std::size_t row,
    std::size_t col) {

  return cmath::element_getter<std::size_t, vecmem::storage_type, scalar_t>()(
      m, row, col);
}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vecmem::storage_type
/// @{

using cmath::cross;
using cmath::dot;
using cmath::normalize;

/// @}

}  // namespace vector
}  // namespace algebra
