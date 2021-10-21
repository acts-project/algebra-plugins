/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/storage/scalar.hpp"

// VecMem include(s).
#include <vecmem/container/static_array.hpp>

/// Main algebra namespace
namespace algebra {
/// @c vecmem::static_array algebra definitions
namespace vecmem_array {

/// Array type used in this storage model
template <typename T, std::size_t N>
using storage_type = vecmem::static_array<T, N>;

/// 3-element "vector" type, using @c std::array
using vector3 = storage_type<scalar, 3>;
/// Point in 3D space, using @c std::array
using point3 = vector3;
/// Point in 2D space, using @c std::array
using point2 = storage_type<scalar, 2>;

}  // namespace vecmem_array

namespace getter {

/** This method retrieves a column from a matrix
 *
 * @param m the input matrix
 **/
template <unsigned int kROWS, typename matrix_type>
ALGEBRA_HOST_DEVICE
inline vecmem_array::storage_type<scalar, kROWS> vector(const matrix_type &m, unsigned int row,
                                                     unsigned int col) noexcept {

  vecmem_array::storage_type<scalar, kROWS> subvector;
  for (unsigned int irow = row; irow < row + kROWS; ++irow) {
    subvector[irow - row] = m[col][irow];
  }
  return subvector;
}

/** This method retrieves a column from a matrix
 *
 * @param m the input matrix
 **/
template <unsigned int kROWS, unsigned int kCOLS, typename matrix_type>
ALGEBRA_HOST_DEVICE inline vecmem_array::storage_type<vecmem_array::storage_type<scalar, kROWS>, kCOLS>
block(const matrix_type &m, unsigned int row, unsigned int col) noexcept {

  vecmem_array::storage_type<vecmem_array::storage_type<scalar, kROWS>, kCOLS> submatrix;
  for (unsigned int icol = col; icol < col + kCOLS; ++icol) {
    for (unsigned int irow = row; irow < row + kROWS; ++irow) {
      submatrix[icol - col][irow - row] = m[icol][irow];
    }
  }
  return submatrix;
}

}  // namespace getter
}  // namespace algebra
