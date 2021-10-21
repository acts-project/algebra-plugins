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

namespace getter {

/** This method retrieves a column from a matrix
 *
 * @tparam matrix_type generic input matrix type
 *
 * @param m the input matrix
 **/
template <unsigned int kROWS, typename matrix_type>
ALGEBRA_HOST_DEVICE inline vc::storage_type<scalar, kROWS> vector(
    const matrix_type &m, unsigned int row, unsigned int col) noexcept {

  vc::storage_type<scalar, kROWS> subvector;
  for (unsigned int irow = row; irow < row + kROWS; ++irow) {
    subvector[irow - row] = m[col][irow];
  }
  return subvector;
}

/** This method retrieves a submatrix
 *
 * @param m the input matrix
 **/
template <unsigned int kROWS, unsigned int kCOLS, typename matrix_type>
ALGEBRA_HOST_DEVICE inline auto block(const matrix_type &m, unsigned int row,
                                      unsigned int col) noexcept {

  Vc::array<vc::storage_type<scalar, kROWS>, kCOLS> submatrix;

  // vc_array::storage_type<scalar, 16> submatrix;

  for (unsigned int icol = col; icol < col + kCOLS; ++icol) {
    for (unsigned int irow = row; irow < row + kROWS; ++irow) {
      submatrix[icol - col][irow - row] = m[icol][irow];
    }
  }
  return submatrix;


  // std::vector<decltype(m.x), Vc::Allocator<decltype(m.x)> >
  //     vecs = {std::move(m.x), std::move(m.y),
  //                                              std::move(m.z), std::move(m.t)};

  // for (unsigned int irow = row; irow < row + kROWS; ++irow) {
  //   for (unsigned int icol = col; icol < col + kCOLS; ++icol) {
  //     submatrix[(icol - col) * (irow - row) + (icol - col)] =
  //         m[irow /*- row*/][icol /*- col*/];
  //   }
  // }
  // return submatrix;
}

}  // namespace getter
}  // namespace algebra
