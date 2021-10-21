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

// Eigen include(s).
#include <Eigen/Core>
#include <Eigen/Dense>

// System include(s).
#include <cstddef>

namespace algebra {
namespace eigen {

/// Array type used in this storage model
template <typename T, std::size_t N>
using storage_type = Eigen::Matrix<T, N, 1>;

/// 3-element "vector" type, using @c Eigen::Matrix
using vector3 = storage_type<scalar, 3>;
/// Point in 3D space, using @c Eigen::Matrix
using point3 = vector3;
/// Point in 2D space, using @c Eigen::Matrix
using point2 = storage_type<scalar, 2>;

/**
 * Method retrieving a single element from a matrix
 */
struct element_getter {
  template <typename matrix_type>
  ALGEBRA_HOST_DEVICE inline scalar&
  operator()(matrix_type &m, unsigned int row, unsigned int col) const {

    return m(row, col);
  }
  template <typename matrix_type>
  ALGEBRA_HOST_DEVICE inline scalar
  operator()(const matrix_type &m, unsigned int row, unsigned int col) const {

    return m(row, col);
  }
};  // struct element_getter

/** This method retrieves a column from a matrix
 *
 * @param m the input matrix
 **/
template <std::size_t kROWS>
struct vector_getter {
  template <typename matrix_type>
  ALGEBRA_HOST_DEVICE
  storage_type<scalar, kROWS> operator()(const matrix_type &m,
                                         std::size_t row, std::size_t col) const {

    return m.template block<kROWS, 1>(row, col);
  }
};  // struct vector_getter

/** This method retrieves a column from a matrix
 *
 * @param m the input matrix
 **/
struct block_getter {
  template <std::size_t kROWS, std::size_t kCOLS, typename matrix_type>
  ALGEBRA_HOST_DEVICE
  Eigen::Matrix<scalar, kROWS, kCOLS> operator()(const matrix_type &m,
                                                 std::size_t row, std::size_t col) const {

    return m.template block<kROWS, kCOLS>(row, col);
  }
};  // struct block_getter

}  // namespace eigen
}  // namespace algebra
