/** Algebra plugins, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/impl/fastor_matrix.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

// System include(s).
#include <cstddef>  // for the std::size_t type
#include <type_traits>

namespace algebra::fastor::storage {

/// Functor used to access elements of Fastor matrices
struct element_getter {

  template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE inline scalar_t &operator()(
      Fastor::Tensor<scalar_t, ROWS, COLS> &m, std::size_t row,
      std::size_t col) const {

    assert(row < ROWS);
    assert(col < COLS);
    return m(row, col);
  }

  template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const Fastor::Tensor<scalar_t, ROWS, COLS> &m, std::size_t row,
      std::size_t col) const {

    assert(row < ROWS);
    assert(col < COLS);
    return m(row, col);
  }

  template <std::size_t N, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE inline scalar_t &operator()(
      Fastor::Tensor<scalar_t, N, 1> &m, std::size_t row) const {

    assert(row < N);
    return m(row);
  }

  template <std::size_t N, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const Fastor::Tensor<scalar_t, N, 1> &m, std::size_t row) const {

    assert(row < N);
    return m(row);
  }

  template <std::size_t N, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE inline scalar_t &operator()(
      Fastor::Tensor<scalar_t, N> &m, std::size_t row) const {

    assert(row < N);
    return m(row);
  }

  template <std::size_t N, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE inline scalar_t operator()(
      const Fastor::Tensor<scalar_t, N> &m, std::size_t row) const {

    assert(row < N);
    return m(row);
  }

};  // element_getter

/// Function extracting an element from a matrix (const)
template <concepts::scalar scalar_t, std::size_t ROWS, std::size_t COLS>
ALGEBRA_HOST_DEVICE inline scalar_t element(
    const Fastor::Tensor<scalar_t, ROWS, COLS> &m, std::size_t row,
    std::size_t col) {
  return element_getter()(m, row, col);
}

/// Function extracting an element from a matrix (non-const)
template <concepts::scalar scalar_t, std::size_t ROWS, std::size_t COLS>
ALGEBRA_HOST_DEVICE inline scalar_t &element(
    Fastor::Tensor<scalar_t, ROWS, COLS> &m, std::size_t row, std::size_t col) {
  return element_getter()(m, row, col);
}

/// Function extracting an element from a matrix (const)
template <concepts::scalar scalar_t, std::size_t N>
ALGEBRA_HOST_DEVICE inline scalar_t element(
    const Fastor::Tensor<scalar_t, N, 1> &m, std::size_t row) {
  return element_getter()(m, row);
}

/// Function extracting an element from a matrix (non-const)
template <concepts::scalar scalar_t, std::size_t N>
ALGEBRA_HOST_DEVICE inline scalar_t &element(Fastor::Tensor<scalar_t, N, 1> &m,
                                             std::size_t row) {
  return element_getter()(m, row);
}

/// Function extracting an element from a vector (const)
template <concepts::scalar scalar_t, std::size_t N>
ALGEBRA_HOST_DEVICE inline scalar_t element(
    const Fastor::Tensor<scalar_t, N> &m, std::size_t row) {
  return element_getter()(m, row);
}

/// Function extracting an element from a vector (non-const)
template <concepts::scalar scalar_t, std::size_t N>
ALGEBRA_HOST_DEVICE inline scalar_t &element(Fastor::Tensor<scalar_t, N> &m,
                                             std::size_t row) {
  return element_getter()(m, row);
}

/// Functor used to extract a block from Fastor matrices
struct block_getter {

  template <std::size_t ROWS, std::size_t COLS, std::size_t oROWS,
            std::size_t oCOLS, concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE algebra::fastor::Matrix<scalar_t, ROWS, COLS> operator()(
      const Fastor::Tensor<scalar_t, oROWS, oCOLS> &m, std::size_t row,
      std::size_t col) const {

    return m(Fastor::seq(static_cast<int>(row), static_cast<int>(row + ROWS)),
             Fastor::seq(static_cast<int>(col), static_cast<int>(col + COLS)));
  }

  template <std::size_t SIZE, std::size_t oROWS, std::size_t oCOLS,
            concepts::scalar scalar_t>
  ALGEBRA_HOST_DEVICE Fastor::Tensor<scalar_t, SIZE> vector(
      const Fastor::Tensor<scalar_t, oROWS, oCOLS> &m, std::size_t row,
      std::size_t col) const {

    return Fastor::Tensor<scalar_t, SIZE>(
        m(Fastor::seq(static_cast<int>(row), static_cast<int>(row + SIZE)),
          static_cast<int>(col)));
  }

};  // struct block_getter

/// Operator getting a block of a const matrix
template <std::size_t ROWS, std::size_t COLS, std::size_t oROWS,
          std::size_t oCOLS, concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE decltype(auto) block(
    const Fastor::Tensor<scalar_t, oROWS, oCOLS> &m, std::size_t row,
    std::size_t col) {
  return block_getter{}.template operator()<ROWS, COLS>(m, row, col);
}

/// Function extracting a slice from the matrix
template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS,
          concepts::scalar scalar_t>
ALGEBRA_HOST_DEVICE inline decltype(auto) vector(
    const Fastor::Tensor<scalar_t, ROWS, COLS> &m, std::size_t row,
    std::size_t col) {

  return block_getter{}.template vector<SIZE>(m, row, col);
}

/// Operator setting a block with a matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t,
          concepts::matrix input_matrix_type>
ALGEBRA_HOST_DEVICE void set_block(
    input_matrix_type &m, const Fastor::Tensor<scalar_t, ROWS, COLS> &b,
    std::size_t row, std::size_t col) {
  m(Fastor::seq(static_cast<int>(row), static_cast<int>(row + ROWS)),
    Fastor::seq(static_cast<int>(col), static_cast<int>(col + COLS))) = b;
}

/// Operator setting a block with a vector
template <std::size_t ROWS, concepts::scalar scalar_t,
          concepts::matrix input_matrix_type>
ALGEBRA_HOST_DEVICE void set_block(input_matrix_type &m,
                                   const Fastor::Tensor<scalar_t, ROWS> &b,
                                   std::size_t row, std::size_t col) {
  m(Fastor::seq(static_cast<int>(row), static_cast<int>(row + ROWS)),
    static_cast<int>(col)) = b;
}

}  // namespace algebra::fastor::storage
