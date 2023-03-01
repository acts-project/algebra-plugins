/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Fastor include(s).
#include <Fastor/Fastor.h>

// System include(s).
#include <cstddef>

namespace algebra::fastor {

/// Fastor matrix type
template <typename T, std::size_t M1, std::size_t N>
class Matrix : public Fastor::Tensor<T, M1, N> {

 public:
  /// Inherit all constructors from the base class
  using Fastor::Tensor<T, M1, N>::Tensor;

  template <typename U, std::size_t M2,
  			  std::enable_if_t<std::is_convertible_v<U, T>, bool> = true
    >
  inline Matrix<T, M1, M2> operator*(const Matrix<U, N, M2>& other) {
  	return Fastor::matmul(static_cast<Fastor::Tensor<T, M1, N>>(*this), static_cast<Fastor::Tensor<T, N, M2>>(other));
  }

  template <typename U, std::size_t M2,
              std::enable_if_t<std::is_convertible_v<U, T>, bool> = true
    >
  inline Matrix<T, M1, M2> operator*(const Fastor::Tensor<U, N, M2>& other) {
  	return Fastor::matmul(static_cast<Fastor::Tensor<T, M1, N>>(*this), static_cast<Fastor::Tensor<T, N, M2>>(other));
  }

  template <typename U,
              std::enable_if_t<std::is_convertible_v<U, T>, bool> = true
    >
  inline auto operator*(const Fastor::Tensor<U, N>& other) {
  	return Fastor::matmul(*this, other);
  }

  template <typename U,
              std::enable_if_t<std::is_convertible_v<U, T>, bool> = true
    >
  inline Fastor::Tensor<T, M1, N> operator*(const U scalar) {
  	return static_cast<T>(scalar) * (*this);
  }

};  // class Matrix 

}  // namespace algebra::fastor
