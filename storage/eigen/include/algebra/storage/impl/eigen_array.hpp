/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Eigen include(s).
#include <Eigen/Core>

// System include(s).
#include <cstddef>

namespace algebra::eigen {

/// Eigen array type
template <typename T, std::size_t N>
class array : public Eigen::Matrix<T, N, 1> {

 public:
  /// Inherit all constructors from the base class
  using Eigen::Matrix<T, N, 1>::Matrix;

};  // class array

}  // namespace algebra::eigen
