/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

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
