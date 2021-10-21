/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/storage/eigen.hpp"
#include "algebra/math/cmath.hpp"

// Eigen include(s).
#include <Eigen/Geometry>

namespace algebra {
namespace eigen {

/// Functor used to access elements of Eigen matrices
struct element_getter {
  /// Get non-const access to a matrix element
  template <typename matrix_type>
  ALGEBRA_HOST_DEVICE inline scalar&
  operator()(matrix_type &m, unsigned int row, unsigned int col) const {

    return m(row, col);
  }
  /// Get const access to a matrix element
  template <typename matrix_type>
  ALGEBRA_HOST_DEVICE inline scalar
  operator()(const matrix_type &m, unsigned int row, unsigned int col) const {

    return m(row, col);
  }
};  // struct element_getter

/// Functor used to extract a slice from Eigen matrices
template <std::size_t kROWS>
struct vector_getter {
  template <typename matrix_type>
  ALGEBRA_HOST_DEVICE
  storage_type<scalar, kROWS> operator()(const matrix_type &m,
                                         std::size_t row, std::size_t col) const {

    return m.template block<kROWS, 1>(row, col);
  }
};  // struct vector_getter

/// Functor used to extract a block from Eigen matrices
struct block_getter {
  template <std::size_t kROWS, std::size_t kCOLS, typename matrix_type>
  ALGEBRA_HOST_DEVICE
  Eigen::Matrix<scalar, kROWS, kCOLS> operator()(const matrix_type &m,
                                                 std::size_t row, std::size_t col) const {

    return m.template block<kROWS, kCOLS>(row, col);
  }
};  // struct block_getter

using transform3 = cmath::transform3<storage_type, scalar, Eigen::Transform<scalar, 3, Eigen::Affine>::MatrixType, algebra::eigen::element_getter, algebra::eigen::block_getter >;
using cartesian2 = cmath::cartesian2<storage_type, scalar, transform3>;
using polar2 = cmath::polar2<storage_type, scalar, transform3>;
using cylindrical2 = cmath::cylindrical2<storage_type, scalar, transform3>;

}  // namespace vc

namespace getter {

auto phi = [](const auto& a) { return cmath::phi<eigen::storage_type>(a); };
auto theta = [](const auto& a) { return cmath::theta<eigen::storage_type>(a); };
auto perp = [](const auto& a) { return cmath::perp<eigen::storage_type>(a); };
auto norm = [](const auto& a) { return cmath::norm<eigen::storage_type>(a); };
auto eta = [](const auto& a) { return cmath::eta<eigen::storage_type>(a); };

}  // namespace getter

namespace vector {

auto cross = [](const auto& a, const auto& b) { return cmath::cross<eigen::storage_type>(a, b); };
auto dot = [](const auto& a, const auto& b) { return cmath::dot<eigen::storage_type>(a, b); };
auto normalize = [](const auto& a) { return cmath::normalize<eigen::storage_type>(a); };

}  // namespace vector
}  // namespace algebra
