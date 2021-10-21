/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/common/algebra_qualifiers.hpp"
#include "algebra/math/impl/eigen_transform3.hpp"

// Eigen include(s).
#include <Eigen/Core>

namespace algebra::eigen::math {

/** Local frame projection into a cartesian coordinate frame
 */
template <typename transform3_t>
struct cartesian2 {

  /// @name Type definitions for the struct
  /// @{

  /// Transformation matching this struct
  using transform3_type = transform3_t;

  /// Point in 2D space
  using point2 = typename transform3_type::point2;
  /// Point in 3D space
  using point3 = typename transform3_type::point3;

  /// @}

  /** This method transform from a point from the global 3D cartesian frame to
   * the local 2D cartesian frame
   *
   * @param v the point in local frame
   *
   * @return a local point2
   */
  template <typename derived_type>
  ALGEBRA_HOST_DEVICE inline point2 operator()(
      const Eigen::MatrixBase<derived_type> &v) const {
    constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
    constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
    static_assert(rows == 3 and cols == 1,
                  "transform::point3_topoint2(v) requires a (3,1) matrix");
    return (v.template segment<2>(0)).eval();
  }

  /** This method transform from a point from the global 3D cartesian frame to
   *the local 2D cartesian frame
   *
   * @param trf the transform from global to local thredimensional frame
   * @param p the point in global frame
   *
   * @return a local point2
   **/
  ALGEBRA_HOST_DEVICE
  inline point2 operator()(const transform3_type &trf, const point3 &p) const {
    return operator()(trf.point_to_local(p));
  }
};  // struct cartesian2

}  // namespace algebra::eigen::math
