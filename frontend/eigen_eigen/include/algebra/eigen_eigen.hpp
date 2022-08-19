/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "algebra/math/eigen.hpp"
#include "algebra/storage/eigen.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace algebra {

namespace getter {

/// Function extracting a slice from the matrix used by
/// @c algebra::eigen::transform3
template <unsigned int SIZE, typename derived_type>
ALGEBRA_HOST_DEVICE inline auto vector(const Eigen::MatrixBase<derived_type>& m,
                                       std::size_t row, std::size_t col) {

  return m.template block<SIZE, 1>(row, col);
}

/// @name Getter functions on @c algebra::eigen::matrix_type
/// @{

using eigen::math::element;

/// @}

}  // namespace getter

namespace vector {

using actor = eigen::vector::actor;

}  // namespace vector

namespace matrix {

template <typename scalar_t>
using actor = eigen::matrix::actor<scalar_t>;

}  // namespace matrix

namespace eigen {

/// @name Eigen based transforms on @c algebra::eigen::storage_type
/// @{

template <typename T>
using transform3 =
    math::transform3<T, algebra::matrix::actor<T>, algebra::vector::actor>;

/// @}

/// @name Eigen based track indices

using track_indices = eigen::index::track_indices;

/// @}

/// @name Eigen based common algebras
/// @{

template <typename T>
using cartesian2 = common::cartesian2<transform3<T>, track_indices>;
template <typename T>
using cartesian3 = common::cartesian3<transform3<T>, track_indices>;
template <typename T>
using polar2 = common::polar2<transform3<T>, track_indices>;
template <typename T>
using cylindrical2 = common::cylindrical2<transform3<T>, track_indices>;
template <typename T>
using cylindrical3 = common::cylindrical3<transform3<T>, track_indices>;
template <typename T>
using line2 = common::line2<transform3<T>, track_indices>;

template <typename T>
using bound_track_parameters =
    common::bound_track_parameters<algebra::matrix::actor<T>,
                                   algebra::vector::actor, track_indices>;
template <typename T>
using free_track_parameters =
    common::free_track_parameters<algebra::matrix::actor<T>,
                                  algebra::vector::actor, track_indices>;

template <typename T>
using column_wise_operator =
    common::column_wise_operator<algebra::matrix::actor<T>,
                                 algebra::vector::actor>;

/// @}

}  // namespace eigen

}  // namespace algebra
