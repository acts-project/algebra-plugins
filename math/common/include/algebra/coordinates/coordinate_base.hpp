/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"
#include "algebra/tracks/detail/track_helper.hpp"

namespace algebra::common {

/** Coordinate base struct
 */
template <template <class, class> class Derived, typename transform3_t,
          typename track_indices_t>
struct coordinate_base {

  /// @name Type definitions for the struct
  /// @{

  /// Transformation matching this struct
  using transform3_type = transform3_t;
  /// Scalar type
  using scalar_type = typename transform3_type::scalar_type;
  /// Point in 2D space
  using point2 = typename transform3_type::point2;
  /// Point in 3D space
  using point3 = typename transform3_type::point3;
  /// Vector in 3D space
  using vector3 = typename transform3_type::vector3;
  /// Vector actor
  using vector_actor = typename transform3_type::vector_actor;
  /// Matrix actor
  using matrix_actor = typename transform3_type::matrix_actor;
  /// Matrix size type
  using size_type = typename matrix_actor::size_ty;
  /// 2D matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = typename matrix_actor::matrix_type<ROWS, COLS>;

  /// Track indices
  using E = track_indices_t;

  /// Shorthand vector/matrix types related to bound track parameters.
  using bound_vector = matrix_type<E::bound_size, 1>;
  using bound_matrix = matrix_type<E::bound_size, E::bound_size>;

  /// Mapping from bound track parameters.
  using bound_to_free_matrix = matrix_type<E::free_size, E::bound_size>;

  // Shorthand vector/matrix types related to free track parameters.
  using free_vector = matrix_type<E::free_size, 1>;
  using free_matrix = matrix_type<E::free_size, E::free_size>;

  // Mapping from free track parameters.
  using free_to_bound_matrix = matrix_type<E::bound_size, E::free_size>;
  using free_to_path_matrix = matrix_type<1, E::free_size>;

  // Track helper
  using track_helper = detail::track_helper<matrix_actor, E>;

  /// @}

  ALGEBRA_HOST_DEVICE
  inline bound_vector free_to_bound_vector(const transform3_type& trf3,
                                           const free_vector& free_vec) const {
    const point3 pos = track_helper().pos(free_vec);
    const vector3 dir = track_helper().dir(free_vec);

    const point2 local =
        Derived<transform3_t, E>().global_to_local(trf3, pos, dir);

    bound_vector bound_vec;
    matrix_actor().element(bound_vec, E::bound_loc0, 0) = local[0];
    matrix_actor().element(bound_vec, E::bound_loc1, 0) = local[1];
    matrix_actor().element(bound_vec, E::bound_phi, 0) =
        vector_actor().phi(dir);
    matrix_actor().element(bound_vec, E::bound_theta, 0) =
        vector_actor().theta(dir);
    matrix_actor().element(bound_vec, E::bound_time, 0) =
        matrix_actor().element(free_vec, E::free_time, 0);
    matrix_actor().element(bound_vec, E::bound_qoverp, 0) =
        matrix_actor().element(free_vec, E::free_qoverp, 0);

    return bound_vec;
  }

  template <typename mask_t>
  ALGEBRA_HOST_DEVICE inline free_vector bound_to_free_vector(
      const transform3_type& trf3, const mask_t& mask,
      const bound_vector& bound_vec) const {

    const point2 local = track_helper().local(bound_vec);
    const vector3 dir = track_helper().dir(bound_vec);

    const auto pos =
        Derived<transform3_t, E>().local_to_global(trf3, mask, local, dir);

    free_vector free_vec;
    matrix_actor().element(free_vec, E::free_pos0, 0) = pos[0];
    matrix_actor().element(free_vec, E::free_pos1, 0) = pos[1];
    matrix_actor().element(free_vec, E::free_pos2, 0) = pos[2];
    matrix_actor().element(free_vec, E::free_time, 0) =
        matrix_actor().element(bound_vec, E::bound_time, 0);
    matrix_actor().element(free_vec, E::free_dir0, 0) = dir[0];
    matrix_actor().element(free_vec, E::free_dir1, 0) = dir[1];
    matrix_actor().element(free_vec, E::free_dir2, 0) = dir[2];
    matrix_actor().element(free_vec, E::free_qoverp, 0) =
        matrix_actor().element(bound_vec, E::bound_qoverp, 0);

    return free_vec;
  }

  ALGEBRA_HOST_DEVICE inline bound_to_free_matrix bound_to_free_jacobian(
      const transform3_type& trf3, const bound_vector& bound_vec) {

    // Declare jacobian for bound to free coordinate transform
    bound_to_free_matrix jac_to_global =
        matrix_actor().template zero<E::free_size, E::bound_size>();

    // Get trigonometric values
    const scalar_type theta =
        matrix_actor().element(bound_vec, E::bound_theta, 0);
    const scalar_type phi = matrix_actor().element(bound_vec, E::bound_phi, 0);
    const scalar_type cos_theta = std::cos(theta);
    const scalar_type sin_theta = std::sin(theta);
    const scalar_type cos_phi = std::cos(phi);
    const scalar_type sin_phi = std::sin(phi);

    // Get d(x_glo,y_glo,z_glo)/d(x_loc, y_loc)
    const matrix_type<3, 2> bound_to_free_rotation =
        Derived<transform3_t, E>().bound_to_free_rotation(trf3, bound_vec);

    matrix_actor().template set_block(jac_to_global, bound_to_free_rotation,
                                      E::free_pos0, E::bound_loc0);

    matrix_actor().element(jac_to_global, E::free_time, E::bound_time) = 1;
    matrix_actor().element(jac_to_global, E::free_dir0, E::bound_phi) =
        -1 * sin_theta * sin_phi;
    matrix_actor().element(jac_to_global, E::free_dir0, E::bound_theta) =
        cos_theta * cos_phi;
    matrix_actor().element(jac_to_global, E::free_dir1, E::bound_phi) =
        sin_theta * cos_phi;
    matrix_actor().element(jac_to_global, E::free_dir1, E::bound_theta) =
        cos_theta * sin_phi;
    matrix_actor().element(jac_to_global, E::free_dir2, E::bound_theta) =
        -1 * sin_theta;
    matrix_actor().element(jac_to_global, E::free_qoverp, E::bound_qoverp) = 1;

    return jac_to_global;
  }

  ALGEBRA_HOST_DEVICE
  inline free_to_bound_matrix free_to_bound_jacobian(
      const transform3_type& trf3, const free_vector& free_vec) {

    // Declare jacobian for bound to free coordinate transform
    free_to_bound_matrix jac_to_local =
        matrix_actor().template zero<E::bound_size, E::free_size>();

    // Free direction
    const vector3 dir = track_helper().dir(free_vec);

    // Get trigonometric values
    const scalar_type theta = vector_actor().theta(dir);
    const scalar_type phi = vector_actor().phi(dir);
    const scalar_type cos_theta = std::cos(theta);
    const scalar_type sin_theta = std::sin(theta);
    const scalar_type inv_sin_theta = 1. / sin_theta;
    const scalar_type cos_phi = std::cos(phi);
    const scalar_type sin_phi = std::sin(phi);

    // Get d(x_loc, y_loc)/d(x_glo,y_glo,z_glo)
    const matrix_type<2, 3> free_to_bound_rotation =
        Derived<transform3_t, E>().free_to_bound_rotation(trf3, free_vec);

    matrix_actor().template set_block(jac_to_local, free_to_bound_rotation,
                                      E::bound_loc0, E::free_pos0);

    // Set d(Free time)/d(Bound time)
    matrix_actor().element(jac_to_local, E::bound_time, E::free_time) = 1;

    // Set d(phi, theta)/d(free dir)
    matrix_actor().element(jac_to_local, E::bound_phi, E::free_dir0) =
        -1. * sin_phi * inv_sin_theta;
    matrix_actor().element(jac_to_local, E::bound_phi, E::free_dir1) =
        cos_phi * inv_sin_theta;
    matrix_actor().element(jac_to_local, E::bound_theta, E::free_dir0) =
        cos_phi * cos_theta;
    matrix_actor().element(jac_to_local, E::bound_theta, E::free_dir1) =
        sin_phi * cos_theta;
    matrix_actor().element(jac_to_local, E::bound_theta, E::free_dir2) =
        -1 * sin_theta;

    // Set d(Free Qop)/d(Bound Qop)
    matrix_actor().element(jac_to_local, E::bound_qoverp, E::free_qoverp) = 1;

    return jac_to_local;
  }
};

}  // namespace algebra::common