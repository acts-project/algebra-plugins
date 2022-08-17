/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

namespace algebra {

/** Coordinate base struct
 */
template <typename transform3_t>
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

  /// Components of a bound track parameters vector.
  ///
  enum bound_indices : size_type {
    // Local position on the reference surface.
    // This is intentionally named different from the position components in
    // the other data vectors, to clarify that this is defined on a surface
    // while the others are defined in free space.
    e_bound_loc0 = 0,
    e_bound_loc1 = 1,
    // Direction angles
    e_bound_phi = 2,
    e_bound_theta = 3,
    // Global inverse-momentum-like parameter, i.e. q/p or 1/p
    // The naming is inconsistent for the case of neutral track parameters where
    // the value is interpreted as 1/p not as q/p. This is intentional to avoid
    // having multiple aliases for the same element and for lack of an
    // acceptable
    // common name.
    e_bound_qoverp = 4,
    e_bound_time = 5,
    // Last uninitialized value contains the total number of components
    e_bound_size,
  };

  /// Components of a free track parameters vector.
  ///
  /// To be used to access components by named indices instead of just numbers.
  /// This must be a regular `enum` and not a scoped `enum class` to allow
  /// implicit conversion to an integer. The enum value are thus visible
  /// directly in `namespace Acts` and are prefixed to avoid naming collisions.
  enum free_indices : size_type {
    // Spatial position
    // The spatial position components must be stored as one continous block.
    e_free_pos0 = 0u,
    e_free_pos1 = e_free_pos0 + 1u,
    e_free_pos2 = e_free_pos0 + 2u,
    // Time
    e_free_time = 3u,
    // (Unit) direction
    // The direction components must be stored as one continous block.
    e_free_dir0 = 4u,
    e_free_dir1 = e_free_dir0 + 1u,
    e_free_dir2 = e_free_dir0 + 2u,
    // Global inverse-momentum-like parameter, i.e. q/p or 1/p
    // See BoundIndices for further information
    e_free_qoverp = 7u,
    // Last uninitialized value contains the total number of components
    e_free_size,
  };

  // Shorthand vector/matrix types related to bound track parameters.
  using bound_vector = matrix_type<e_bound_size, 1>;
  using bound_matrix = matrix_type<e_bound_size, e_bound_size>;

  // Mapping from bound track parameters.
  using bound_to_free_matrix = matrix_type<e_free_size, e_bound_size>;

  // Shorthand vector/matrix types related to free track parameters.
  using free_vector = matrix_type<e_free_size, 1>;
  using free_matrix = matrix_type<e_free_size, e_free_size>;

  // Mapping from free track parameters.
  using free_to_bound_matrix = matrix_type<e_bound_size, e_free_size>;
  using free_to_path_matrix = matrix_type<1, e_free_size>;

  /*
  ALGEBRA_HOST_DEVICE
  inline bound_to_free_matrix bound_to_free_jacobian(
      const transform3_type& trf3, const bound_vector& bound_vec) {

    // Declare jacobian for bound to free coordinate transform
    bound_to_free_matrix jac_to_global =
        matrix_actor().template zero<e_free_size, e_bound_size>();

    // Get trigonometric values
    const scalar_type theta =
        matrix_actor().element(bound_vec, e_bound_theta, 0);
    const scalar_type phi = matrix_actor().element(bound_vec, e_bound_phi, 0);
    const scalar_type cos_theta = std::cos(theta);
    const scalar_type sin_theta = std::sin(theta);
    const scalar_type cos_phi = std::cos(phi);
    const scalar_type sin_phi = std::sin(phi);

    // Set d(x,y,z)/d(loc0,loc1)
    const matrix_type<3, 2> bound_to_free_rotation =
        matrix_actor().template block<3, 2>(trf3.matrix(), 0, 0);

    matrix_actor().template set_block(jac_to_global, bound_to_free_rotation,
                                      e_free_pos0, e_bound_loc0);

    matrix_actor().element(jac_to_global, e_free_time, e_bound_time) = 1;
    matrix_actor().element(jac_to_global, e_free_dir0, e_bound_phi) =
        -1 * sin_theta * sin_phi;
    matrix_actor().element(jac_to_global, e_free_dir0, e_bound_theta) =
        cos_theta * cos_phi;
    matrix_actor().element(jac_to_global, e_free_dir1, e_bound_phi) =
        sin_theta * cos_phi;
    matrix_actor().element(jac_to_global, e_free_dir1, e_bound_theta) =
        cos_theta * sin_phi;
    matrix_actor().element(jac_to_global, e_free_dir2, e_bound_theta) =
        -1 * sin_theta;
    matrix_actor().element(jac_to_global, e_free_qoverp, e_bound_qoverp) = 1;

    return jac_to_global;
  }
  */
  /// @}
};

}  // namespace algebra