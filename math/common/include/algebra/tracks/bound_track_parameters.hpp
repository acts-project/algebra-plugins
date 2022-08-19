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

template <typename matrix_actor_t, typename vector_actor_t,
          typename track_indices_t>
struct bound_track_parameters {

  /// @name Type definitions for the struct
  /// @{

  /// Matrix actor
  using matrix_actor = matrix_actor_t;
  /// Vector actor
  using vector_actor = vector_actor_t;
  /// Size type
  using size_type = typename matrix_actor_t::size_ty;
  /// Scalar type
  using scalar_type = typename matrix_actor_t::scalar_type;
  /// 2D Matrix type
  template <size_type ROWS, size_type COLS>
  using matrix_type = typename matrix_actor::template matrix_type<ROWS, COLS>;
  /// Array type
  template <size_type N>
  using array_type = typename matrix_actor::template array_type<N>;
  /// 3-element "vector" type
  using vector3 = array_type<3>;
  /// Point in 3D space
  using point3 = vector3;
  /// Point in 2D space
  using point2 = array_type<2>;

  /// Track indices
  using E = track_indices_t;

  // Shorthand vector/matrix types related to bound track parameters.
  using vector_type = matrix_type<E::bound_size, 1>;
  using covariance_type = matrix_type<E::bound_size, E::bound_size>;

  // Track helper
  using track_helper = detail::track_helper<matrix_actor, vector_actor, E>;

  /// @}

  ALGEBRA_HOST_DEVICE
  bound_track_parameters()
      : m_vector(matrix_actor().template zero<E::bound_size, 1>()),
        m_covariance(
            matrix_actor().template zero<E::bound_size, E::bound_size>()) {}

  ALGEBRA_HOST_DEVICE
  bound_track_parameters(const std::size_t sf_idx, const vector_type& vec,
                         const covariance_type& cov)
      : m_surface_link(sf_idx), m_vector(vec), m_covariance(cov) {}

  ALGEBRA_HOST_DEVICE
  const std::size_t& surface_link() const { return m_surface_link; }

  ALGEBRA_HOST_DEVICE
  const vector_type& vector() const { return m_vector; }

  ALGEBRA_HOST_DEVICE
  void set_vector(const vector_type& v) { m_vector = v; }

  ALGEBRA_HOST_DEVICE
  const covariance_type& covariance() const { return m_covariance; }

  ALGEBRA_HOST_DEVICE
  void set_covariance(const covariance_type& c) { m_covariance = c; }

  ALGEBRA_HOST_DEVICE
  point3 local() const { return track_helper().local(m_vector); }

  ALGEBRA_HOST_DEVICE
  scalar_type phi() const {
    return matrix_actor().element(m_vector, E::bound_phi, 0);
  }

  ALGEBRA_HOST_DEVICE
  scalar_type theta() const {
    return matrix_actor().element(m_vector, E::bound_theta, 0);
  }

  ALGEBRA_HOST_DEVICE
  vector3 dir() const { return track_helper().dir(m_vector); }

  ALGEBRA_HOST_DEVICE
  scalar_type time() const {
    return matrix_actor().element(m_vector, E::bound_time, 0);
  }

  ALGEBRA_HOST_DEVICE
  scalar_type charge() const {
    if (matrix_actor().element(m_vector, E::bound_qoverp, 0) < 0) {
      return -1.;
    } else {
      return 1.;
    }
  }

  ALGEBRA_HOST_DEVICE
  scalar_type qop() const {
    return matrix_actor().element(m_vector, E::bound_qoverp, 0);
  }

 private:
  std::size_t m_surface_link;
  vector_type m_vector;
  covariance_type m_covariance;
};

}  // namespace algebra::common