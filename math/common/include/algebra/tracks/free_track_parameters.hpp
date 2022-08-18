/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

namespace algebra::common {

template <typename matrix_actor_t, typename vector_actor_t,
          typename track_indices_t>
struct free_track_parameters {

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

  // Shorthand vector/matrix types related to free track parameters.
  using vector_type = matrix_type<E::free_size, 1>;
  using covariance_type = matrix_type<E::free_size, E::free_size>;

  // Track helper
  using track_helper = detail::track_helper<matrix_actor, E>;

  /// @}

  ALGEBRA_HOST_DEVICE
  free_track_parameters()
      : m_vector(matrix_actor().template zero<E::free_size, 1>()),
        m_covariance(
            matrix_actor().template zero<E::free_size, E::free_size>()){};

  ALGEBRA_HOST_DEVICE
  free_track_parameters(const vector_type& vec, const covariance_type& cov)
      : m_vector(vec), m_covariance(cov) {}

  ALGEBRA_HOST_DEVICE
  free_track_parameters(const point3& pos, const scalar_type time,
                        const vector3& mom, const scalar_type q) {

    matrix_actor().element(m_vector, E::free_pos0, 0) = pos[0];
    matrix_actor().element(m_vector, E::free_pos1, 0) = pos[1];
    matrix_actor().element(m_vector, E::free_pos2, 0) = pos[2];
    matrix_actor().element(m_vector, E::free_time, 0) = time;

    scalar_type p = vector_actor().norm(mom);
    auto mom_norm = vector_actor().normalize(mom);
    matrix_actor().element(m_vector, E::free_dir0, 0) = mom_norm[0];
    matrix_actor().element(m_vector, E::free_dir1, 0) = mom_norm[1];
    matrix_actor().element(m_vector, E::free_dir2, 0) = mom_norm[2];
    matrix_actor().element(m_vector, E::free_qoverp, 0) = q / p;
  }

  ALGEBRA_HOST_DEVICE
  const vector_type& vector() const { return m_vector; }

  ALGEBRA_HOST_DEVICE
  void set_vector(const vector_type& v) { m_vector = v; }

  ALGEBRA_HOST_DEVICE
  const covariance_type& covariance() const { return m_covariance; }

  ALGEBRA_HOST_DEVICE
  void set_covariance(const covariance_type& c) { m_covariance = c; }

  ALGEBRA_HOST_DEVICE
  scalar_type overstep_tolerance() const { return m_overstep_tolerance; }

  ALGEBRA_HOST_DEVICE
  void set_overstep_tolerance(const scalar_type tolerance) {
    m_overstep_tolerance = tolerance;
  }

  ALGEBRA_HOST_DEVICE
  point3 pos() const { return track_helper().pos(m_vector); }

  ALGEBRA_HOST_DEVICE
  void set_pos(const vector3& pos) { track_helper().set_pos(m_vector, pos); }

  ALGEBRA_HOST_DEVICE
  vector3 dir() const { return track_helper().dir(m_vector); }

  ALGEBRA_HOST_DEVICE
  void set_dir(const vector3& dir) { track_helper().set_dir(m_vector, dir); }

  ALGEBRA_HOST_DEVICE
  vector3 mom() const { return 1. / std::abs(this->qop()) * this->dir(); }

  ALGEBRA_HOST_DEVICE
  scalar_type time() const {
    return matrix_actor().element(m_vector, E::free_time, 0);
  }

  ALGEBRA_HOST_DEVICE
  scalar_type charge() const {
    if (matrix_actor().element(m_vector, E::free_qoverp, 0) < 0) {
      return -1.;
    } else {
      return 1.;
    }
  }

  ALGEBRA_HOST_DEVICE
  scalar_type qop() const {
    return matrix_actor().element(m_vector, E::free_qoverp, 0);
  }

  ALGEBRA_HOST_DEVICE
  scalar_type pT() const {
    auto dir = this->dir();
    return std::abs(1. / this->qop() * vector_actor().perp(dir));
  }

  ALGEBRA_HOST_DEVICE
  void flip() {
    matrix_actor().element(m_vector, E::free_dir0, 0) *= -1.;
    matrix_actor().element(m_vector, E::free_dir1, 0) *= -1.;
    matrix_actor().element(m_vector, E::free_dir2, 0) *= -1.;
    matrix_actor().element(m_vector, E::free_qoverp, 0) *= -1.;
  }

 private:
  vector_type m_vector;
  covariance_type m_covariance;
  scalar_type m_overstep_tolerance = -1e-4;
};

}  // namespace algebra::common