/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/qualifiers.hpp"

// Test include(s).
#include "algebra/test/framework/types.hpp"
#include "device_fixture.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cmath>
#include <memory>

namespace algebra::test {

/// Data fixture for device test cases
template <algebra::concepts::algebra A>
class vector_fixture : public device_fixture<A, algebra::get_scalar_t<A>> {

  using scalar_t = algebra::get_scalar_t<A>;
  using point2_t = algebra::get_point2D_t<A>;
  using vector3_t = algebra::get_vector3D_t<A>;

  using result_t = scalar_t;
  using base_fixture = device_fixture<A, result_t>;

 public:
  /// Constructor, setting up the inputs for all of the tests
  vector_fixture(vecmem::memory_resource& mr) : base_fixture(mr) {}

 protected:
  /// Function setting things up for a test
  virtual void SetUp() override {

    // Set up the result vectors
    base_fixture::SetUp();

    // Set up the input vectors.
    m_p1 = std::make_unique<vecmem::vector<point2_t>>(this->size(),
                                                      &this->resource());
    m_p2 = std::make_unique<vecmem::vector<point2_t>>(this->size(),
                                                      &this->resource());

    m_v1 = std::make_unique<vecmem::vector<vector3_t>>(this->size(),
                                                       &this->resource());
    m_v2 = std::make_unique<vecmem::vector<vector3_t>>(this->size(),
                                                       &this->resource());

    // Initialise the input and output vectors.
    for (std::size_t i = 0; i < this->size(); ++i) {

      m_p1->at(i) = {static_cast<scalar_t>(i * 0.5),
                     static_cast<scalar_t>((i + 1) * 1.0)};
      m_p2->at(i) = {static_cast<scalar_t>((i + 2) * 1.2),
                     static_cast<scalar_t>(i * 0.6)};

      m_v1->at(i) = {static_cast<scalar_t>(i * 0.6),
                     static_cast<scalar_t>((i + 1) * 1.2),
                     static_cast<scalar_t>((i + 2) * 1.3)};
      m_v2->at(i) = {static_cast<scalar_t>((i + 1) * 1.8),
                     static_cast<scalar_t>(i * 2.3),
                     static_cast<scalar_t>((i + 2) * 3.4)};
    }
  }

  /// Function tearing things down after the test
  virtual void TearDown() override {
    // Delete the vectors.
    m_p1.reset();
    m_p2.reset();
    m_v1.reset();
    m_v2.reset();
  }

  /// @name Inputs for the tests
  /// @{

  std::unique_ptr<vecmem::vector<point2_t>> m_p1, m_p2;
  std::unique_ptr<vecmem::vector<vector3_t>> m_v1, m_v2;

  /// @}

};  // class device_fixture

/// Functor running @c test_device_basics::vector_2d_ops
template <algebra::concepts::algebra A>
class vector_2d_ops_functor {

  using scalar_t = algebra::get_scalar_t<A>;
  using point2_t = algebra::get_point2D_t<A>;

 public:
  ALGEBRA_HOST_DEVICE void operator()(
      std::size_t i, vecmem::data::vector_view<const point2_t> a,
      vecmem::data::vector_view<const point2_t> b,
      vecmem::data::vector_view<scalar_t> output) const {

    // Create the VecMem vector(s).
    vecmem::device_vector<const point2_t> vec_a(a), vec_b(b);
    vecmem::device_vector<scalar_t> vec_output(output);

    // Perform the operation.
    auto ii = static_cast<typename decltype(vec_output)::size_type>(i);
    vec_output[ii] = vector_2d_ops(vec_a[ii], vec_b[ii]);
  }

 private:
  /// Perform various 2D vector operations, and produce a scalar output
  ALGEBRA_HOST_DEVICE
  scalar_t vector_2d_ops(point2_t a, point2_t b) const {
    using namespace algebra;

    point2_t c = a + b;
    point2_t c2 = c * 2.0;

    scalar_t phi = algebra::vector::phi(c2);
    scalar_t perp = algebra::vector::perp(c2);
    scalar_t norm1 = algebra::vector::norm(c2);

    scalar_t dot = algebra::vector::dot(a, b);
    point2_t norm2 = algebra::vector::normalize(c);
    scalar_t norm3 = algebra::vector::norm(norm2);

    return (phi + perp + norm1 + dot + norm3);
  }
};

/// Functor running @c test_device_basics::vector_3d_ops
template <algebra::concepts::algebra A>
class vector_3d_ops_functor {

  using scalar_t = algebra::get_scalar_t<A>;
  using vector3_t = algebra::get_vector3D_t<A>;

 public:
  ALGEBRA_HOST_DEVICE void operator()(
      std::size_t i, vecmem::data::vector_view<const vector3_t> a,
      vecmem::data::vector_view<const vector3_t> b,
      vecmem::data::vector_view<scalar_t> output) const {

    // Create the VecMem vector(s).
    vecmem::device_vector<const vector3_t> vec_a(a), vec_b(b);
    vecmem::device_vector<scalar_t> vec_output(output);

    // Perform the operation.
    auto ii = static_cast<typename decltype(vec_output)::size_type>(i);
    vec_output[ii] = vector_3d_ops(vec_a[ii], vec_b[ii]);
  }

 private:
  /// Perform various 3D vector operations, and produce a scalar output
  ALGEBRA_HOST_DEVICE
  scalar_t vector_3d_ops(vector3_t a, vector3_t b) const {
    using namespace algebra;

    vector3_t c = a + b;
    vector3_t c2 = c * 2.0;

    scalar_t phi = algebra::vector::phi(c2);
    scalar_t perp = algebra::vector::perp(c2);
    scalar_t norm1 = algebra::vector::norm(c2);

    vector3_t d = algebra::vector::cross(a, b);

    scalar_t dot = algebra::vector::dot(a, d);
    vector3_t norm2 = algebra::vector::normalize(c);
    scalar_t norm3 = algebra::vector::norm(norm2);

    return (phi + perp + norm1 + dot + norm3);
  }
};

}  // namespace algebra::test
