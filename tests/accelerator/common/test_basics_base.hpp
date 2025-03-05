/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "test_base.hpp"

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cmath>
#include <memory>

/// Test case class, to be specialised for the different plugins and languages
template <typename T>
class test_basics_base : public testing::Test, public test_base<T> {

 public:
  /// Constructor, setting up the inputs for all of the tests
  test_basics_base() = default;

 protected:
  virtual void fill_pointers() final {
    // Initialise the input and output vectors.
    for (std::size_t i = 0; i < s_arraySize; ++i) {

      m_t1[i] = {static_cast<typename T::scalar>(1.1),
                 static_cast<typename T::scalar>(2.2),
                 static_cast<typename T::scalar>(3.3)};
      m_t2[i] = {static_cast<typename T::scalar>(4.4),
                 static_cast<typename T::scalar>(5.5),
                 static_cast<typename T::scalar>(6.6)};
      m_t3[i] = {static_cast<typename T::scalar>(7.7),
                 static_cast<typename T::scalar>(8.8),
                 static_cast<typename T::scalar>(9.9)};

      m_p1[i] = {static_cast<typename T::scalar>(i * 0.5),
                 static_cast<typename T::scalar>((i + 1) * 1.0)};
      m_p2[i] = {static_cast<typename T::scalar>((i + 2) * 1.2),
                 static_cast<typename T::scalar>(i * 0.6)};

      m_v1[i] = {static_cast<typename T::scalar>(i * 0.6),
                 static_cast<typename T::scalar>((i + 1) * 1.2),
                 static_cast<typename T::scalar>((i + 2) * 1.3)};
      m_v2[i] = {static_cast<typename T::scalar>((i + 1) * 1.8),
                 static_cast<typename T::scalar>(i * 2.3),
                 static_cast<typename T::scalar>((i + 2) * 3.4)};

      for (typename T::size_type j = 0; j < 6; ++j) {
        for (typename T::size_type k = 0; k < 4; ++k) {
          algebra::getter::element(m_m1[i], j, k) =
              static_cast<typename T::scalar>(j * 20.3 + k * 10.5);
        }
      }

      algebra::getter::element(m_m2[i], 0, 0) = 4;
      algebra::getter::element(m_m2[i], 0, 1) = 3;
      algebra::getter::element(m_m2[i], 1, 0) = 12;
      algebra::getter::element(m_m2[i], 1, 1) = 13;

      m_output_host[i] = 0;
      m_output_device[i] = 0;
    }
  }

  /// Compare the outputs, after the data processing is finished.
  void compareOutputs() const {
    for (std::size_t i = 0; i < this->s_arraySize; ++i) {
      EXPECT_NEAR(m_output_host[i], m_output_device[i],
                  std::abs(0.001 * m_output_host[i]));
    }
  }

  /// Size for the tested arrays.
  static constexpr std::size_t s_arraySize = 5000;

  /// @name Inputs for the tests
  /// @{

  typename T::vector3 *m_t1, *m_t2, *m_t3;
  typename T::point2 *m_p1, *m_p2;
  typename T::vector3 *m_v1, *m_v2;
  typename T::template matrix<6, 4> *m_m1;
  typename T::template matrix<2, 2> *m_m2;

  /// @}

  /// @name Outputs for the tests
  /// @{

  typename T::scalar *m_output_host, *m_output_device;

  /// @}

};  // class test_basics_base
