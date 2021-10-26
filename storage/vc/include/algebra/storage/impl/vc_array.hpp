/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Vc include(s).
#include <Vc/Vc>

// System include(s).
#include <cassert>
#include <cstddef>

namespace algebra::vc {

/// Vc array type
template <typename T, std::size_t N>
class array : public Vc::SimdArray<T, N> {

 public:
  /// Inherit all constructors from the base class
  using Vc::SimdArray<T, N>::SimdArray;

};  // class array

/// Special Vc array type for size-3 arrays
template <typename T>
class array<T, 3> : public Vc::SimdArray<T, 4> {

 public:
  /// Inherit all constructors from the base class
  using Vc::SimdArray<T, 4>::SimdArray;

  /// Define a constructor that could receive just 3 elements
  array(const std::initializer_list<T>& init)
      : Vc::SimdArray<T, 4>(
            {*(init.begin()), *(init.begin() + 1), *(init.begin() + 2), 0}) {
    assert(init.size() == 3);
  }

};  // class array

}  // namespace algebra::vc
