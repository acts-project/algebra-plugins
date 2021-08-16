/** Algebra plugins, part of the ACTS project
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <Vc/Vc>


#ifdef ALGEBRA_CUSTOM_SCALARTYPE
using algebra_scalar = ALGEBRA_CUSTOM_SCALARTYPE;
using algebra_scalar_v = Vc::Vector<ALGEBRA_CUSTOM_SCALARTYPE>;
#else
using algebra_scalar = double;
using algebra_scalar_v = Vc::double_v;
#endif

namespace algebra {

using scalar   = algebra_scalar;

// Namespace for simd types
namespace simd {
  using scalar_v = algebra_scalar_v;

  // Default alignment
  constexpr size_t alignment = alignof(scalar_v);

  // Namespace for aligned containers
  namespace aligned {
    // Allow Vc to get alignment right
    //template <typename T, typename Allocator = std::allocator<T> >
    //template <typename T, typename Allocator = Eigen::aligned_allocator<T> >
    template <typename T, typename Allocator = Vc::Allocator<T> >
    // Add subscript operator to allow for gather operations from AoS
    //using vector = Vc::Common::AdaptSubscriptOperator<std::vector<T, Allocator> >;
    using vector = std::vector<T, Allocator>;

    template<size_t kDIM>
    using mem_t = Vc::Memory<scalar_v, kDIM>;
  } //namespace aligned

  template<typename scalar, unsigned int kDIM>
  using array = Vc::SimdArray<scalar, kDIM>;

  /** 2D Structured vector data. This structure is designed to also hold data for 
   *  horizontal vectorizarion or matrices
   * 
   * @tparam The type of data for the respective coordinate
   * 
   * @note The structure is aligned to vector types, so padding has to be considerd
  **/
  template<typename data_t>
  struct alignas(data_t) Vector2
  {
    data_t x, y;
  };

  /** 3D Structured vector data. This structure is designed to also hold data for 
   *  horizontal vectorizarion or matrices
   * 
   * @tparam The type of data for the respective coordinate
   * 
   * @note 
  **/
  template<typename data_t>
  struct Vector3
  {
    data_t x, y, z;
  };

  /** 4D Structured vector data. This structure is designed to also hold data for 
   *  horizontal vectorizarion or tranformation matrices
   * 
   * @tparam The type of data for the respective coordinate
   * 
   * @note The structures define custom alignment, so watch for padding
  **/
  template<typename data_t>
  struct Vector4
  {
    data_t x, y, z, t;
  };

  /** Structured vector data
   * 
   * @tparam The type of data for the respective coordinate
   * 
   * @note This structure is designed to also hold data for 
   *       horizontal vectorizarion or tranformation matrices
  **/
  template<typename data_t, unsigned int kDIM>
  struct alignas(data_t) VectorV
  {
    data_t x[kDIM];
  };

  /** Equality operator definition for Vc::SimdArray in a simd Vector4 struct.
   *
   * @tparam scalar_t value type used in the wrapper
   *
   * @return Result of comparision
   */
  template<typename scalar_t, size_t kDIM>
  bool operator==(const Vector4<simd::array<scalar_t, kDIM>> &lhs, 
                  const Vector4<simd::array<scalar_t, kDIM>> &rhs) {
    bool is_eq_x = (lhs.x == rhs.x).isFull();
    bool is_eq_y = (lhs.y == rhs.y).isFull();
    bool is_eq_z = (lhs.z == rhs.z).isFull();
    bool is_eq_t = (lhs.t == rhs.t).isFull();

    return is_eq_x && is_eq_y && is_eq_z && is_eq_t;
  }

} //namespace simd

} //namespace algebra