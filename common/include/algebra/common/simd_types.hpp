/** Algebra plugins, part of the ACTS project
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
 
#include <limits>
#include <type_traits>

#include <Vc/Vc>

//#include <vecmem/containers/Vector.hpp>
//#include <vecmem/memory/host_memory_resource.hpp>

#pragma once


#ifdef ALGEBRA_CUSTOM_SCALARTYPE
using algebra_scalar = ALGEBRA_CUSTOM_SCALARTYPE;
using algebra_scalar_v = Vc::Vector<ALGEBRA_CUSTOM_SCALARTYPE>;
#else
using algebra_scalar = double;
using algebra_scalar_v = Vc::double_v;
#endif

namespace algebra {

using scalar   = algebra_scalar;

namespace simd {
  using scalar_v = algebra_scalar_v;

  // Check for max alignment!!e.g. Vc::VectorAlignment 
  //constexpr size_t alignment = std::hardware_constructive_interference_size;
  constexpr size_t alignment = alignof(scalar_v);
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

    // TODO use placement new rather
    template<typename vector_s, typename mem_v>
    using storage_t = std::aligned_union_t<1, vector_s, mem_v>;

    template<typename vector_t, typename memory_t>
    union data_container {vector_t Vector;
                          memory_t memory;};
  } //namespace aligned


  template<typename data_t>
  struct Vector3;

  template<typename scalar, unsigned int kDIM>
  using array        = Vc::SimdArray<scalar, kDIM>;
  using vector_s     = aligned::vector<scalar>;
  using vector3_v    = Vector3<scalar_v>;

  // Needed in AoS
  template<typename data_t>
  struct alignas(alignment) Vector2
  {
    data_t x, y;
  };

  template<typename data_t>
  struct alignas(alignment) Vector3
  {
    data_t x, y, z;
  };

  template<typename data_t>
  struct alignas(alignment) Vector4
  {
    data_t x, y, z, t;
  };

  template<typename data_t, unsigned int kDIM>
  struct alignas(alignment) VectorV
  {
    data_t x[kDIM];
  };


  // Need a special ==operator for Vc types
  template<typename scalar, size_t kDIM>
  struct alignas(alignment) Vector4 <simd::array<scalar, kDIM>> {
    simd::array<scalar, kDIM> x, y, z, t;
  };

  template<typename scalar, size_t kDIM>
  bool operator==(const Vector4<simd::array<scalar, kDIM>> &lhs, 
                  const Vector4<simd::array<scalar, kDIM>> &rhs) {
    bool is_eq_x = (lhs.x == rhs.x).isFull();
    bool is_eq_y = (lhs.y == rhs.y).isFull();
    bool is_eq_z = (lhs.z == rhs.z).isFull();
    bool is_eq_t = (lhs.t == rhs.t).isFull();

    return is_eq_x && is_eq_y && is_eq_z && is_eq_t;
  }

  // Convenience types
  /*template<typename data_t>
  struct alignas(alignment) ray_data
  {
    data_t direction, point;
  };

  template<typename data_t>
  struct alignas(alignment) plane_data
  {
    data_t normals, points;
  };

  // Define vector types
  template<typename data_t>
  struct vector_trait {
    using type  = data_t;
    using scalar_type = typename scalar;
    using vec_type = typename scalar_v;
    static const size_t size = sizeof(type)/sizeof(scalar_type);

    static constexpr bool is_std_layout = std::is_standard_layout_v<type>;
    static constexpr bool is_same_type  = std::is_same_v<scalar_type, scalar_v::value_type>;
    static constexpr bool is_vec_dim    = !(scalar_v::Size > size ? scalar_v::Size % size : size % scalar_v::Size);
    static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
  };

  template<>
  struct vector_trait<Vector3<scalar_v> > {
    using type  = Vector3<scalar_v>;
    using scalar_type = scalar_v::value_type;
    using vec_type = typename scalar_v;
    static const size_t size = sizeof(type)/sizeof(scalar_type);

    static constexpr bool is_std_layout = std::is_standard_layout_v<type>;
    static constexpr bool is_same_type  = std::is_same_v<scalar_type, scalar_v::value_type>;
    static constexpr bool is_vec_dim    = scalar_v::Size > size ? scalar_v::Size % size == 0 : size % scalar_v::Size == 0;
    static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
  };*/

} //namespace simd
} //namespace algebra