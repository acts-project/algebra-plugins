 /**
 * author: joana.niermann@cern.ch
 **/
#include <limits>
#include <type_traits>

//#include <Eigen/Dense>
//#include <Eigen/Geometry>

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
    using vector = Vc::Common::AdaptSubscriptOperator<std::vector<T, Allocator> >;

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

/*using Vector3_s    = Eigen::matrix<scalar, 3, 1>;
using Vector4_s    = Eigen::matrix<scalar, 4, 1>;
using VectorV_s    = Eigen::matrix<scalar, scalar_v::Size, 1>;
using transform4_s = Eigen::transform<scalar, 4, Eigen::Affine>;*/
template<typename scalar, unsigned int kDIM>
using array        = Vc::SimdArray<scalar, kDIM>;
using vector_s     = aligned::vector<scalar>;
using vectorV_s    = std::array<scalar, scalar_v::Size>;
using vector3_v    = Vector3<scalar_v>;
using transform4_s = std::array<scalar, 16>;

// Needed in AoS
template<typename data_t>
struct Vector2
{
  data_t x, y;
};

template<typename data_t>
struct Vector3
{
  data_t x, y, z;
};

template<typename data_t>
struct Vector4
{
  data_t x, y, z, t;
};

template<typename data_t, unsigned int kDIM>
struct VectorV
{
  data_t x[kDIM];
};

template<typename data_t>
bool operator==(const Vector4<data_t> &lhs, const Vector4<data_t> &rhs) {
    bool is_eq_x = (lhs.x == rhs.x).isFull();
    bool is_eq_y = (lhs.y == rhs.y).isFull();
    bool is_eq_z = (lhs.z == rhs.z).isFull();
    bool is_eq_t = (lhs.t == rhs.t).isFull();

    return is_eq_x && is_eq_y && is_eq_z && is_eq_t;
  }

/*template<typename data_t>
using point2 = Vector2<data_t>;

template<typename data_t>
using point3 = Vector2<data_t>;*/

// Assume 4x4 transformation matrix (placement) as input
template<typename data_t>
struct transform4_v
{
  // the rest of the 4x4 matrix (avoid padding)
  Vector4<data_t> vec0, vec1;
  // plane normal
  Vector4<data_t> normal;
  // plane translation
  Vector4<data_t> translation;
};

// Assume 4x4 transformation (placement) as input 
// Eigen compact format (3x3 in memory)
template<typename data_t>
struct transform3_v
{
  // the rest of the 3x3 matrix (avoid padding)
  Vector3<data_t> vec0, vec1;
  // plane normal
  Vector3<data_t> normal;
  // plane translation
  Vector3<data_t> translation;
};

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
};*/

// detray style intersection structure to be filled as result type
//template <typename scalar_t, typename vector_t = vector3_s>
//struct /*alignas(alignment)*/ intersection {
//  vector_t path;
//  scalar_t dist = std::numeric_limits<scalar_t>::infinity();
//};

// define type that holds input data: has to be memory layout compatible with the Vc Vector type, i.e. 
//                                    has to have same ABI as array of T with rows*cols many entries
// To be specialized by the LA plugins
/*template<typename T, class Enable = void>
struct data_trait {
  using type = void;
  using value_type = void;

  static constexpr bool is_vec_layout = false;
  static constexpr bool is_std_layout = false;
  static constexpr bool is_same_type  = false;
  static constexpr bool is_vec_dim    = false;
};*/

//---------------------------------------------------Plugin specific

// Small-time wrapper for a uniform interface, will be taken care of by the LA plugin
// later
/*template<typename Derived>
struct alignas(alignment) matrixV {
  using scalar_type = typename Eigen::DenseBase<Derived>::scalar;
  using vec_type = typename Eigen::DenseBase<Derived>::scalar;
  using type   = Derived;

  type obj;

  type& operator()() {return obj;}

  //void init() {obj = decltype(obj)::Random();}
  constexpr size_t n_elemts() {return obj.rows() * obj.cols();}
  const scalar_type* data() {return obj.data();}
  constexpr size_t padding() {return (alignment - (n_elemts() * sizeof(scalar_type)) % alignment) % alignment / sizeof(scalar_type);}
};

// Affine transform is not derived from Eigen::DenseBase
template<>
struct alignas(alignment) matrixV<transform4_s> {
  using scalar_type = typename transform4_s::scalar;
  using vec_type = typename transform4_s::scalar;
  using type   = transform4_s;

  transform4_s obj;

  type& operator()() {return obj;}
   
  const size_t n_elemts() {return obj.rows() * obj.cols();}
  const scalar_type* data() {return obj.data();}
  const size_t padding() {return (alignment - (n_elemts() * sizeof(scalar_type)) % alignment) % alignment / sizeof(scalar_type);}
};

// 3D Structure of Vector data
template<>
struct alignas(alignment) matrixV<Vector3<scalar_v> > {
  using scalar_type = scalar_v::value_type;
  using vec_type = scalar_v;
  using type   = Vector3<scalar_v>;

  Vector3<scalar_v> obj;

  //matrixV() {}
  type& operator()() {return obj;}

  const size_t n_elemts() {return 3*scalar_v::Size;}
  const scalar_type* data() {return reinterpret_cast<const scalar_type*>(&obj.x);}
  const size_t padding() {return (alignment - (n_elemts() * sizeof(scalar_type)) % alignment) % alignment / sizeof(scalar_type);}
};

// Eigen specific types
template<typename Derived>
struct data_trait {
  using type  = matrixV<Derived>;
  using scalar_type = typename type::scalar_type;
  using vec_type = typename type::vec_type;
  static const size_t size = sizeof(type)/sizeof(scalar_type);

  static constexpr bool is_std_layout = std::is_standard_layout_v<type>;
  static constexpr bool is_same_type  = std::is_same_v<scalar_type, scalar_v::value_type>;
  static constexpr bool is_vec_dim    = !(scalar_v::Size > size ? scalar_v::Size % size : size % scalar_v::Size);
  static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
};

// Affine transform is not derived from Eigen::DenseBase
template<>
struct data_trait<transform4_s> {
  using type  = matrixV<transform4_s>;
  using scalar_type = transform4_s::scalar;
  using vec_type = typename type::vec_type;
  static const size_t size = sizeof(type)/sizeof(scalar_type);

  static constexpr bool is_std_layout = std::is_standard_layout_v<type> 
                                        && std::is_standard_layout_v<transform4_s::VectorType>;
  static constexpr bool is_same_type  = std::is_same_v<scalar_type, scalar_v::value_type>;
  static constexpr bool is_vec_dim    = scalar_v::Size > size ? scalar_v::Size % size == 0 : size % scalar_v::Size == 0;
  static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
};

template<>
struct data_trait<Vector3<scalar_v> > {
  using type  = matrixV<Vector3<scalar_v> >;
  using scalar_type = scalar_v::value_type;
  using vec_type = typename type::vec_type;
  static const size_t size = sizeof(type)/sizeof(scalar_type);

  static constexpr bool is_std_layout = std::is_standard_layout_v<type>;
  static constexpr bool is_same_type  = std::is_same_v<scalar_type, scalar_v::value_type>;
  static constexpr bool is_vec_dim    = scalar_v::Size > size ? scalar_v::Size % size == 0 : size % scalar_v::Size == 0;
  static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
};*/
  } //namespace simd
} //namespace algebra