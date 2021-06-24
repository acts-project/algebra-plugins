  #include "common/simd_types.hpp"
  
  namespace algebra {

  namespace simd {
  
  template<typename scalar_t>
  class array4_wrapper /*: public simd::array<scalar_t, 4>*/ {
    public:

    simd::array<scalar_t, 4> _array;

    using mask_type = typename decltype(_array)::mask_type;
    
    array4_wrapper() : _array() {}

    array4_wrapper(const array4_wrapper& array4) : _array(array4._array) {}

    array4_wrapper(array4_wrapper&& array4) : _array(std::move(array4._array)) {}

    array4_wrapper(scalar_t value) : _array(value) {}

    array4_wrapper(scalar_t v1, scalar_t v2, scalar_t v3, scalar_t v4 = scalar_t{0.0}) 
      : _array{v1, v2, v3, v4} {}


    array4_wrapper(simd::array<scalar_t, 4> base) : _array(std::move(base)) {}

    template<typename other_type>
    inline const array4_wrapper<scalar_t> & operator=(const other_type& lhs)
    {
        _array = lhs;
        return *this;
    }

    inline array4_wrapper<scalar_t>& operator=(const array4_wrapper<scalar_t>& other)
    {
        _array = other._array;
        return *this;
    }

    inline array4_wrapper<scalar_t>& operator=(std::initializer_list<scalar_t> &list)
    {
        _array = simd::array<scalar_t, 4>(list);
        return *this;
    }

    operator simd::array<scalar_t, 4>&() { return _array; }
    operator const simd::array<scalar_t, 4>&() const { return _array; }

    inline auto operator[](int i) {return _array[i];}
    inline auto operator[](int i) const {return _array[i];}

    inline auto operator*=(scalar_t factor) ->decltype(_array*=factor) {return _array*=factor;}
  };

  template<typename scalar_t>
  inline auto operator+(const array4_wrapper<scalar_t> &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs._array + rhs._array)
  {
    return (lhs._array + rhs._array);
  }

  template<typename scalar_t, typename other_type>
  inline auto operator+(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) ->decltype(lhs._array + rhs)
  {
    return (lhs._array + rhs);
  }

  template<typename scalar_t, typename other_type>
  inline auto operator+(const other_type &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs + rhs._array)
  {
    return (lhs + rhs._array);
  }

  template<typename scalar_t>
  inline auto operator-(const array4_wrapper<scalar_t> &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs._array - rhs._array)
  {
    return (lhs._array - rhs._array);
  }

  template<typename scalar_t, typename other_type>
  inline auto operator-(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) ->decltype(lhs._array - rhs)
  {
    return (lhs._array - rhs);
  }

  template<typename scalar_t, typename other_type>
  inline auto operator-(const other_type  &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs - rhs._array)
  {
    return (lhs - rhs._array);
  }

  template<typename scalar_t>
  inline auto operator*(const array4_wrapper<scalar_t> &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs._array * rhs._array)
  {
    return (lhs._array * rhs._array);
  }

  template<typename scalar_t, typename other_type>
  inline auto operator*(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) ->decltype(lhs._array * rhs)
  {
    return (lhs._array * rhs);
  }

  template<typename scalar_t, typename other_type>
  inline auto operator*(const other_type &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs * rhs._array)
  {
    return (lhs * rhs._array);
  }

  template<typename scalar_t, typename other_type>
  inline auto operator/(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) ->decltype(lhs._array / rhs)
  {
    return (lhs._array / rhs);
  }

  template<typename scalar_t, typename other_type>
  inline auto operator==(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) {
    return (lhs._array == rhs);
  }

  template<typename scalar_t>
  inline auto operator==(const array4_wrapper<scalar_t> &lhs, const array4_wrapper<scalar_t> &rhs) {
    return (lhs._array == rhs._array);
  }

  // Need a special ==operator for Vc types
  template<typename scalar_t>
  struct alignas(alignment) Vector4 <simd::array4_wrapper<scalar_t>> {
    simd::array4_wrapper<scalar> x, y, z, t;
  };

  template<typename scalar_t>
  bool operator==(const Vector4<simd::array4_wrapper<scalar_t>> &lhs, 
                  const Vector4<simd::array4_wrapper<scalar_t>> &rhs) {
    bool is_eq_x = (lhs.x == rhs.x).isFull();
    bool is_eq_y = (lhs.y == rhs.y).isFull();
    bool is_eq_z = (lhs.z == rhs.z).isFull();
    bool is_eq_t = (lhs.t == rhs.t).isFull();

    return is_eq_x && is_eq_y && is_eq_z && is_eq_t;
  }

  } // namespace simd

} //namespace algebra

