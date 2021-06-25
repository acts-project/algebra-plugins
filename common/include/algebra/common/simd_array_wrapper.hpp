/** Algebra plugins, part of the ACTS project
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "common/simd_types.hpp"

#include <Vc/Vc>

namespace algebra {

namespace simd {
  
  /** Wrapper structure around the Vc::SimdArray class, that allows correct initialization
   *  from std::initializer_list with a different length than the size of the array.
   *
   * @tparam scalar_t scalar type used in the SimdArray
   */
  template<typename scalar_t>
  struct array4_wrapper {

    // Instance of SimdArray that is being wrapped
    simd::array<scalar_t, 4> _array;

    // Mask type that is returned from boolean operations on vector types
    using mask_type = typename decltype(_array)::mask_type;
    
    /** Default constructor
     */
    array4_wrapper() : _array() {}

    /** const copy constructor
     * 
     * @param array wrapper object to copy the data from
     */
    array4_wrapper(const array4_wrapper& array) : _array(array._array) {}

    /** Copy constructor
     * 
     * @param array wrapper object to copy the data from
     */
    array4_wrapper(array4_wrapper&& array) : _array(std::move(array._array)) {}

    /** Initialization from a single value. The Vc type broadcasts this into 
     *  vector data.
     *
     * @param value The value type to construct an array from
     */
    array4_wrapper(scalar_t value) : _array(value) {}

    /** Parametrized constructor that takes a number of values and makes sure that
     *  also all remaining values in the array are properly initialized.
     *
     * @param v1 value 1. dimension
     * @param v1 value 2. dimension
     * @param v1 value 3. dimension
     * @param v1 optional value 4. dimension. if not passed, set to 0.0
     *
     */
    array4_wrapper(scalar_t v1, scalar_t v2, scalar_t v3, scalar_t v4 = scalar_t{0.0}) 
      : _array{v1, v2, v3, v4} {}

    /** Constructor from wrapped class. Used for conversions into wrapped type.
     *
     * @param base data to be wrapped
     */
    array4_wrapper(simd::array<scalar_t, 4> base) : _array(std::move(base)) {}

    /** Generig assignment operator
     * 
     * @param lhs wrap a copy of this data
     */
    template<typename other_type>
    inline const array4_wrapper<scalar_t> & operator=(const other_type& lhs)
    {
        _array = lhs;
        return *this;
    }

    /** Assignment operator from another wrapper
     * 
     * @param other wrap a copy of this wrapped data
     */
    inline array4_wrapper<scalar_t>& operator=(const array4_wrapper<scalar_t>& other)
    {
        _array = other._array;
        return *this;
    }

    /** Assignment operator from std::initializer_list
     * 
     * @param list wrap an array of this data
     */
    inline array4_wrapper<scalar_t>& operator=(std::initializer_list<scalar_t> &list)
    {
        _array = simd::array<scalar_t, 4>(list);
        return *this;
    }

    /** Conversion operator from wrapper to SimdArray.
     */
    operator simd::array<scalar_t, 4>&() { return _array; }
    operator const simd::array<scalar_t, 4>&() const { return _array; }

    /** Operator[] overload from SimdArray for simd array wrapper.
     *
     * @return Value at given index
     */
    inline auto operator[](const unsigned int i) {return _array[i];}
    inline auto operator[](const unsigned int i) const {return _array[i];}

    /** Operator*= overload from SimdArray for simd array wrapper.
     *
     * @return Vector expression/ return type according to the operation
     */
    inline auto operator*=(scalar_t factor) ->decltype(_array*=factor) 
    {
      return _array*=factor;
    }
  };

  /** Operator overload from SimdArray for simd array wrapper that handles the default
   *  floating point literal conversion.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator+(const array4_wrapper<scalar_t> &lhs, const double &rhs) /*->decltype(lhs._array + scalar_t{rhs})*/
  {
    return (lhs._array + scalar_t{rhs});
  }

  /** Operator overload from SimdArray for simd array wrapper that handles the default
   *  floating point literal conversion.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator+(const double &lhs, const array4_wrapper<scalar_t> &rhs) /*->decltype(scalar_t{lhs} + rhs._array)*/
  {
    return (scalar_t{lhs} + rhs._array);
  }
  
  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator+(const array4_wrapper<scalar_t> &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs._array + rhs._array)
  {
    return (lhs._array + rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t, typename other_type>
  inline auto operator+(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) ->decltype(lhs._array + rhs)
  {
    return (lhs._array + rhs);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t, typename other_type>
  inline auto operator+(const other_type &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs + rhs._array)
  {
    return (lhs + rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper that handles the default
   *  floating point literal conversion.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator-(const array4_wrapper<scalar_t> &lhs, const double &rhs) /*->decltype(lhs._array - scalar_t{rhs})*/
  {
    return (lhs._array - scalar_t{rhs});
  }

  /** Operator overload from SimdArray for simd array wrapper that handles the default
   *  floating point literal conversion.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator-(const double &lhs, const array4_wrapper<scalar_t> &rhs) /*->decltype(scalar_t{lhs} - rhs._array)*/
  {
    return (scalar_t{lhs} - rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator-(const array4_wrapper<scalar_t> &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs._array - rhs._array)
  {
    return (lhs._array - rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t, typename other_type>
  inline auto operator-(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) ->decltype(lhs._array - rhs)
  {
    return (lhs._array - rhs);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t, typename other_type>
  inline auto operator-(const other_type &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs - rhs._array)
  {
    return (lhs- rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper that handles the default
   *  floating point literal conversion.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator*(const array4_wrapper<scalar_t> &lhs, const double &rhs) /*->decltype(lhs._array * scalar_t{rhs})*/
  {
    return (lhs._array * scalar_t{rhs});
  }

  /** Operator overload from SimdArray for simd array wrapper that handles the default
   *  floating point literal conversion.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator*(const double &lhs, const array4_wrapper<scalar_t> &rhs) /*->decltype(scalar_t{lhs} * rhs._array)*/
  {
    return (scalar_t{lhs} * rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator*(const array4_wrapper<scalar_t> &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs._array * rhs._array)
  {
    return (lhs._array * rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t, typename other_type>
  inline auto operator*(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) ->decltype(lhs._array * rhs)
  {
    return (lhs._array * rhs);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t, typename other_type>
  inline auto operator*(const other_type &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs * rhs._array)
  {
    return (lhs * rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper that handles the default
   *  floating point literal conversion.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator/(const array4_wrapper<scalar_t> &lhs, const double &rhs) /*->decltype(lhs._array / scalar_t{rhs})*/
  {
    return (lhs._array / scalar_t{rhs});
  }

  /** Operator overload from SimdArray for simd array wrapper that handles the default
   *  floating point literal conversion.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator/(const double &lhs, const array4_wrapper<scalar_t> &rhs) /*->decltype(scalar_t{lhs} / rhs._array)*/
  {
    return (scalar_t{lhs} / rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t custom scalar type used in the wrapper
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator/(const array4_wrapper<scalar_t> &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs._array / rhs._array)
  {
    return (lhs._array / rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t, typename other_type>
  inline auto operator/(const other_type &lhs, const array4_wrapper<scalar_t> &rhs) ->decltype(lhs / rhs._array)
  {
    return (lhs / rhs._array);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t, typename other_type>
  inline auto operator/(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) ->decltype(lhs._array / rhs)
  {
    return (lhs._array / rhs);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t, typename other_type>
  inline auto operator==(const array4_wrapper<scalar_t> &lhs, const other_type &rhs) {
    return (lhs._array == rhs);
  }

  /** Operator overload from SimdArray for simd array wrapper.
   *
   * @tparam scalar_t value type used in the wrapper
   * @tparam other_type generic type of the second operand
   *
   * @return Vector expression/ return type according to the operation
   */
  template<typename scalar_t>
  inline auto operator==(const array4_wrapper<scalar_t> &lhs, const array4_wrapper<scalar_t> &rhs) {
    return (lhs._array == rhs._array);
  }

  /** Equality operator definition for simd array wrapper in a simd Vector4 struct.
   *
   * @tparam scalar_t value type used in the wrapper
   *
   * @return Result of comparision
   */
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

