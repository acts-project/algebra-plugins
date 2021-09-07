/** Algebra plugins, part of the ACTS project
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "common/types.hpp"
#include "common/simd_array_wrapper.hpp"

#include <any>
#include <cmath>

// namespace of the algebra object definitions
#ifndef __plugin
#define __plugin algebra::vc_array
// Name of the plugin
#define ALGEBRA_PLUGIN vc_array

#define __plugin_without_matrix_element_accessor 1

namespace algebra
{
    // Define scalar array operator for the 2dim point types
    inline std::array<scalar, 2> operator*(const std::array<scalar, 2> &a, const scalar s)
    {
        return {a[0] * s, a[1] * s};
    }

    inline std::array<scalar, 2> operator*(const scalar s, const std::array<scalar, 2> &a)
    {
        return {s * a[0], s * a[1]};
    }

    inline std::array<scalar, 2> operator/(const std::array<scalar, 2> &a, const scalar s)
    {
        return {a[0] / s, a[1] / s};
    }

    inline std::array<scalar, 2> operator-(const std::array<scalar, 2> &a, const std::array<scalar, 2> &b)
    {
        return {a[0] - b[0], a[1] - b[1]};
    }

    inline std::array<scalar, 2> operator+(const std::array<scalar, 2> &a, const std::array<scalar, 2> &b)
    {
        return {a[0] + b[0], a[1] + b[1]};
    }

    namespace vector
    {
        /** Dot product between two input vectors
         * 
         * @tparam vector_type generic input vector type
         * 
         * @param a the first input vector
         * @param b the second input vector
         * 
         * @return the scalar dot product value 
         **/
        template <typename vector_type>
        inline scalar dot(const vector_type &a, const vector_type &b)
        {
            return (a*b).sum();
        }

        /** Dot product between two input vectors
         *
         * @tparam vec_expr1 vector or vector expression type
         * @tparam vec_expr2 vector or vector expression type
         * 
         * @param a the first input vector/expression
         * @param b the second input vector/expression
         * 
         * @return the scalar dot product value 
         **/
        template <typename vec_expr1, typename vec_expr2>
        inline scalar dot(const vec_expr1 &a, const vec_expr2 &b)
        {
            return (a*b).sum();
        }

        /** Dot product between two input vectors - 2 Dim
         * 
         * @param a the first input vector
         * @param b the second input vector
         * 
         * @return the scalar dot product value 
         **/
        inline scalar dot(const std::array<scalar, 2> &a, const std::array<scalar, 2> &b)
        {
            return (a[0]*b[0] + a[1]*b[1]);
        }

        /** Get a normalized version of the input vector
         * 
         * @tparam vector_type generic input vector type
         *
         * @param v the input vector
         **/
        template <typename vector_type>
        inline vector_type normalize(const vector_type &v)
        {
            return v / std::sqrt(dot(v, v));
        }

        /** Cross product between two input vectors - 3 Dim
         * 
         * @tparam vector_type generic input vector type
         *           
         * @param a the first input vector
         * @param b the second input vector
         * 
         * @return a vector representing the cross product
         **/
        template <typename vector_type>
        inline vector_type cross(const vector_type &a, const vector_type &b)
        {
            return {a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0], a[0] * b[1] - b[0] * a[1], 0};
        }

        /** Cross product between two input vectors - 3 Dim
         * 
         * @tparam vec_expr1 vector or vector expression type
         * @tparam vec_expr2 vector or vector expression type
         *           
         * @param a the first input vector
         * @param b the second input vector
         * 
         * @return a vector (expression) representing the cross product
         **/
        template <typename vec_expr1, typename vec_expr2>
        inline auto cross(const vec_expr1 &a, const vec_expr2 &b) ->decltype(a * b - a * b)
        {
            return {a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0], a[0] * b[1] - b[0] * a[1], 0};
        }

    } // namespace vector

    // array getter methdos
    namespace getter
    {
        /** This method retrieves phi from a vector, vector base with rows > 2
         * 
         * @tparam vector_type generic input vector type
         * 
         * @param v the input vector 
         **/
        template <typename vector_type>
        inline auto phi(const vector_type &v) noexcept
        {
            return std::atan2(v[1], v[0]);
        }

        /** This method retrieves theta from a vector, vector base with rows >= 3
         * 
         * @tparam vector_type generic input vector type
         * 
         * @param v the input vector 
         **/
        template <typename vector_type>
        inline auto theta(const vector_type &v) noexcept
        {
            return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
        }

        /** This method retrieves the perpenticular magnitude of a vector with 2 rows
         * 
         * @tparam vector_type generic input vector type
         * 
         * @param v the input vector 
         **/
        template <typename vector_type>
        inline auto perp(const vector_type &v) noexcept
        {
            return std::sqrt(v[0] * v[0] + v[1] * v[1]);
        }

        /** This method retrieves the norm of a vector, no dimension restriction
         * 
         * @tparam vector_type generic input vector type
         * 
         * @param v the input vector 
         **/
        template <typename vector_type>
        inline auto norm(const vector_type &v)
        {
            return std::sqrt(vector::dot(v, v));
        }

        /** This method retrieves the pseudo-rapidity from a vector with rows >= 3
         * 
         * @tparam vector_type generic input vector type
         * 
         * @param v the input vector 
         **/
        template <typename vector_type>
        inline auto eta(const vector_type &v) noexcept
        {
            return std::atanh(v[2] / norm(v));
        }

        /** This method retrieves a column from a matrix
         * 
         * @tparam matrix_type generic input matrix type
         * 
         * @param m the input matrix 
         **/
        template <unsigned int kROWS = 4, typename matrix_type>
        inline simd::array<scalar, 4> vector(const matrix_type m, unsigned int row, unsigned int col) noexcept
        {
            switch (col) {
                case 2: return m.z;
                case 3: return m.t;
                case 0: return m.x;
                case 1: return m.y;
            }
            return {};
        }

        /** This method retrieves a submatrix
         * 
         * @param m the input matrix 
         **/
        template <unsigned int kROWS, unsigned int kCOLS, typename matrix_type>
        inline auto block(const matrix_type &m, unsigned int row, unsigned int col) noexcept
        {
            simd::array<scalar, 16> submatrix;
            simd::aligned::vector<decltype(m.x)> vecs = {std::move(m.x), std::move(m.y), std::move(m.z), std::move(m.t)};
            for (unsigned int irow = row; irow < row + kROWS; ++irow) {
                for (unsigned int icol = col; icol < col + kCOLS; ++icol) {
                    submatrix[(icol - col) * (irow - row) + (icol - col)] = vecs[irow - row][icol - col];
                }
            }
            return submatrix;
        }

    } // namespace getter

    // array definitions
    namespace vc_array
    {
        // Needed only for correct value initialization
        using vector3 = simd::array4_wrapper<scalar>;
        using point3  = vector3;
        // Don't use vectorization on potentially half-filled vectors
        using vector2 = std::array<scalar, 2>;
        using point2  = vector2;

        /** Transform wrapper class to ensure standard API within differnt plugins
         **/
        struct transform3
        {
            // Keep 4 simd vector for easy handling
            using matrix44 = simd::Vector4<simd::array4_wrapper<scalar>>;

            matrix44 _data;
            matrix44 _data_inv;

            /** Contructor with arguments: t, z, x
             * 
             * @param t the translation (or origin of the new frame)
             * @param z the z axis of the new frame, normal vector for planes
             * @param x the x axis of the new frame
             * 
             * @note y will be constructed by cross product
             * 
             **/
            transform3(const vector3 &t, const vector3 &z, const vector3 &x)
            {
                auto y = vector::cross(z, x);
                _data.x = {x[0], x[1], x[2], 0.};
                _data.y = {y[0], y[1], y[2], 0.};
                _data.z = {z[0], z[1], z[2], 0.};
                _data.t = {t[0], t[1], t[2], 1.};

                _data_inv = invert(_data);
            }

            /** Constructor with arguments: translation
             *
             * @param t is the transform
             **/
            transform3(const vector3 &t)
            {
                _data.x = {1., 0., 0., 0.};
                _data.y = {0., 1., 0., 0.};
                _data.z = {0., 0., 1., 0.};
                _data.t = {t[0], t[1], t[2], 1.};

                _data_inv = invert(_data);
            }

            /** Constructor with arguments: matrix 
             * 
             * @param m is the full 4x4 matrix 
             **/
            transform3(const matrix44 &m)
            {
                _data = m;
            }

            /** Constructor with arguments: matrix as std::aray of scalar
             * 
             * @param ma is the full 4x4 matrix 16 array
             **/
            transform3(const array_s<scalar, 16> &ma)
            {
                _data.x = {ma[0], ma[4], ma[8], ma[12]};
                _data.y = {ma[1], ma[5], ma[9], ma[13]};
                _data.z = {ma[2], ma[6], ma[10], ma[14]};
                _data.t = {ma[3], ma[7], ma[11], ma[15]};

                _data_inv = invert(_data);
            }

            /** Constructor with arguments: identity
             *
             **/
            transform3()
            {
                _data.x = {1., 0., 0., 0.};
                _data.y = {0., 1., 0., 0.};
                _data.z = {0., 0., 1., 0.};
                _data.t = {0., 0., 0., 1.};

                _data_inv = _data;
            }

            /** Default contructors */
            transform3(const transform3 &rhs) = default;
            ~transform3() = default;

            /** Equality operator */
            inline bool operator==(const transform3 &rhs) const
            {
                return (_data == rhs._data);
            }

            /** The determinant of a 4x4 matrix
             * 
             * @param m is the matrix
             *
             * @return a sacalar determinant - no checking done 
             */
            static inline scalar determinant(const matrix44 &m)
            {
                return m.t[0] * m.z[1] * m.y[2] * m.x[3] - m.z[0] * m.t[1] * m.y[2] * m.x[3] - m.t[0] * m.y[1] * m.z[2] * m.x[3] + m.y[0] * m.t[1] * m.z[2] * m.x[3] +
                       m.z[0] * m.y[1] * m.t[2] * m.x[3] - m.y[0] * m.z[1] * m.t[2] * m.x[3] - m.t[0] * m.z[1] * m.x[2] * m.y[3] + m.z[0] * m.t[1] * m.x[2] * m.y[3] +
                       m.t[0] * m.x[1] * m.z[2] * m.y[3] - m.x[0] * m.t[1] * m.z[2] * m.y[3] - m.z[0] * m.x[1] * m.t[2] * m.y[3] + m.x[0] * m.z[1] * m.t[2] * m.y[3] +
                       m.t[0] * m.y[1] * m.x[2] * m.z[3] - m.y[0] * m.t[1] * m.x[2] * m.z[3] - m.t[0] * m.x[1] * m.y[2] * m.z[3] + m.x[0] * m.t[1] * m.y[2] * m.z[3] +
                       m.y[0] * m.x[1] * m.t[2] * m.z[3] - m.x[0] * m.y[1] * m.t[2] * m.z[3] - m.z[0] * m.y[1] * m.x[2] * m.t[3] + m.y[0] * m.z[1] * m.x[2] * m.t[3] +
                       m.z[0] * m.x[1] * m.y[2] * m.t[3] - m.x[0] * m.z[1] * m.y[2] * m.t[3] - m.y[0] * m.x[1] * m.z[2] * m.t[3] + m.x[0] * m.y[1] * m.z[2] * m.t[3];
            }

            /** The inverse of a 4x4 matrix
             * 
             * @param m is the matrix
             *
             * @return an inverse matrix 
             */
            static inline matrix44 invert(const matrix44 &m)
            {
                matrix44 i;
                i.x[0] = m.z[1] * m.t[2] * m.y[3] - m.t[1] * m.z[2] * m.y[3] + m.t[1] * m.y[2] * m.z[3] - m.y[1] * m.t[2] * m.z[3] - m.z[1] * m.y[2] * m.t[3] + m.y[1] * m.z[2] * m.t[3];
                i.x[1] = m.t[1] * m.z[2] * m.x[3] - m.z[1] * m.t[2] * m.x[3] - m.t[1] * m.x[2] * m.z[3] + m.x[1] * m.t[2] * m.z[3] + m.z[1] * m.x[2] * m.t[3] - m.x[1] * m.z[2] * m.t[3];
                i.x[2] = m.y[1] * m.t[2] * m.x[3] - m.t[1] * m.y[2] * m.x[3] + m.t[1] * m.x[2] * m.y[3] - m.x[1] * m.t[2] * m.y[3] - m.y[1] * m.x[2] * m.t[3] + m.x[1] * m.y[2] * m.t[3];
                i.x[3] = m.z[1] * m.y[2] * m.x[3] - m.y[1] * m.z[2] * m.x[3] - m.z[1] * m.x[2] * m.y[3] + m.x[1] * m.z[2] * m.y[3] + m.y[1] * m.x[2] * m.z[3] - m.x[1] * m.y[2] * m.z[3];
                i.y[0] = m.t[0] * m.z[2] * m.y[3] - m.z[0] * m.t[2] * m.y[3] - m.t[0] * m.y[2] * m.z[3] + m.y[0] * m.t[2] * m.z[3] + m.z[0] * m.y[2] * m.t[3] - m.y[0] * m.z[2] * m.t[3];
                i.y[1] = m.z[0] * m.t[2] * m.x[3] - m.t[0] * m.z[2] * m.x[3] + m.t[0] * m.x[2] * m.z[3] - m.x[0] * m.t[2] * m.z[3] - m.z[0] * m.x[2] * m.t[3] + m.x[0] * m.z[2] * m.t[3];
                i.y[2] = m.t[0] * m.y[2] * m.x[3] - m.y[0] * m.t[2] * m.x[3] - m.t[0] * m.x[2] * m.y[3] + m.x[0] * m.t[2] * m.y[3] + m.y[0] * m.x[2] * m.t[3] - m.x[0] * m.y[2] * m.t[3];
                i.y[3] = m.y[0] * m.z[2] * m.x[3] - m.z[0] * m.y[2] * m.x[3] + m.z[0] * m.x[2] * m.y[3] - m.x[0] * m.z[2] * m.y[3] - m.y[0] * m.x[2] * m.z[3] + m.x[0] * m.y[2] * m.z[3];
                i.z[0] = m.z[0] * m.t[1] * m.y[3] - m.t[0] * m.z[1] * m.y[3] + m.t[0] * m.y[1] * m.z[3] - m.y[0] * m.t[1] * m.z[3] - m.z[0] * m.y[1] * m.t[3] + m.y[0] * m.z[1] * m.t[3];
                i.z[1] = m.t[0] * m.z[1] * m.x[3] - m.z[0] * m.t[1] * m.x[3] - m.t[0] * m.x[1] * m.z[3] + m.x[0] * m.t[1] * m.z[3] + m.z[0] * m.x[1] * m.t[3] - m.x[0] * m.z[1] * m.t[3];
                i.z[2] = m.y[0] * m.t[1] * m.x[3] - m.t[0] * m.y[1] * m.x[3] + m.t[0] * m.x[1] * m.y[3] - m.x[0] * m.t[1] * m.y[3] - m.y[0] * m.x[1] * m.t[3] + m.x[0] * m.y[1] * m.t[3];
                i.z[3] = m.z[0] * m.y[1] * m.x[3] - m.y[0] * m.z[1] * m.x[3] - m.z[0] * m.x[1] * m.y[3] + m.x[0] * m.z[1] * m.y[3] + m.y[0] * m.x[1] * m.z[3] - m.x[0] * m.y[1] * m.z[3];
                i.t[0] = m.t[0] * m.z[1] * m.y[2] - m.z[0] * m.t[1] * m.y[2] - m.t[0] * m.y[1] * m.z[2] + m.y[0] * m.t[1] * m.z[2] + m.z[0] * m.y[1] * m.t[2] - m.y[0] * m.z[1] * m.t[2];
                i.t[1] = m.z[0] * m.t[1] * m.x[2] - m.t[0] * m.z[1] * m.x[2] + m.t[0] * m.x[1] * m.z[2] - m.x[0] * m.t[1] * m.z[2] - m.z[0] * m.x[1] * m.t[2] + m.x[0] * m.z[1] * m.t[2];
                i.t[2] = m.t[0] * m.y[1] * m.x[2] - m.y[0] * m.t[1] * m.x[2] - m.t[0] * m.x[1] * m.y[2] + m.x[0] * m.t[1] * m.y[2] + m.y[0] * m.x[1] * m.t[2] - m.x[0] * m.y[1] * m.t[2];
                i.t[3] = m.y[0] * m.z[1] * m.x[2] - m.z[0] * m.y[1] * m.x[2] + m.z[0] * m.x[1] * m.y[2] - m.x[0] * m.z[1] * m.y[2] - m.y[0] * m.x[1] * m.z[2] + m.x[0] * m.y[1] * m.z[2];
                scalar idet = 1. / determinant(i);

                i.x *= idet;
                i.y *= idet;
                i.z *= idet;
                i.t *= idet;

                return i;
            }

            /** Rotate a vector into / from a frame 
             * 
             * @param m is the rotation matrix
             * @param v is the vector to be rotated
             */
            static inline vector3 rotate(const matrix44 &m, const vector3 &v)
            {
                return m.x*v[0] + m.y*v[1] + m.z*v[2];
            }

            /** This method retrieves the rotation of a transform */
            inline auto rotation() const
            {
                return getter::block<3, 3>(_data, 0, 0);
            }

            /** This method retrieves the translation of a transform */
            inline point3 translation() const
            {
                return _data.t;
            }

            /** This method retrieves the 4x4 matrix of a transform */
            inline const matrix44 &matrix() const
            {
                return _data;
            }

            /** This method transform from a point from the local 3D cartesian frame 
             *  to the global 3D cartesian frame 
             *
             * @tparam point_type 3D point
             *
             * @param v is the point to be transformed
             *
             * @return a global point
             */
            template <typename point_type>
            inline const point_type point_to_global(const point_type &v) const
            {
                return _data.x*v[0] + _data.y*v[1] + _data.z*v[2] + _data.t;
            }

            /** This method transform from a vector from the global 3D cartesian frame 
             *  into the local 3D cartesian frame
             *
             * @tparam point_type 3D point
             *
             * @param v is the point to be transformed
             *
             * @return a local point
             */
            template <typename point_type>
            inline const point_type point_to_local(const point_type &v) const
            {
                return _data_inv.x*v[0] + _data_inv.y*v[1] + _data_inv.z*v[2] + _data_inv.t;
            }

            /** This method transform from a vector from the local 3D cartesian frame 
             *  to the global 3D cartesian frame
             *
             * @tparam vector_type 3D vector 
             *
             * @param v is the vector to be transformed
             *
             * @return a vector in global coordinates
             */
            template <typename vector_type>
            inline const vector_type vector_to_global(const vector_type &v) const
            {
                return rotate(_data, v);
            }

            /** This method transform from a vector from the global 3D cartesian frame
             *  into the local 3D cartesian frame
             *
             * @tparam vector_type 3D vector 
             *
             * @param v is the vector to be transformed
             *
             * @return a vector in global coordinates
             */
            template <typename vector_type>
            inline const auto vector_to_local(const vector_type &v) const
            {
                return rotate(_data_inv, v);
            }
        };

        /** Frame projection into a cartesian coordinate frame
         */
        struct cartesian2
        {
            /** This method transform from a point from the global 3D cartesian frame 
             *  to the local 2D cartesian frame 
             * 
             * @param trf the transform from global to local thredimensional frame
             * @param p the point in global frame
             * 
             * @return a local point2
             **/
            inline point2 operator()(const transform3 &trf,
                                  const point3 &p) const
            {
                return operator()(trf.point_to_local(p));
            }

            /** This method transform from a point from the global 3D cartesian 
             *  frame to the local 2D cartesian frame
             *
             * @param v the point in local frame
             * 
             * @return a local point2
             */
            inline point2 operator()(const point3 &v) const
            {
                return {v[0], v[1]};
            }
        };

        /** Local frame projection into a polar coordinate frame */
        struct polar2
        {
            /** This method transform from a point from the global 3D cartesian 
             *  frame to the local 2D cartesian frame 
             * 
             * @param trf the transform from global to local thredimensional frame
             * @param p the point in global frame
             * 
             * @return a local point2
             **/
            inline point2 operator()(const transform3 &trf,
                                  const point3 &p) const
            {
                return operator()(trf.point_to_local(p));
            }

            /** This method transform from a point from 2D or 3D cartesian frame 
                to a 2D polar point */
            template <typename point_type>
            inline point2 operator()(const point_type &v) const
            {
                return point2{getter::perp(v), getter::phi(v)};
            }
        };

        /** Local frame projection into a polar coordinate frame */
        struct cylindrical2
        {
             /** This method transform from a point from the global 3D cartesian 
              *  frame to the local 2D cartesian frame 
              * 
              * @param trf the transform from global to local thredimensional frame
              * @param p the point in global frame
              * 
              * @return a local point2
              **/
            inline point2 operator()(const transform3 &trf,
                                  const point3 &p) const
            {
                return operator()(trf.point_to_local(p));
            }

            /** This method transform from a point from 2 3D cartesian frame to 
                a 2D cylindrical point */
            inline point2 operator()(const point3 &v) const
            {
                return point2{getter::perp(v) * getter::phi(v), v[2]};
            }
        };

    } // namespace vc_array

} // namespace algebra

#endif
