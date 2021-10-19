/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "algebra/common/types.hpp"
#include "algebra/common/algebra_qualifiers.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <any>
#include <tuple>
#include <cmath>

#ifdef ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE
using algebra_scalar = ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE;
#else
using algebra_scalar = double;
#endif

// namespace of the algebra object definitions
#ifndef __plugin
#define __plugin algebra::eigen
// Name of the plugin
#define ALGEBRA_PLUGIN eigen

namespace algebra
{
    using scalar = algebra_scalar;

    // eigen getter methdos
    namespace getter
    {
        /** This method retrieves phi from a vector, vector base with rows > 2
         *
         * @param v the input vector
         **/
        template <typename derived_type>
	ALGEBRA_HOST_DEVICE	
	inline auto phi(const Eigen::MatrixBase<derived_type> &v) noexcept
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            static_assert(rows >= 2, "vector::phi() required rows >= 2.");
            return std::atan2(v[1], v[0]);
        }

        /** This method retrieves theta from a vector, vector base with rows >= 3
         *
         * @param v the input vector
         **/
        template <typename derived_type>
	ALGEBRA_HOST_DEVICE	
        inline auto theta(const Eigen::MatrixBase<derived_type> &v) noexcept
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            static_assert(rows >= 2, "vector::theta() required rows >= 3.");
            return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
        }

        /** This method retrieves the pseudo-rapidity from a vector or vector base with rows >= 3
         *
         * @param v the input vector
         **/
        template <typename derived_type>
	ALGEBRA_HOST_DEVICE	
        inline auto eta(const Eigen::MatrixBase<derived_type> &v) noexcept
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            static_assert(rows >= 2, "vector::eta() required rows >= 3.");
            return std::atanh(v[2] / v.norm());
        }

        /** This method retrieves the perpenticular magnitude of a vector with rows >= 2
         *
         * @param v the input vector
         **/
        template <typename derived_type>
	ALGEBRA_HOST_DEVICE	
        inline auto perp(const Eigen::MatrixBase<derived_type> &v) noexcept
        {
            constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
            static_assert(rows >= 2, "vector::perp() required rows >= 2.");
            return std::sqrt(v[0] * v[0] + v[1] * v[1]);
        }

        /** This method retrieves the norm of a vector, no dimension restriction
         *
         * @param v the input vector
         **/
        template <typename derived_type>
	ALGEBRA_HOST_DEVICE	
        inline auto norm(const Eigen::MatrixBase<derived_type> &v)
        {
            return v.norm();
        }

        /** This method retrieves a column from a matrix
         *
         * @param m the input matrix
         **/
        template <unsigned int kROWS, typename derived_type>
	ALGEBRA_HOST_DEVICE	
        inline auto vector(const Eigen::MatrixBase<derived_type> &m, unsigned int row, unsigned int col)
        {
            return m.template block<kROWS, 1>(row, col);
        }

        /** This method retrieves a column from a matrix
         *
         * @param m the input matrix
         **/
        template <unsigned int kROWS, unsigned int kCOLS, typename derived_type>
	ALGEBRA_HOST_DEVICE	
        inline auto block(const Eigen::MatrixBase<derived_type> &m, unsigned int row, unsigned int col)
        {
            return m.template block<kROWS, kCOLS>(row, col);
        }

    } // namespace getter

    // eigen definitions
    namespace eigen
    {
        using vector3 = Eigen::Matrix<scalar, 3, 1>;
        using point3  = vector3;
        using point2  = Eigen::Matrix<scalar, 2, 1>;

        /** Transform wrapper class to ensure standard API within differnt plugins */
        struct transform3
        {
            Eigen::Transform<scalar, 3, Eigen::Affine> _data =
                Eigen::Transform<scalar, 3, Eigen::Affine>::Identity();

            Eigen::Transform<scalar, 3, Eigen::Affine> _data_inv =
                Eigen::Transform<scalar, 3, Eigen::Affine>::Identity();

            using matrix44 = Eigen::Transform<scalar, 3, Eigen::Affine>::MatrixType;

            /** Contructor with arguments: t, z, x
             *
             * @param t the translation (or origin of the new frame)
             * @param z the z axis of the new frame, normal vector for planes
             * @param x the x axis of the new frame
             *
             **/
	    ALGEBRA_HOST_DEVICE
            transform3(const vector3 &t, const vector3 &z, const vector3 &x)
            {
                auto y = z.cross(x);

                auto &matrix = _data.matrix();
                matrix.block<3, 1>(0, 0) = x;
                matrix.block<3, 1>(0, 1) = y;
                matrix.block<3, 1>(0, 2) = z;
                matrix.block<3, 1>(0, 3) = t;

                _data_inv = _data.inverse();
            }

            /** Constructor with arguments: translation
             *
             * @param t is the transform
             **/
	    ALGEBRA_HOST_DEVICE
            transform3(const vector3 &t)
            {
                auto &matrix = _data.matrix();
                matrix.block<3, 1>(0, 3) = t;

                _data_inv = _data.inverse();
            }

            /** Constructor with arguments: matrix
             *
             * @param m is the full 4x4 matrix
             **/
	    ALGEBRA_HOST_DEVICE
            transform3(const matrix44 &m)
            {
                _data.matrix() = m;

                _data_inv = _data.inverse();
            }

            /** Constructor with arguments: matrix as std::aray of scalar
             *
             * @param ma is the full 4x4 matrix asa 16 array
             **/
	    ALGEBRA_HOST_DEVICE
            transform3(const array_s<scalar, 16> &ma)
            {
                _data.matrix() << ma[0], ma[1], ma[2], ma[3], ma[4], ma[5], ma[6], ma[7],
                    ma[8], ma[9], ma[10], ma[11], ma[12], ma[13], ma[14], ma[15];

                _data_inv = _data.inverse();
            }

            /** Default contructors */
            transform3() = default;
            transform3(const transform3 &rhs) = default;
            ~transform3() = default;

            /** Equality operator */
	    ALGEBRA_HOST_DEVICE
            inline bool operator==(const transform3 &rhs) const
            {
                return (_data.isApprox(rhs._data));
            }

            /** Rotate a vector into / from a frame
             *
             * @param m is the rotation matrix
             * @param v is the vector to be rotated
             */
            template <typename derived_type>
	    ALGEBRA_HOST_DEVICE
            static inline auto rotate(const Eigen::Transform<scalar, 3, Eigen::Affine> &m, const Eigen::MatrixBase<derived_type> &v)
            {
                constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
                constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
                static_assert(rows == 3 and cols == 1, "transform::rotate(m,v) requires a (3,1) matrix");
                return m.matrix().block<3,3>(0,0)*v;
            }

            /** This method retrieves the rotation of a transform  **/
	    ALGEBRA_HOST_DEVICE
            inline auto rotation() const
            {
                return _data.matrix().block<3, 3>(0, 0);
            }

            /** This method retrieves the translation of a transform **/
	    ALGEBRA_HOST_DEVICE
            inline auto translation() const
            {
                return _data.matrix().block<3, 1>(0, 3);
            }

            /** This method retrieves the 4x4 matrix of a transform */
	    ALGEBRA_HOST_DEVICE
            inline const auto &matrix() const
            {
                return _data.matrix();
            }

            /** This method transform from a point from the local 3D cartesian frame to the global 3D cartesian frame */
            template <typename derived_type>
	    ALGEBRA_HOST_DEVICE	    
            inline auto point_to_global(const Eigen::MatrixBase<derived_type> &v) const
            {
                constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
                constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
                static_assert(rows == 3 and cols == 1, "transform::point_to_global(v) requires a (3,1) matrix");
                return (_data * v);
            }

            /** This method transform from a vector from the global 3D cartesian frame into the local 3D cartesian frame */
            template <typename derived_type>
	    ALGEBRA_HOST_DEVICE	    
            inline auto point_to_local(const Eigen::MatrixBase<derived_type> &v) const
            {
                constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
                constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
                static_assert(rows == 3 and cols == 1, "transform::point_to_local(v) requires a (3,1) matrix");
                return (_data_inv * v);
            }

            /** This method transform from a vector from the local 3D cartesian frame to the global 3D cartesian frame */
            template <typename derived_type>
	    ALGEBRA_HOST_DEVICE
            inline auto vector_to_global(const Eigen::MatrixBase<derived_type> &v) const
            {
                constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
                constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
                static_assert(rows == 3 and cols == 1, "transform::vector_to_global(v) requires a (3,1) matrix");
                return (_data.linear() * v);
            }

            /** This method transform from a vector from the global 3D cartesian frame into the local 3D cartesian frame */
            template <typename derived_type>
	    ALGEBRA_HOST_DEVICE
            inline auto vector_to_local(const Eigen::MatrixBase<derived_type> &v) const
            {
                constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
                constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
                static_assert(rows == 3 and cols == 1, "transform::vector_to_local(v) requires a (3,1) matrix");
                return (_data_inv.linear() * v);
            }
        };

        /** Local frame projection into a cartesian coordinate frame
         */
        struct cartesian2
        {
            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             *
             * @param v the point in local frame
             *
             * @return a local point2
             */
            template <typename derived_type>
	    ALGEBRA_HOST_DEVICE	    
            inline auto operator()(const Eigen::MatrixBase<derived_type> &v) const
            {
                constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
                constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
                static_assert(rows == 3 and cols == 1, "transform::point3_topoint2(v) requires a (3,1) matrix");
                return (v.template segment<2>(0)).eval();
            }

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             *
             * @param trf the transform from global to local thredimensional frame
             * @param p the point in global frame
             *
             * @return a local point2
             **/
	    ALGEBRA_HOST_DEVICE	    
            inline auto operator()(const transform3 &trf,
                            const point3 &p) const
            {
                return operator()(trf.point_to_local(p));
            }
        };

        /** Local frame projection into a polar coordinate frame */
        struct polar2
        {
            /** This method transform from a point from 2D or 3D cartesian frame to a 2D polar point
             *
             * @param v the point in local frame
             *
             * @return a local point2
             */
            template <typename derived_type>
	    ALGEBRA_HOST_DEVICE	    
            inline auto operator()(const Eigen::MatrixBase<derived_type> &v) const
            {
                constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
                constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
                static_assert(rows >= 2 and cols == 1, "transform::point_topoint2pol(v) requires a (>2,1) matrix");
                return point2{getter::perp(v), getter::phi(v)};
            }

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             *
             * @param trf the transform from global to local thredimensional frame
             * @param p the point in global frame
             *
             * @return a local point2
             **/
	    ALGEBRA_HOST_DEVICE
            inline auto operator()(const transform3 &trf,
				   const point3 &p) const
            {
                return operator()(trf.point_to_local(p));
            }
        };

        /** Local frame projection into a polar coordinate frame */
        struct cylindrical2
        {
            /** This method transform from a point from 2D or 3D cartesian frame to a 2D polar point
             *
             * @param v the point in local frame
             *
             * @return a local point2
             */
            template <typename derived_type>
	    ALGEBRA_HOST_DEVICE	    
            inline auto operator()(const Eigen::MatrixBase<derived_type> &v) const
            {

                constexpr int rows = Eigen::MatrixBase<derived_type>::RowsAtCompileTime;
                constexpr int cols = Eigen::MatrixBase<derived_type>::ColsAtCompileTime;
                static_assert(rows == 3 and cols == 1, "transform::point3_topoint2cyl(v) requires a a (3,1) matrix");
                return point2{getter::perp(v) * getter::phi(v), v[2]};
            }

            /** This method transform from a point from the global 3D cartesian frame to the local 2D cartesian frame
             *
             * @param trf the transform from global to local thredimensional frame
             * @param p the point in global frame
             *
             * @return a local point2
             **/
	    ALGEBRA_HOST_DEVICE
            inline auto operator()(const transform3 &trf,
				   const point3 &p) const
            {
                return operator()(trf.point_to_local(p));
            }
        };

    } // namespace eigen

    // Vector transfroms
    namespace vector
    {

        /** Get a normalized version of the input vector
         *
         * @tparam derived_type is the matrix template
         *
         * @param v the input vector
         **/
        template <typename derived_type>
	ALGEBRA_HOST_DEVICE
        inline auto normalize(const Eigen::MatrixBase<derived_type> &v)
        {
            return v.normalized();
        }

        /** Dot product between two input vectors
         *
         * @tparam derived_type_lhs is the first matrix (epresseion) template
         * @tparam derived_type_rhs is the second matrix (epresseion) template
         *
         * @param a the first input vector
         * @param b the second input vector
         *
         * @return the scalar dot product value
         **/
        template <typename derived_type_lhs, typename derived_type_rhs>
	ALGEBRA_HOST_DEVICE	
        inline auto dot(const Eigen::MatrixBase<derived_type_lhs> &a, const Eigen::MatrixBase<derived_type_rhs> &b)
        {
            return a.dot(b);
        }

        /** Cross product between two input vectors
         *
         * @tparam derived_type_lhs is the first matrix (epresseion) template
         * @tparam derived_type_rhs is the second matrix (epresseion) template
         *
         * @param a the first input vector
         * @param b the second input vector
         *
         * @return a vector (expression) representing the cross product
         **/
        template <typename derived_type_lhs, typename derived_type_rhs>
	ALGEBRA_HOST_DEVICE	
        inline auto cross(const Eigen::MatrixBase<derived_type_lhs> &a, const Eigen::MatrixBase<derived_type_rhs> &b)
        {
            return a.cross(b);
        }
    } // namespace vector

} // namespace algebra

#endif
