/** Algebra plugins, part of the ACTS project
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
 
#pragma once

#include "common/types.hpp"

#include <array>
#include <cmath>

// namespace of the algebra object definitions
#ifndef __plugin
#define __plugin algebra::std_array
// Name of the plugin
#define ALGEBRA_PLUGIN array

// Use the std array type
#define __plugin_array std::array
#include "common/definitions/array.hpp"

namespace algebra {
    // array definitions
    namespace std_array
    {
        using vector3 = std::array<scalar, 3>;
        using point3 = vector3;
        using point2 = std::array<scalar, 2>;

        /** Transform wrapper class to ensure standard API within differnt plugins
         **/
        using transform3 = array::transform3;

        using cartesian2 = array::cartesian2;
        using polar2 = array::polar2;
        using cylindrical2 = array::cylindrical2;

    } // namespace std_array
} // namespace algebra

#undef __plugin_array
#endif
