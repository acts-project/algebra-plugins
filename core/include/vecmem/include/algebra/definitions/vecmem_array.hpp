/** Algebra plugins, part of the ACTS project
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
 
#pragma once

#include "common/types.hpp"

#include <vecmem/containers/static_array.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

#include <array>
#include <cmath>

// namespace of the algebra object definitions
#define __plugin algebra::vecmem_array
// Name of the plugin
#define ALGEBRA_PLUGIN vecmem_array
// Use the vecmem array type
#define __plugin_array vecmem::static_array
#include "common/definitions/array.hpp"

namespace algebra {
    // array definitions
    namespace vecmem_array
    {
        using vector3 = vecmem::static_array<scalar, 3>;
        using point3 = vector3;
        using point2 = vecmem::static_array<scalar, 2>;

        /** Transform wrapper class to ensure standard API within differnt plugins
         **/
        using transform3 = algebra::array::transform3;

        using cartesian2 = algebra::array::cartesian2;
        using polar2 = algebra::array::polar2;
        using cylindrical2 = algebra::array::cylindrical2;

    } // namespace std_array
} // namespace algebra
