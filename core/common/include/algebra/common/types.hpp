/** Algebra plugins, part of the ACTS project
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>
#include <vector>
#include <map>
#include <tuple>

#ifdef ALGEBRA_PLUGIN_USE_VECMEM
#include <vecmem/containers/vector.hpp>
#include <vecmem/containers/static_array.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#endif

#ifdef ALGEBRA_PLUGIN_INCLUDE_VC
#include "simd_types.hpp"
#endif

namespace algebra
{
    // General algebra container types to be used with plugin
    #ifdef ALGEBRA_PLUGIN_USE_VECMEM

    template <typename value_type, unsigned int kDIM>
    using array_s = vecmem::static_array<value_type, kDIM>;

    template <typename value_type>
    using vector_s = vecmem::vector<value_type>;

    #else

    template <typename value_type, unsigned int kDIM>
    using array_s = std::array<value_type, kDIM>;

    template <typename value_type>
    using vector_s = std::vector<value_type>;

    #endif


    #ifdef ALGEBRA_PLUGIN_INCLUDE_VC

    template <typename value_type, unsigned int kDIM>
    using array_v = simd::array<value_type, kDIM>;

    template <typename value_type>
    using vector_v = simd::aligned::vector<value_type>;

    #endif


    template <typename key_type, typename value_type>
    using map_s = std::map<key_type, value_type>;

    template< class... types>
    using tuple_s = std::tuple<types ...>;


} // namespace algebra

