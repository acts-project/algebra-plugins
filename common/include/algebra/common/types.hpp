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

#ifdef ALGEBRA_PLUGIN_INCLUDE_VC
#include "simd_types.hpp"
#endif

namespace algebra
{
    template <typename value_type, unsigned int kDIM>
    using array_s = std::array<value_type, kDIM>;

    template <typename value_type>
    using vector_s = std::vector<value_type>;

    template <typename key_type, typename value_type>
    using map_s = std::map<key_type, value_type>;

    template< class... types>
    using tuple_s = std::tuple<types ...>;


    #ifdef ALGEBRA_PLUGIN_INCLUDE_VC

    template <typename value_type, unsigned int kDIM>
    using array_v = simd::array<value_type, kDIM>;

    template <typename value_type>
    using vector_v = simd::aligned::vector<value_type>;

    #endif


} // namespace algebra

