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

namespace algebra
{
    template <typename value_type, unsigned int kDIM>
    using array_t = std::array<value_type, kDIM>;

    template <typename value_type>
    using vector_t = std::vector<value_type>;

    template <typename key_type, typename value_type>
    using map_t = std::map<key_type, value_type>;

    template< class... types>
    using tuple_t = std::tuple<types ...>;


} // namespace algebra

