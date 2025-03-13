/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

 #pragma once

 // Project include(s)
 #include "algebra/fastor_fastor.hpp"

namespace algebra {

/// Fill an @c Fastor based vector with random values
template <concepts::vector vector_t>
inline void fill_random_vec(std::vector<vector_t> &collection) {

auto rand_obj = [](vector_t& v) { return v.random(); };

collection.resize(collection.capacity());
std::ranges::for_each(collection, rand_obj);
}

/// Fill a @c Fastor based matrix with random values
template <concepts::matrix matrix_t>
inline void fill_random_matrix(std::vector<matrix_t> &collection) {

auto rand_obj = [](matrix_t& m) { return m.random(); };

collection.resize(collection.capacity());
std::ranges::for_each(collection, rand_obj);
}

}
