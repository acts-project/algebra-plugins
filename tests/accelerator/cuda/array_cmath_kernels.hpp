/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/array_cmath.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

/// Function performing vector operation tests
void array_cmath_2d_vector_ops(
    vecmem::data::vector_view<const algebra::array::point2> a,
    vecmem::data::vector_view<const algebra::array::point2> b,
    vecmem::data::vector_view<algebra::scalar> output_host,
    vecmem::data::vector_view<algebra::scalar> output_device);
