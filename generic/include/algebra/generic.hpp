/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Impl include(s).
#include "algebra/boolean.hpp"
#include "algebra/impl/generic_matrix.hpp"
#include "algebra/impl/generic_transform3.hpp"
#include "algebra/impl/generic_vector.hpp"
#include "algebra/math.hpp"

// Algorithms include(s).
#include "algebra/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
#include "algebra/algorithms/matrix/determinant/cofactor.hpp"
#include "algebra/algorithms/matrix/determinant/hard_coded.hpp"
#include "algebra/algorithms/matrix/determinant/partial_pivot_lud.hpp"
#include "algebra/algorithms/matrix/inverse/cofactor.hpp"
#include "algebra/algorithms/matrix/inverse/hard_coded.hpp"
#include "algebra/algorithms/matrix/inverse/partial_pivot_lud.hpp"
#include "algebra/algorithms/utils/algorithm_selector.hpp"
