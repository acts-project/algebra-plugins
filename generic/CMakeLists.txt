# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the generic algebra library.
algebra_add_library(algebra_generic generic
   "include/algebra/generic.hpp"
   # impl include
   "include/algebra/impl/generic_matrix.hpp"
   "include/algebra/impl/generic_transform3.hpp"
   "include/algebra/impl/generic_vector.hpp"
   # algorithms include
   "include/algebra/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
   "include/algebra/algorithms/matrix/determinant/cofactor.hpp"
   "include/algebra/algorithms/matrix/determinant/hard_coded.hpp"
   "include/algebra/algorithms/matrix/determinant/partial_pivot_lud.hpp"
   "include/algebra/algorithms/matrix/inverse/cofactor.hpp"
   "include/algebra/algorithms/matrix/inverse/hard_coded.hpp"
   "include/algebra/algorithms/matrix/inverse/partial_pivot_lud.hpp"
   "include/algebra/algorithms/utils/algorithm_selector.hpp"
)
target_link_libraries(algebra_generic INTERFACE algebra::common algebra::utils)
algebra_test_public_headers( algebra_generic
   "algebra/generic.hpp"
)
