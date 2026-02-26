# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_smatrix smatrix
   "include/algebra/smatrix.hpp"
)

target_link_libraries(
    algebra_smatrix
    INTERFACE
        algebra::common
        algebra::smatrix_storage
        algebra::smatrix_math
        algebra::generic_math
)

algebra_test_public_headers( algebra_smatrix
   "algebra/smatrix.hpp"
)
