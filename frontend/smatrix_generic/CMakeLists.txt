# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Find Smatrix.
find_package( ROOT COMPONENTS Smatrix REQUIRED)

# Set up the library.
algebra_add_library( algebra_smatrix_generic smatrix_generic
   "include/algebra/smatrix_generic.hpp" )
target_link_libraries( algebra_smatrix_generic
   INTERFACE algebra::common algebra::smatrix_storage algebra::generic_math
             algebra::smatrix_math ROOT::Smatrix )
algebra_test_public_headers( algebra_smatrix_generic
   "algebra/smatrix_generic.hpp" )
