# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_utils utils
   "include/algebra/utils/print.hpp" )
target_link_libraries( algebra_utils INTERFACE algebra::common )
algebra_test_public_headers( algebra_utils
   "algebra/utils/print.hpp" )
