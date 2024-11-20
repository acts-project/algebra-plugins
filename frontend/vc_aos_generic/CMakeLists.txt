# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_vc_aos_generic vc_aos_generic
   "include/algebra/vc_aos_generic.hpp" )
target_link_libraries( algebra_vc_aos_generic
   INTERFACE algebra::common algebra::vc_aos_storage algebra::generic_math
             algebra::vc_aos_math algebra::cmath_math algebra::array_storage )
algebra_test_public_headers( algebra_vc_aos_generic
   "algebra/vc_aos_generic.hpp" )
