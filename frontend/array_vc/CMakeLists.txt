# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_array_vc array_vc
   "include/algebra/array_vc.hpp" )
target_link_libraries( algebra_array_vc
   INTERFACE algebra::common algebra::array_storage algebra::cmath_math
             algebra::vc_math )
