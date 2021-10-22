# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_vecmem_vc vecmem_vc
   "include/algebra/vecmem_vc.hpp" )
target_link_libraries( algebra_vecmem_vc
   INTERFACE algebra::common algebra::vecmem_storage algebra::cmath_math
             algebra::vc_math )
