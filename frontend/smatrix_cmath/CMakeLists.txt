# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Find Smatrix.
find_package( ROOT COMPONENTS Smatrix REQUIRED)

# Set up the library.
algebra_add_library( algebra_smatrix_cmath smatrix_cmath
   "include/algebra/smatrix_cmath.hpp" )
target_link_libraries( algebra_smatrix_cmath
   INTERFACE algebra::common algebra::smatrix_storage algebra::cmath_math
             algebra::smatrix_math ROOT::Smatrix )
