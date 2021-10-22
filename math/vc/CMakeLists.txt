# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_vc_math vc_math
   "include/algebra/math/vc.hpp"
   "include/algebra/math/impl/vc_cartesian2.hpp"
   "include/algebra/math/impl/vc_cylindrical2.hpp"
   "include/algebra/math/impl/vc_polar2.hpp"
   "include/algebra/math/impl/vc_transform3.hpp"
   "include/algebra/math/impl/vc_vector.hpp" )
target_link_libraries( algebra_vc_math
   INTERFACE algebra::common algebra::cmath_math algebra::vc_storage Vc )
