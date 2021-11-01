# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_cmath_math cmath_math
   "include/algebra/math/cmath.hpp"
   "include/algebra/math/impl/cmath_cartesian2.hpp"
   "include/algebra/math/impl/cmath_cylindrical2.hpp"
   "include/algebra/math/impl/cmath_getter.hpp"
   "include/algebra/math/impl/cmath_operator.hpp"
   "include/algebra/math/impl/cmath_polar2.hpp"
   "include/algebra/math/impl/cmath_transform3.hpp"
   "include/algebra/math/impl/cmath_vector.hpp" )
target_link_libraries( algebra_cmath_math
   INTERFACE algebra::common )
