# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_utils utils
   "include/algebra/utils/approximately_equal.hpp"
   "include/algebra/utils/casts.hpp"
   "include/algebra/utils/print.hpp" )
target_link_libraries( algebra_utils INTERFACE algebra::common algebra_common_math )
algebra_test_public_headers( algebra_utils
   "algebra/utils/approximately_equal.hpp"
   "algebra/utils/casts.hpp"
   "algebra/utils/print.hpp" )
