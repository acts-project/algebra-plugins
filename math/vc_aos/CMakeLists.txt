# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_vc_aos_math vc_aos_math
   "include/algebra/math/vc_aos.hpp"
   "include/algebra/math/impl/vc_aos_matrix.hpp"
   "include/algebra/math/impl/vc_aos_transform3.hpp"
   "include/algebra/math/impl/vc_aos_vector.hpp" )
target_link_libraries( algebra_vc_aos_math
   INTERFACE Vc::Vc algebra::common algebra::common_math algebra::generic_math algebra::vc_aos_storage )
algebra_test_public_headers( algebra_vc_aos_math
   "algebra/math/vc_aos.hpp" )
