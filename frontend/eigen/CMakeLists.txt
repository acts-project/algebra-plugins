# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_eigen eigen
   "include/algebra/eigen.hpp" )

target_link_libraries( algebra_eigen
   INTERFACE algebra::common algebra::eigen_storage algebra::eigen_math algebra::generic_math
             Eigen3::Eigen )

algebra_test_public_headers( algebra_eigen
   "algebra/eigen.hpp" )
