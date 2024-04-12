# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the library.
algebra_add_library( algebra_eigen_generic eigen_generic
   "include/algebra/eigen_generic.hpp" )
target_link_libraries( algebra_eigen_generic
   INTERFACE algebra::common algebra::eigen_storage algebra::generic_math
             algebra::eigen_math Eigen3::Eigen )
algebra_test_public_headers( algebra_eigen_generic
   "algebra/eigen_generic.hpp" )
