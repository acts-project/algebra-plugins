# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building the SMatrix algebra plugin")

find_package(ROOT COMPONENTS Core MathCore Smatrix REQUIRED)

# Set up the library.
algebra_add_library( algebra_smatrix smatrix
   "include/algebra/smatrix.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/impl/detail/*.hpp"
)

target_link_libraries(
    algebra_smatrix
    INTERFACE
        ROOT::Core
        ROOT::MathCore
        ROOT::Smatrix
        algebra::common
        algebra::generic
)

algebra_test_public_headers( algebra_smatrix
   "algebra/smatrix.hpp"
)
