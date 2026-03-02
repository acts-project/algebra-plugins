# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2022-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building the Fastor algebra plugin")

# Set up the library.
algebra_add_library( algebra_fastor fastor
   "include/algebra/fastor.hpp"
   "include/algebra/impl/*.hpp"
)

target_link_libraries(
    algebra_fastor
    INTERFACE Fastor::Fastor algebra::common algebra::generic_math
)

algebra_test_public_headers( algebra_fastor
   "algebra/fastor.hpp"
)
