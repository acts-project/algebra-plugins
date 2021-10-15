# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up an INTERFACE library for the common header(s).
add_library( algebra_tests_common INTERFACE )
target_include_directories( algebra_tests_common INTERFACE
   "${CMAKE_CURRENT_SOURCE_DIR}" )
target_link_libraries( algebra_tests_common INTERFACE
   algebra::common )
add_library( algebra::tests_common ALIAS algebra_tests_common )
