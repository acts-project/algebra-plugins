# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the test(s).
algebra_add_test( vecmem_array "vecmem_array_algebra_plugin.cpp"
   LINK_LIBRARIES GTest::gtest_main algebra::tests_common
                  algebra::vecmem_storage algebra::cmath_math )
