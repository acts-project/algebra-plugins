# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable SYCL as a language.
vecmem_check_language( SYCL )

enable_language(SYCL)

# Set the default SYCL compiler flags.
include(algebra-plugins-compiler-options-sycl)

# Macro setting up the SYCL tests for a specific algebra plugin.
macro(algebra_add_test_suite plugin)
    string(TOUPPER ${plugin} PLUGIN)

    algebra_add_test(${plugin}_sycl
        "algebra.sycl"
        LINK_LIBRARIES GTest::gtest_main algebra::test_sycl algebra::${plugin}
    )

    target_compile_definitions(
        algebra_test_${plugin}_sycl
        PUBLIC ALGEBRA_TEST_${PLUGIN}
    )
endmacro()

# Set up all of the (available) SYCL tests.
if(ALGEBRA_PLUGINS_INCLUDE_ARRAY)
    algebra_add_test_suite(array)
endif()
