# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Macro setting up the CPU tests for a specific algebra plugin.
macro(algebra_add_test_suite plugin)
    string(TOUPPER ${plugin} PLUGIN)

    algebra_add_test(${plugin}
        "algebra/algebra.cpp"
        LINK_LIBRARIES GTest::gtest_main algebra::test_cpu algebra::utils
                        algebra::${plugin}
    )

    target_compile_definitions(
        algebra_test_${plugin}
        PUBLIC ALGEBRA_TEST_${PLUGIN}
    )
endmacro()

# Set up all of the (available) "host" tests.
if(ALGEBRA_PLUGINS_INCLUDE_ARRAY)
    algebra_add_test_suite(array)
endif()

if(ALGEBRA_PLUGINS_INCLUDE_EIGEN)
    algebra_add_test_suite(eigen)
endif()

if(ALGEBRA_PLUGINS_INCLUDE_FASTOR)
    algebra_add_test_suite(fastor)
endif()

if(ALGEBRA_PLUGINS_INCLUDE_SMATRIX)
    algebra_add_test_suite(smatrix)
endif()

if(ALGEBRA_PLUGINS_INCLUDE_VC)
    algebra_add_test_suite(vc_aos)

    algebra_add_test(vc_soa
        "algebra/vc_soa.cpp"
        LINK_LIBRARIES GTest::gtest_main algebra::test_cpu algebra::utils
                        algebra::vc_soa
    )

    target_compile_definitions(algebra_test_vc_soa PUBLIC ALGEBRA_TEST_VC_SOA)
endif()
