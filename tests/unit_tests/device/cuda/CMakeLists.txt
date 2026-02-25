# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable CUDA as a language.
enable_language(CUDA)
find_package(CUDAToolkit REQUIRED)

# Set the default CUDA compiler flags.
include(algebra-plugins-compiler-options-cuda)

# Macro setting up the CUDA tests for a specific algebra plugin.
macro(algebra_add_test_suite plugin)
    string(TOUPPER ${plugin} PLUGIN)

    algebra_add_test(${plugin}_cuda
        "algebra.cu"
        LINK_LIBRARIES GTest::gtest_main algebra::test_cuda algebra::${plugin}
    )

    target_compile_options(
        algebra_test_${plugin}_cuda
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
    )

    target_compile_definitions(
        algebra_test_${plugin}_cuda
        PUBLIC ALGEBRA_TEST_${PLUGIN}
    )
endmacro()

# Set up all of the (available) CUDA tests.
if(ALGEBRA_PLUGINS_INCLUDE_ARRAY)
    algebra_add_test_suite(array)
endif()

if(ALGEBRA_PLUGINS_INCLUDE_EIGEN)
    algebra_add_test_suite(eigen)
endif()
