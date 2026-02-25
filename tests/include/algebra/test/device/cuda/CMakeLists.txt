# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building 'algebra::test_cuda' component")

# CMake version requirement.
cmake_minimum_required(VERSION 3.18)

# Enable CUDA as a language.
enable_language(CUDA)
find_package(CUDAToolkit REQUIRED)

# Set the default CUDA compiler flags.
include(algebra-plugins-compiler-options-cuda)

if(${CMAKE_CUDA_STANDARD} LESS 20)
    message(
        SEND_ERROR
        "CMAKE_CUDA_STANDARD=${CMAKE_CUDA_STANDARD}, but algebra-plugins requires C++>=20"
    )
endif()

# Helper library for the CUDA tests.
add_library(
    algebra_test_cuda
    STATIC
    "cuda_error_check.hpp"
    "cuda_error_check.cpp"
    "cuda_test_fixture.hpp"
    "execute_cuda_test.cuh"
    "algebra_test_suite.cuh"
)
target_include_directories(
    algebra_test_cuda
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(
    algebra_test_cuda
    PUBLIC CUDA::cudart
    INTERFACE
        GTest::gtest
        algebra::test_framework
        algebra::test_device
        vecmem::core
        vecmem::cuda
)

target_compile_options(
    algebra_test_cuda
    INTERFACE
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe
        --diag_suppress=177>
        $<$<COMPILE_LANGUAGE:CUDA>:--ftz=false
        --prec-div=true>
        $<$<COMPILE_LANGUAGE:CUDA>:--prec-sqrt=true
        --fmad=false>
)

if("${CUDAToolkit_VERSION}" VERSION_GREATER_EQUAL "11.5")
    # Replace all "#pragma diag_suppress" calls in the Eigen code with
    # "#pragma nv_diag_suppress" calls. To be removed once the Eigen
    # code starts behaving better. Note that this is only needed for
    # Windows, as these warnings are suppressed on Linux by declaring
    # Eigen as a "system include".
    target_compile_definitions(
        algebra_test_cuda
        INTERFACE diag_suppress=nv_diag_suppress
    )
endif()

add_library(algebra::test_cuda ALIAS algebra_test_cuda)
