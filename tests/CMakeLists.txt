# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the default C++ compiler flags.
include(algebra-plugins-compiler-options-cpp)

# Set up an INTERFACE library for the common header(s).
add_library(algebra_tests_common INTERFACE)
target_include_directories(
    algebra_tests_common
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/common"
)
target_link_libraries(algebra_tests_common INTERFACE algebra::utils)
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    target_compile_definitions(
        algebra_tests_common
        INTERFACE -D_USE_MATH_DEFINES
    )
endif()
add_library(algebra::tests_common ALIAS algebra_tests_common)

# Set up all of the (available) "host" tests.
algebra_add_test( array
   "array/array.cpp"
   LINK_LIBRARIES GTest::gtest_main algebra::tests_common
                  algebra::array
)

if(ALGEBRA_PLUGINS_INCLUDE_EIGEN)
    algebra_add_test(eigen
      "eigen/eigen.cpp"
       LINK_LIBRARIES GTest::gtest_main algebra::tests_common
                      algebra::eigen
    )
endif()

if(ALGEBRA_PLUGINS_INCLUDE_SMATRIX)
    algebra_add_test(smatrix
      "smatrix/smatrix.cpp" "smatrix/smatrix.cpp"
       LINK_LIBRARIES GTest::gtest_main algebra::tests_common
                      algebra::smatrix algebra::smatrix
    )
endif()

if(ALGEBRA_PLUGINS_INCLUDE_VC)
    algebra_add_test(vc_aos
      "vc_aos/vc_aos.cpp"
       LINK_LIBRARIES GTest::gtest_main algebra::tests_common
                      algebra::vc_aos
    )

    algebra_add_test(vc_soa
      "vc_soa/vc_soa.cpp"
       LINK_LIBRARIES GTest::gtest_main algebra::utils algebra::vc_soa
    )
endif()

if(ALGEBRA_PLUGINS_INCLUDE_FASTOR)
    algebra_add_test(fastor
      "fastor/fastor.cpp"
       LINK_LIBRARIES GTest::gtest_main algebra::tests_common
                      algebra::fastor
    )
endif()

# Set up all of the (available) "device" tests.
add_subdirectory(accelerator)
