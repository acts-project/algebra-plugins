# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the default C++ compiler flags.
include( algebra-plugins-compiler-options-cpp )

# Set up an INTERFACE library for the common header(s).
add_library( algebra_bench_common INTERFACE )
target_include_directories( algebra_bench_common INTERFACE
   "${CMAKE_CURRENT_SOURCE_DIR}/common/include" )
if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )
   target_compile_definitions( algebra_bench_common INTERFACE
      -D_USE_MATH_DEFINES )
endif()
add_library( algebra::bench_common ALIAS algebra_bench_common )

# Set up all of the (available) "host" benchmarks.
add_library( algebra_bench_array INTERFACE )
target_include_directories( algebra_bench_array INTERFACE
   "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/array/include>" )
target_link_libraries(algebra_bench_array INTERFACE algebra::array_cmath
                                                    algebra::common_storage)
algebra_add_benchmark( array_getter
   "array/array_getter.cpp"
   LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                   algebra_bench_array algebra::array_cmath )
algebra_add_benchmark( array_vector
   "array/array_vector.cpp"
   LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                   algebra_bench_array algebra::array_cmath )
algebra_add_benchmark( array_transform3
   "array/array_transform3.cpp"
   LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                   algebra_bench_array algebra::array_cmath )
algebra_add_benchmark( array_matrix
   "array/array_matrix.cpp"
   LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                   algebra_bench_array algebra::array_cmath )

if( ALGEBRA_PLUGINS_INCLUDE_EIGEN )
   add_library( algebra_bench_eigen INTERFACE )
   target_include_directories( algebra_bench_eigen INTERFACE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/eigen/include>" )
   target_link_libraries(algebra_bench_eigen INTERFACE algebra::eigen_eigen
                                                    algebra::common_storage)
   algebra_add_benchmark( eigen_getter
      "eigen/eigen_getter.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_eigen algebra::eigen_eigen )
   algebra_add_benchmark( eigen_vector
      "eigen/eigen_vector.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_eigen algebra::eigen_eigen )
   algebra_add_benchmark( eigen_transform3
      "eigen/eigen_transform3.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_eigen algebra::eigen_eigen )
   algebra_add_benchmark( eigen_matrix
      "eigen/eigen_matrix.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_eigen algebra::eigen_eigen )
endif()

if( ALGEBRA_PLUGINS_INCLUDE_VC )
   add_library( algebra_bench_vc_aos INTERFACE )
   target_include_directories( algebra_bench_vc_aos INTERFACE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vc_aos/include>" )
   target_link_libraries(algebra_bench_vc_aos INTERFACE
                                                   algebra::vc_aos
                                                   algebra::common_storage)
   algebra_add_benchmark( vc_aos_getter
      "vc_aos/vc_aos_getter.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_vc_aos algebra::vc_aos )
   algebra_add_benchmark( vc_aos_vector
      "vc_aos/vc_aos_vector.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_vc_aos algebra::vc_aos )
   algebra_add_benchmark( vc_aos_transform3
      "vc_aos/vc_aos_transform3.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_vc_aos algebra::vc_aos )
   algebra_add_benchmark( vc_aos_matrix
      "vc_aos/vc_aos_matrix.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_vc_aos algebra::vc_aos )

   add_library( algebra_bench_vc_soa INTERFACE )
   target_include_directories( algebra_bench_vc_soa INTERFACE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vc_soa/include>" )
   target_link_libraries(algebra_bench_vc_soa INTERFACE
                                                   algebra::vc_soa
                                                   algebra::common_storage)
   algebra_add_benchmark( vc_soa_getter
      "vc_soa/vc_soa_getter.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_vc_soa algebra::vc_soa )
   algebra_add_benchmark( vc_soa_vector
      "vc_soa/vc_soa_vector.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_vc_soa algebra::vc_soa )
   algebra_add_benchmark( vc_soa_transform3
      "vc_soa/vc_soa_transform3.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_vc_soa algebra::vc_soa )
   algebra_add_benchmark( vc_soa_matrix
      "vc_soa/vc_soa_matrix.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_vc_soa algebra::vc_soa )
endif()

if( ALGEBRA_PLUGINS_INCLUDE_FASTOR )
   add_library( algebra_bench_fastor INTERFACE )
   message(STATUS "Include directory: ${CMAKE_CURRENT_SOURCE_DIR}/fastor/include")
   target_include_directories( algebra_bench_fastor INTERFACE
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fastor/include>" )
   target_link_libraries(algebra_bench_fastor INTERFACE algebra::fastor_fastor
                                                    algebra::common_storage)
   algebra_add_benchmark( fastor_getter
      "fastor/fastor_getter.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_fastor algebra::fastor_fastor )
   algebra_add_benchmark( fastor_vector
      "fastor/fastor_vector.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_fastor algebra::fastor_fastor )
   algebra_add_benchmark( fastor_transform3
      "fastor/fastor_transform3.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_fastor algebra::fastor_fastor )
   algebra_add_benchmark( fastor_matrix
      "fastor/fastor_matrix.cpp"
      LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                     algebra_bench_fastor algebra::fastor_fastor )
endif()
