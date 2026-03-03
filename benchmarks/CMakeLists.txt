# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2023-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the default C++ compiler flags.
include(algebra-plugins-compiler-options-cpp)

# Set up an INTERFACE library for the common header(s).
add_library(algebra_bench_common INTERFACE)

target_include_directories(
    algebra_bench_common
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/common/include"
)

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    target_compile_definitions(
        algebra_bench_common
        INTERFACE -D_USE_MATH_DEFINES
    )
endif()

add_library(algebra::bench_common ALIAS algebra_bench_common)

# Macro setting up the CPU benchmarks for a specific algebra plugin.
macro(algebra_add_benchmarks plugin)
    string(TOUPPER ${plugin} PLUGIN)

    add_library(algebra_bench_${plugin} INTERFACE)

    target_include_directories(
        algebra_bench_${plugin}
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${plugin}/include>"
    )

    target_link_libraries(
        algebra_bench_${plugin}
        INTERFACE algebra::${plugin} algebra::common
    )

    algebra_add_benchmark( ${plugin}_getter
    "common/src/benchmark_getter.cpp"
    LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                    algebra_bench_${plugin} algebra::${plugin}
    )
    target_compile_definitions(
        algebra_benchmark_${plugin}_getter
        PUBLIC ALGEBRA_BENCHMARK_${PLUGIN}
    )

    algebra_add_benchmark( ${plugin}_vector
    "common/src/benchmark_vector.cpp"
    LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                    algebra_bench_${plugin} algebra::${plugin}
    )
    target_compile_definitions(
        algebra_benchmark_${plugin}_vector
        PUBLIC ALGEBRA_BENCHMARK_${PLUGIN}
    )

    algebra_add_benchmark( ${plugin}_transform3
    "common/src/benchmark_transform3.cpp"
    LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                    algebra_bench_${plugin} algebra::${plugin}
    )
    target_compile_definitions(
        algebra_benchmark_${plugin}_transform3
        PUBLIC ALGEBRA_BENCHMARK_${PLUGIN}
    )

    algebra_add_benchmark( ${plugin}_matrix
    "common/src/benchmark_matrix.cpp"
    LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                    algebra_bench_${plugin} algebra::${plugin}
    )
    target_compile_definitions(
        algebra_benchmark_${plugin}_matrix
        PUBLIC ALGEBRA_BENCHMARK_${PLUGIN}
    )
endmacro()

# Set up all of the (available) "host" benchmarks.
if(ALGEBRA_PLUGINS_INCLUDE_ARRAY)
    algebra_add_benchmarks(array)
endif()

if(ALGEBRA_PLUGINS_INCLUDE_EIGEN)
    algebra_add_benchmarks(eigen)
endif()

if(ALGEBRA_PLUGINS_INCLUDE_FASTOR)
    algebra_add_benchmarks(fastor)
endif()

if(ALGEBRA_PLUGINS_INCLUDE_SMATRIX)
    algebra_add_benchmarks(smatrix)
endif()

if(ALGEBRA_PLUGINS_INCLUDE_VC)
    algebra_add_benchmarks(vc_aos)

    add_library(algebra_bench_vc_soa INTERFACE)
    target_include_directories(
        algebra_bench_vc_soa
        INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vc_soa/include>"
    )
    target_link_libraries(
        algebra_bench_vc_soa
        INTERFACE algebra::vc_soa algebra::common
    )
    algebra_add_benchmark(vc_soa_getter
      "vc_soa/vc_soa_getter.cpp"
       LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                      algebra_bench_vc_soa algebra::vc_soa
    )
    algebra_add_benchmark(vc_soa_vector
      "vc_soa/vc_soa_vector.cpp"
       LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                      algebra_bench_vc_soa algebra::vc_soa
    )
    algebra_add_benchmark(vc_soa_transform3
      "vc_soa/vc_soa_transform3.cpp"
       LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                      algebra_bench_vc_soa algebra::vc_soa
    )
    algebra_add_benchmark(vc_soa_matrix
      "vc_soa/vc_soa_matrix.cpp"
       LINK_LIBRARIES benchmark::benchmark algebra::bench_common
                      algebra_bench_vc_soa algebra::vc_soa
    )
endif()
