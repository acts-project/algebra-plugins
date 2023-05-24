# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2022-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.24 )
include( FetchContent )

# Declare where to get GoogleTest from, in case it is not found with
# find_package.
set( ALGEBRA_PLUGINS_BENCHMARK_SOURCE
   "URL;https://github.com/google/benchmark/archive/refs/tags/v1.6.1.tar.gz;URL_MD5;8c33c51f9b7154e6c290df3750081c87"
   CACHE STRING "Source for Google Benchmark, when built as part of this project" )
mark_as_advanced( ALGEBRA_PLUGINS_BENCHMARK_SOURCE )
FetchContent_Declare( benchmark ${ALGEBRA_PLUGINS_BENCHMARK_SOURCE}
   FIND_PACKAGE_ARGS )

# Options used in the build of Google Benchmark.
set( BENCHMARK_ENABLE_TESTING OFF CACHE BOOL
   "Don't build the Google Benchmark tests" )
set( BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL
   "Turn off the installation of Google Benchmark" )
set( BENCHMARK_INSTALL_DOCS OFF CACHE BOOL
   "Don't install the Google Benchmark docs" )
set( BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL
   "Don't build the unit tests which depend on gtest" )
set( BENCHMARK_ENABLE_WERROR OFF CACHE BOOL
   "Enable/disable using -Werror in the build of Google Benchmark" )

# Set it up.
FetchContent_MakeAvailable( benchmark )
