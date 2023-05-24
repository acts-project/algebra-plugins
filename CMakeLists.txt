# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the project.
cmake_minimum_required( VERSION 3.14 )
project( algebra-plugins VERSION 0.20 LANGUAGES CXX )

# Set up the used C++ standard.
set( CMAKE_CXX_STANDARD 17 CACHE STRING "The (host) C++ standard to use" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "Disable (host) C++ extensions" )
set( CMAKE_CUDA_STANDARD 17 CACHE STRING "The (CUDA) C++ standard to use" )
set( CMAKE_CUDA_EXTENSIONS FALSE CACHE BOOL "Disable (CUDA) C++ extensions" )
set( CMAKE_SYCL_STANDARD 17 CACHE STRING "The (SYCL) C++ standard to use" )

# Standard CMake include(s).
include( CTest )
include( GNUInstallDirs )

# Explicitly set the output directory for the binaries. Such that if this
# project is included by another project, the main project's configuration would
# win out.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}" CACHE PATH
   "Directory for the built binaries" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
   "Directory for the built libraries" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
   "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
   "Directory for the built static libraries" )

# Flags controlling how the Algebra Plugins code should behave.
option( ALGEBRA_PLUGINS_INCLUDE_EIGEN
   "Include Eigen types in Algebra Plugins" FALSE )
option( ALGEBRA_PLUGINS_INCLUDE_SMATRIX
   "Include Smatrix types in Algebra Plugins" FALSE )
option( ALGEBRA_PLUGINS_INCLUDE_VC
   "Include Vc types in Algebra Plugins" FALSE )
option( ALGEBRA_PLUGINS_INCLUDE_FASTOR
   "Include Fastor types in Algebra Plugins" FALSE )
option( ALGEBRA_PLUGINS_INCLUDE_VECMEM
   "Include VecMem types in Algebra Plugins" FALSE )
option( ALGEBRA_PLUGINS_BUILD_TESTING "Build the unit tests of Algebra Plugins"
   TRUE )
option( ALGEBRA_PLUGINS_BUILD_BENCHMARKS "Build the benchmark suite of Algebra Plugins"
	FALSE )
option( ALGEBRA_PLUGINS_FAIL_ON_WARNINGS
   "Make the build fail on compiler/linker warnings" FALSE )

# Include the Algebra Plugins CMake code.
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" )
include( algebra-plugins-functions )

# Set up all "buildable" externals.
if( ALGEBRA_PLUGINS_BUILD_BENCHMARKS )
   add_subdirectory( extern/benchmark )
endif()
if( BUILD_TESTING AND ALGEBRA_PLUGINS_BUILD_TESTING )
   add_subdirectory( extern/googletest )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_EIGEN )
   add_subdirectory( extern/eigen3 )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_FASTOR )
   add_subdirectory( extern/fastor )
endif()
if( ALGEBRA_PLUGINS_INCLUDE_VC )
   add_subdirectory( extern/vc )
   if( NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )
      # Use the preferred compiler flags from Vc for the entire project. Do not
      # set them on the libraries individually, as the clients of the libraries
      # may want to use a different optimisation.
      vc_set_preferred_compiler_flags()
      add_compile_definitions( ${Vc_DEFINITIONS} )
      add_compile_options( ${Vc_COMPILE_FLAGS} ${Vc_ARCHITECTURE_FLAGS} )
   endif()
endif()
add_subdirectory( extern/vecmem ) # VecMem is always needed.

# Set up the Algebra Plugin libraries.
add_subdirectory( common )
add_subdirectory( storage )
add_subdirectory( math )
add_subdirectory( frontend )

# Set up the test(s).
if( BUILD_TESTING AND ALGEBRA_PLUGINS_BUILD_TESTING )
  add_subdirectory( tests )
endif()

# Set up the packaging of the project.
include( algebra-plugins-packaging )
