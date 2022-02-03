# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# FindCUDAToolkit needs at least CMake 3.17.
cmake_minimum_required( VERSION 3.17 )

# Include the helper function(s).
include( algebra-plugins-functions )

# Figure out the properties of CUDA being used.
find_package( CUDAToolkit REQUIRED )

# Set the architecture to build code for.
set( CMAKE_CUDA_ARCHITECTURES "52" CACHE STRING
   "CUDA architectures to build device code for" )

# Make CUDA generate debug symbols for the device code as well in a debug
# build.
algebra_add_flag( CMAKE_CUDA_FLAGS_DEBUG "-G" )

# More rigorous tests for the Debug builds.
if( "${CUDAToolkit_VERSION}" VERSION_GREATER_EQUAL "10.2" )
   algebra_add_flag( CMAKE_CUDA_FLAGS_DEBUG "-Werror all-warnings" )
endif()
