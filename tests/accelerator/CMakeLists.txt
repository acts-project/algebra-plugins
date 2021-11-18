# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Run the CUDA test(s) by default, if CUDA is available.
vecmem_check_language( CUDA )
set( _cudaEnabledDefault FALSE )
if( CMAKE_CUDA_COMPILER )
   set( _cudaEnabledDefault TRUE )
endif()
option( ALGEBRA_PLUGINS_TEST_CUDA "Test the code in CUDA device code"
   ${_cudaEnabledDefault} )
unset( _cudaEnabledDefault )

if( ALGEBRA_PLUGINS_TEST_CUDA )
   add_subdirectory( cuda )
endif()
