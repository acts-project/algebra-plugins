# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Run the CUDA test(s), if CUDA is available.
vecmem_check_language( CUDA )
if( CMAKE_CUDA_COMPILER )
   add_subdirectory( cuda )
endif()
