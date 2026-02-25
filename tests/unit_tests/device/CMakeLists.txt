# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the CUDA tests.
if(ALGEBRA_BUILD_CUDA)
    add_subdirectory(cuda)
endif()

# Set up all of the CUDA tests.
if(ALGEBRA_BUILD_SYCL)
    add_subdirectory(sycl)
endif()
