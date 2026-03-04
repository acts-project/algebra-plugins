# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "host" tests.
message(STATUS "Building algebra-plugins host unit tests")

add_subdirectory(cpu)

# Set up all of the "device" tests.
if(ALGEBRA_BUILD_CUDA OR ALGEBRA_BUILD_SYCL)
    add_subdirectory(device)
endif()
