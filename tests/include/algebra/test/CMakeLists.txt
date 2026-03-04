# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

add_subdirectory(cpu)

if(ALGEBRA_BUILD_CUDA OR ALGEBRA_BUILD_SYCL)
    add_subdirectory(device)
endif()
