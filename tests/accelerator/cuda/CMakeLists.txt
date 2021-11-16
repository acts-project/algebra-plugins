# Algebra plugins library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable CUDA as a language.
enable_language( CUDA )
find_package( CUDAToolkit REQUIRED )

# Test(s) for algebra::array_cmath with CUDA.
algebra_add_test( array_cuda
   "array_cmath_kernels.hpp" "array_cmath_kernels.cu" "array_cmath.cpp"
   "cuda_error_check.hpp" "cuda_error_check.cpp" "execute_cuda_test.cuh"
   LINK_LIBRARIES CUDA::cudart GTest::gtest_main vecmem::core vecmem::cuda
                  algebra::tests_common algebra::array_cmath )
target_compile_options( algebra_test_array_cuda
   PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr> )
