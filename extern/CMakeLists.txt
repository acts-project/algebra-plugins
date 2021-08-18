
if(ALGEBRA_PLUGIN_BUILD_GOOGLE_BENCHMARK)
     message(STATUS "Building googletest")

     include(FetchContent)
     FetchContent_Declare(
       googletest
       GIT_REPOSITORY https://github.com/google/googletest.git
       GIT_TAG        release-1.8.0
     )

     FetchContent_MakeAvailable(googletest)

     if (ALGEBRA_PLUGIN_BENCHMARKS)
     # Also need benchmark
     endif()
elseif(ALGEBRA_PLUGIN_UNIT_TESTS OR ALGEBRA_PLUGIN_BENCHMARKS)
     message(WARNING "Need google test/benchmark to be installed for unit tests/benchmarks to work")
endif()

if(ALGEBRA_PLUGIN_BUILD_VECMEM)
     message(STATUS "Building vecmem")
     include( FetchContent )
     FetchContent_Declare( 
          Vecmem
          GIT_REPOSITORY "https://github.com/acts-project/vecmem.git"
          GIT_TAG "v0.2.0"
     )

     # Do not build google test in vecmem again
     set(BUILD_TESTING OFF CACHE INTERNAL "")

     # Get it into the current directory.
     FetchContent_Populate( vecmem )
     add_subdirectory( "${vecmem_SOURCE_DIR}" "${vecmem_BINARY_DIR}"
     EXCLUDE_FROM_ALL )
endif()
