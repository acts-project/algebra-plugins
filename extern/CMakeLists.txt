
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
     message(VERBOSE "Need google test/benchmark to be intalled for unittests/benchmarks to work")
endif()
