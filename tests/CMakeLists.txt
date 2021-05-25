enable_testing()

add_subdirectory(common)

if(ALGEBRA_PLUGIN_UNIT_TESTS)
    add_subdirectory(unit_tests)
endif()

if(ALGEBRA_PLUGIN_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
