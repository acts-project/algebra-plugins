cmake_minimum_required(VERSION 3.14)

project(algebra-plugins VERSION 0.1 LANGUAGES CXX)

option(ALGEBRA_PLUGIN_BUILD_VC "Download and build local Vc" Off)

option(ALGEBRA_PLUGIN_UNIT_TESTS "Enable unit tests for algebra backends" On)
option(ALGEBRA_PLUGIN_BENCHMARKS "Enable benchmark tests for algebra bakends" Off)

if (NOT EIGEN3_INCLUDE_DIRS)
    find_package(Eigen3 REQUIRED)
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})
endif()

# This is a submodule project
set(ALGEBRA_PLUGIN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
if (NOT CMAKE_INSTALL_INCLUDEDIR)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
endif()
if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
endif()
if (NOT CMAKE_INSTALL_LIBDIR)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
endif()

add_subdirectory(core)
add_subdirectory(extern)
add_subdirectory(tests)
