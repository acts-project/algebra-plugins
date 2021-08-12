cmake_minimum_required(VERSION 3.14)

project(algebra-plugins VERSION 0.1 LANGUAGES CXX)

# This is a submodule project
set(ALGEBRA_PLUGIN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

option(ALGEBRA_PLUGIN_UNIT_TESTS "Enable unit tests for algebra backends" On)
option(ALGEBRA_PLUGIN_BENCHMARKS "Enable benchmark tests for algebra bakends" On)

if(ALGEBRA_PLUGIN_INCLUDE_VC)
     find_package(Vc 1.4.1 REQUIRED)
     # Define an interface target for the Vc library
     add_library(Vc INTERFACE)
     target_include_directories(Vc INTERFACE ${Vc_INCLUDE_DIR})
     target_link_libraries(Vc INTERFACE libVc.a -L${Vc_LIB_DIR})
endif()

option(ALGEBRA_PLUGIN_BUILD_VC "Download and build local Vc" Off)

if(ALGEBRA_PLUGIN_BUILD_VC)
     add_subdirectory(extern/Vc)
endif()

# This is a submodule project
set(ALGEBRA_PLUGIN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
if (NOT CMAKE_INSTALL_INCLUDEDIR)
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
endif()
if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
endif()
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
endif()

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(core)
add_subdirectory(extern)
add_subdirectory(tests)
