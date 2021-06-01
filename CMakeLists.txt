cmake_minimum_required(VERSION 3.11)

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

# To find or download packages
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(ALGEBRA_PLUGIN_BUILD_VC)
     add_subdirectory(externals/Vc)
endif()

add_subdirectory(core)
add_subdirectory(tests)
