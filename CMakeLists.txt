cmake_minimum_required(VERSION 3.11)

project(algebra-plugins VERSION 0.1 LANGUAGES CXX)

option(ALGEBRA_PLUGIN_BUILD_VC "Download and build local Vc" Off)

option(ALGEBRA_PLUGIN_UNIT_TESTS "Enable unit tests for algebra backends" On)
option(ALGEBRA_PLUGIN_BENCHMARKS "Enable benchmark tests for algebra bakends" On)

# To find or download packages
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# This is a submodule project
set(ALGEBRA_PLUGIN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#if(not )

add_subdirectory(core)
add_subdirectory(external)
