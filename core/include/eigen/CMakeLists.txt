message(STATUS "Building 'algebra::eigen' plugin")

add_library(algebra_eigen INTERFACE)

target_include_directories(algebra_eigen
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${ALGEBRA_PLUGIN_SOURCE_DIR}/core/common/include/algebra)

install(
  DIRECTORY include/algebra
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    algebra_eigen
    INTERFACE -DALGEBRA_PLUGIN_CUSTOM_SCALARTYPE=${ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE})
endif()

find_package(Eigen3 REQUIRED)
target_link_libraries(algebra_eigen INTERFACE Eigen3::Eigen)

if(ALGEBRA_PLUGIN_USE_VECMEM)
  target_compile_options(algebra_eigen INTERFACE -DALGEBRA_PLUGIN_USE_VECMEM)
  target_link_libraries(algebra_eigen INTERFACE vecmem::core)
endif()

add_library(algebra::eigen ALIAS algebra_eigen)
