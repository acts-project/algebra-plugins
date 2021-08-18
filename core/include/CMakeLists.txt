

if (ALGEBRA_PLUGIN_INCLUDE_ARRAY)
    add_subdirectory(array)
endif()

if (ALGEBRA_PLUGIN_INCLUDE_EIGEN)
    add_subdirectory(eigen)
endif()

if (ALGEBRA_PLUGIN_INCLUDE_SMATRIX)
    add_subdirectory(smatrix)
endif()

if (ALGEBRA_PLUGIN_INCLUDE_VC)
    add_subdirectory(vc)
endif()

if(ALGEBRA_PLUGIN_USE_VECMEM)
    add_subdirectory(vecmem)
endif()