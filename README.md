# algebra-plugins

This repository provides different algebra plugins with minimal functionality for the R&D projects `detray` and `traccc`.

To build it standalone, run e.g.

mkdir build
cd build

cmake -DALGEBRA_PLUGIN_INCLUDE_ARRAY=ON ..

make -j5

Available options:
ALGEBRA_PLUGIN_INCLUDE_XXX: Build XXX based plugin (e.g. array, eigen)
ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE: Choose a value type (e.g. float, double)
ALGEBRA_PLUGIN_UNIT_TESTS: Build unit tests (needs googletest)
ALGEBRA_PLUGIN_BUILD_GOOGLE_BENCHMARK: Build google benchmark if not provided
ALGEBRA_PLUGIN_USE_VECMEM: Use vecmem container types
ALGEBRA_PLUGIN_BUILD_VECMEM: Build vecmem if not provided (automatically set if vecmem cannot be found)
ALGEBRA_PLUGIN_BUILD_VC: Build Vc if not provided

