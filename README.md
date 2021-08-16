# algebra-plugins

This repository provides different algebra plugins with minimal functionality for the R&D projects `detray` and `traccc`.

To build it standalone, run e.g.

mkdir build
cd build

cmake -DALGEBRA_PLUGIN_UNIT_TESTS=ON -DALGEBRA_PLUGIN_BENCHMARKS=ON -DALGEBRA_PLUGIN_BUILD_VC=ON -DALGEBRA_PLUGIN_BUILD_GOOGLE_BENCHMARK=ON -DALGEBRA_PLUGIN_INCLUDE_ARRAY=ON ..

make -j5

ALGEBRA_PLUGIN_INCLUDE_XXX: Build XXX based plugin (array, eigen, smatric, vc)
ALGEBRA_PLUGIN_USE_VECMEM: Use vecmem container types

