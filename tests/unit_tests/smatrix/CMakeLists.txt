enable_testing()

set(PLUGIN_EXTRA_LIBRARIES ROOT::MathCore)

foreach(etest ${all_unit_tests})
    add_algebra_test(smatrix_algebra_${etest}
                     smatrix_algebra_${etest}.cpp 
                     algebra::smatrix)
endforeach(etest)
