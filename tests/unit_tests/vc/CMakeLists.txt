enable_testing()

set(PLUGIN_EXTRA_LIBRARIES Vc)

foreach(etest ${all_unit_tests})
    add_algebra_test(vc_array_algebra_${etest}
                     vc_array_algebra_${etest}.cpp 
                     algebra::vc_array)
endforeach(etest)

