#!/bin/bash
ls /lib/

cd /src/js_iteration_1




em++ \
    -I "/lib/boost_1_61_0"  \
    -I "/lib/eigen" \
    -O3   \
    -Oz \
    -s OUTLINING_LIMIT=100000 \
    --profiling     \
    -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
    -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size', '_get_pointset_ptr', '_get_pointset_size', '_about' ]" \
    -s NO_EXIT_RUNTIME=1          \
    -Winline         \
    -s TOTAL_MEMORY=30100100    \
    -s ABORTING_MALLOC=0 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s DISABLE_EXCEPTION_CATCHING=0  \
    -s DEMANGLE_SUPPORT=1 \
    -pedantic -std=c++14  \
mcc2.cpp  \
    -o  ../build/mcc2.compiled.js
exit 0


