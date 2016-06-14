em++ ^
        -I C:\sohail\March\emscripten\boost_1_61_0\   ^
        -O3   ^
        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
 ^
        -s EXPORTED_FUNCTIONS="['_produce_object', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry']" ^
        -s NO_EXIT_RUNTIME=1          ^
        --llvm-lto 1     ^
        -Winline         ^
 ^
        -pedantic -std=c++14  ^
    mcc1.cpp  ^
        -o  mcc1.cpp.js

