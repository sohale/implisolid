em++ ^
        -I C:\sohail\March\emscripten\boost_1_61_0\   ^
        -O3   ^
        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
 ^
        -s EXPORTED_FUNCTIONS="['_produce_object', '_main']" ^
        -s NO_EXIT_RUNTIME=1          ^
 ^
        -pedantic -std=c++14  ^
    mcc1.cpp  ^
        -o  mcc1.cpp.js

