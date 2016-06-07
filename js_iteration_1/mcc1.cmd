em++ ^
        -I C:\sohail\March\emscripten\boost_1_61_0\   ^
        -s EXPORTED_FUNCTIONS="['_main']" ^
        -s NO_EXIT_RUNTIME=0          ^
        -s ASSERTIONS=1               ^
        --profiling   ^
                                      ^
        -pedantic -std=c++14  ^
    mcc1.cpp  ^
        -o  mcc1.cpp.js
