
rem set EMCC_DEBUG=1

em++ ^
        -I C:\sohail\March\emscripten\boost_1_61_0\   ^
        -s ASSERTIONS=1               ^
        -s DEMANGLE_SUPPORT=1   ^
                                    ^
        -s EXPORTED_FUNCTIONS="['_produce_object', '_main']" ^
        -s NO_EXIT_RUNTIME=1          ^
                                          ^
        -pedantic -std=c++14  ^
    mcc1.cpp  ^
        -o  mcc1.cpp.js




rem        -s EXPORTED_FUNCTIONS="['_main']" ^
rem        -s NO_EXIT_RUNTIME=0          ^
rem      --profiling   ^
