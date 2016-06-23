em++ ^
        -I C:\sohail\March\emscripten\boost_1_61_0\   ^
        -s EXPORTED_FUNCTIONS="['_main']" ^
        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
        -pedantic -std=c++14  ^
    timer_test.cpp  ^
        -o timer_test.cpp.js

::     -s EXPORTED_FUNCTIONS="['_main', '_make_object']" ^
::     -s EXPORTED_FUNCTIONS="['_main']" ^
:: @rem     -O3   ^
:: @rem     -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^


::     -O3   ^
::     -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
