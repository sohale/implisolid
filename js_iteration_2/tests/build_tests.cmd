@echo off
set OPTIM=0

rem http://stackoverflow.com/questions/5683583/windows-batch-files-if-else

if %OPTIM% == 0 (

    @echo on
    echo "--------------------------------------------------- dev ------------------------------------------------------"

    em++ ^
            -I C:\sohail\March\emscripten\boost_1_61_0\   ^
              ^
            -s EXPORTED_FUNCTIONS="['_main' ]" ^
            -s NO_EXIT_RUNTIME=1          ^
            -Winline         ^
            -s TOTAL_MEMORY=300100    ^
            -s DEMANGLE_SUPPORT=1   ^
             -s ASSERTIONS=1               ^
                                        ^
            -pedantic -std=c++14  ^
        primitives_test.cpp  ^
            -o  ./build/primitives_test.compiled.js

    @rem         -c ^

    @rem         -s TOTAL_MEMORY=30100100    ^

    @rem -g   <------- preserve debug information

    @rem         --llvm-lto 1     ^
    @rem        -O3   ^
    @rem        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^

    @rem      -s ALLOW_MEMORY_GROWTH=1  ^   # This makes it 3X slower!
    @rem     -s TOTAL_MEMORY=16777216

) ELSE (
    @echo
    @echo "--------------------------------------------------- optim ------------------------------------------------------"
    @echo on
    em++ ^
            -I C:\sohail\March\emscripten\boost_1_61_0\   ^
              ^
            -s EXPORTED_FUNCTIONS="['_main' ]" ^
            -s NO_EXIT_RUNTIME=1          ^
            -Winline         ^
            -s TOTAL_MEMORY=30100100    ^
            -s DEMANGLE_SUPPORT=1   ^
             -s ASSERTIONS=1               ^
             --llvm-lto 1     ^
            -O3   ^
            -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
                                        ^
            -pedantic -std=c++14  ^
        primitives_test.cpp  ^
            -o  ./build/primitives_test.compiled.js



    @rem      -s ALLOW_MEMORY_GROWTH=1  ^   # This makes it 3X slower!
    @rem     -s TOTAL_MEMORY=16777216

)

