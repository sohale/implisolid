# Compile tests using GoogleTest, with emscripten.
# ------------------------------------------------
# GTEST_ROOT:  root folder of your local googletest repo
# OPTIM:       switches between developemnt and optimized compilation #TODO optimized

OPTIM=0
GTEST_ROOT=~/googletest

# check if the variable EM_PREPARE is set to 1, which means that we can compile.
if [[ -n $EM_PREPARE ]]
then
    echo ""
else
    echo "You need to first run \"source em_prepare.sh\""
    exit 1
fi

if [ $OPTIM -eq 0 ]
then
    echo " * * * Development Version * * * "
    em++    -I $BOOST_FOLDER -I $GTEST_ROOT/googletest/include \
            -s ASSERTIONS=1 \
            -pedantic -std=c++14 \
            test_crisp_subtract.cpp \
            ${GTEST_ROOT}/build/googlemock/gtest/libgtest.a -o \
            test_crisp_subtract.js \

fi



#
# if [ $OPTIM -eq 1 ]
#     echo " * * * ERROR * * * "
#     # @echo on
#     # em++ ^
#     #         -I C:\sohail\March\emscripten\boost_1_61_0\   ^
#     #           ^
#     #         -s EXPORTED_FUNCTIONS="['_main' ]" ^
#     #         -s NO_EXIT_RUNTIME=1          ^
#     #         -Winline         ^
#     #         -s TOTAL_MEMORY=30100100    ^
#     #         -s DEMANGLE_SUPPORT=1   ^
#     #          -s ASSERTIONS=1               ^
#     #          --llvm-lto 1     ^
#     #         -O3   ^
#     #         -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
#     #                                     ^
#     #         -pedantic -std=c++14  ^
#     #     primitives_test.cpp  ^
#     #         -o  ./build/primitives_test.compiled.js
#     #
#     #
#     #
#     # @rem      -s ALLOW_MEMORY_GROWTH=1  ^   # This makes it 3X slower!
#     # @rem     -s TOTAL_MEMORY=16777216
# fi
