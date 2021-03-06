# File: build_tests.sh
# Desc: Compile and run tests using GoogleTest and Emscripten.
# ------------------------------------------------
# This files accepts a file as input, let it be input_file.cpp and will compile
# it with emscripten to produce a JavaScript file with the name input_file.compiled.js
# and if the compilation is successfull it will also run the tests with node .
#
#
# Parameters:
#
# GTEST_ROOT:  root folder of your local googletest repo
#
# OPTIM:       switches between developemnt and optimized compilation #TODO optimized


if [ $# -eq 0 ]; then
    echo "Usage: bash build_tests.sh input_file.cpp."
    echo "You should provide an input file"
    exit 1
fi

OPTIM=1
GTEST_ROOT=/home/charles/MP5/googletest
fullfilename=$1 # get first argument as input file

extension="${fullfilename##*.}" # seperate extension and filename
filename="${fullfilename%.*}"   # this is needed for the output file name


# check if the variable EM_PREPARE is set to 1, which means that we can compile.
if [[ -n $EM_PREPARE ]]
then
    echo ""
else
    echo "You need to first run \"source em_prepare.sh\""
    exit 1;
fi


if [ $OPTIM -eq 0 ]
then
    cols=$( tput cols )
    tput bold
    echo " * * * Development Version * * * "
    echo Compiling ...
    em++    -I $BOOST_FOLDER -I $GTEST_ROOT/googletest/include \
            -s ASSERTIONS=1 \
            -pedantic -std=c++14 \
            "$1" \
            ${GTEST_ROOT}/build/googlemock/gtest/libgtest.a -o \
            "$filename".compiled.js 2>&1 \
            && js "$filename".compiled.js
    tput sgr0  # reset terminal options

fi

if [ $OPTIM -eq 1 ]
then
    cols=$( tput cols )
    tput bold
    echo " * * * Optimised Version * * * (doesnt work with      --llvm-lto 1  ) "
    echo Compiling ...
    em++    -I $BOOST_FOLDER -I $GTEST_ROOT/googletest/include \
            -s ASSERTIONS=1 \
            -pedantic -std=c++14 \
            -O3   \
            -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
            -s TOTAL_MEMORY=30100100   \
            "$1" \
            ${GTEST_ROOT}/build/googlemock/gtest/libgtest.a -o \
            "$filename".compiled.js 2>&1 \
            && js "$filename".compiled.js
    tput sgr0  # reset terminal options

fi

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
