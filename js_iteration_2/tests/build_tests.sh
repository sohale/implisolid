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
GTEST_ROOT=~/googletest
fullfilename=$1 # get first argument as input file
BOOST_FOLDER=/usr/local/include/boost_1_61_0/
extension="${fullfilename##*.}" # seperate extension and filename
filename="${fullfilename%.*}"   # this is needed for the output file name


# check if the variable EM_PREPARE is set to 1, which means that we can compile.
if [[ -n $EM_PREPARE ]]
then
    echo ""
else
    echo "Error: You need to first run \"source em_prepare.sh\""
    echo
    exit 1;
fi

if [ ! -d "$GTEST_ROOT" ]; then
    echo "Error: googletest directory not found at $GTEST_ROOT"
    echo
    exit 1;
fi

if [ ! -d "$BOOST_FOLDER" ]; then
    echo "Error: Boost 1.61.0+ not found at $BOOST_FOLDER"
    echo
    exit 1;
fi


if [ $OPTIM -eq 0 ]
then
    cols=$( tput cols )
    tput bold
    echo " * * * Development Version * * * "
    echo Compiling ...
    echo Error messages, and output information will be sent to "em_compile.log"
    em++    -I $BOOST_FOLDER -I $GTEST_ROOT/googletest/include \
            -s ASSERTIONS=1 \
            -pedantic -std=c++14 \
            "$1" \
            ${GTEST_ROOT}/build/googlemock/gtest/libgtest.a -o \
            "$filename".compiled.js > em_compile.log 2>&1 \
            && echo "compile success." \
            && node "$filename".compiled.js \
    tput sgr0  # reset terminal options

fi


if [ $OPTIM -eq 1 ]; then
    echo " * * * Optimized Version  * * *  "
    em++    -I $BOOST_FOLDER -s EXPORTED_FUNCTIONS="['_main' ]"  \
            -s NO_EXIT_RUNTIME=1                \
            -Winline    \
            -s TOTAL_MEMORY=301001000           \
            -s DEMANGLE_SUPPORT=1               \
             -s ASSERTIONS=1                    \
             --llvm-lto 1                       \
            -O3                                 \
            -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS \
            -pedantic -std=c++14  \
             "$filename".cpp      \
            -o  "$filename".compiled.js
fi


