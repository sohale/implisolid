#!/bin/bash

cd $(git rev-parse --show-toplevel)
cd ./js_iteration_1/

# python ~/install/styleguide/cpplint/cpplint.py  --linelength=1200  centroids_projection.cpp 2>lint.txt

# echo "@rem em++ -I /usr/local/inlcude/boost_1_61_0/ -s EXPORTED_FUNCTIONS="['_main']" -pedantic -std=c++14  -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  -O3     vect_test_2.cpp    -o vect_test_2.js"


#BOOST_FOLDER=/usr/local/include/boost_1_61_0
#BOOST_FOLDER=./boost_1_61_0


#if [ -z ${EM_PREPARE+x} ]

if [[ -n "$EM_PREPARE" ]]
then
    echo ""
else
    echo "First run \"source em_prepare.sh\""
    exit
fi

OPTIM=0
SAVE_BC=0


TEMP=`getopt -o dob -l dev,opt,bitcode -- "$@"`
eval set -- "$TEMP"
# echo $TEMP
# extract options and their arguments into variables
USAGE=$'Usage: ./build_mcc2.sh [option]
        -d or --dev for development version
        -o or --opt for optimized version
        -b or --bitcode to emit llvm Code'

if [ "$#" -lt 2 ]; then
    echo "$USAGE"
    exit 1
fi

while true; do
    case "$1" in
        -d|--dev)
            OPTIM=0
            shift
            ;;
        -o|--opt)
            OPTIM=1
            shift
            ;;
        -b|--bitcode)
            SAVE_BC=1
            shift
            ;;
        --)
            shift; break ;;
        *) echo "Internal Error"; exit 1 ;;
    esac
done

if [ $SAVE_BC -eq 1 ]
then

    echo "Optimized version, bytecode will be emitted in the file " mcc2.ll
    g++ \
            --save-bc mcc2.bc \
            -I $BOOST_FOLDER   \
            -I $EIGEN_LIB_FOLDER \
            -O3   \
            --profiling     \
            -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
            -Winline         \
            -s TOTAL_MEMORY=30100100    \
            -s ABORTING_MALLOC=0 \
            -pedantic -std=c++14  \
        mcc2.cpp  \
            -o  mcc2.compiled.js
    echo "Producing human readable LLVM ... "
    llvm-dis mcc2.bc -o mcc2.ll && rm mcc2.bc # transform to human-readable form and delete bytecode file.
    exit 0
fi

if [ $OPTIM -eq 0 ]
then

    echo "** dev version **"

    EMCC_DEBUG=1

    g++ -I $BOOST_FOLDER  \
        -I $EIGEN_LIB_FOLDER \
        -pedantic -std=c++14 \
        mcc2.cpp -o mcc2.compiled.js
#        --profiling \
    exit 0
fi


if [ $OPTIM -eq 1 ]
then
    echo "** optimised version **"

    g++ \
        -I $BOOST_FOLDER   \
        -I $EIGEN_LIB_FOLDER \
        -O3   \
        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
        -Winline         \
        -pedantic -std=c++14  \
    mcc2.cpp  \
        -o  mcc2.compiled.js
    exit 0
fi

