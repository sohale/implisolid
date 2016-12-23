#!/bin/bash

cd $(git rev-parse --show-toplevel)
cd ./js_iteration_1/
mkdir -p ../build

# python ~/install/styleguide/cpplint/cpplint.py  --linelength=1200  centroids_projection.cpp 2>lint.txt

# echo "@rem em++ -I /usr/local/inlcude/boost_1_61_0/ -s EXPORTED_FUNCTIONS="['_main']" -pedantic -std=c++14  -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  -O3     vect_test_2.cpp    -o vect_test_2.js"


#BOOST_FOLDER=/usr/local/include/boost_1_61_0
#BOOST_FOLDER=./boost_1_61_0


#if [ -z ${EM_PREPARE+x} ]

if [[ -n "$EM_PREPARE" ]]
then
    echo -n ""
else
    echo "First run \"source em_prepare.sh\""
    exit
fi

OPTIM=0
SAVE_BC=0
WORKER=0

TEMP=`getopt -o dobw -l dev,opt,bitcode,worker -- "$@"`
eval set -- "$TEMP"
# echo $TEMP
# extract options and their arguments into variables
USAGE=$'Usage: ./build_mcc2.sh [option]
        -d or --dev for development version
        -o or --opt for optimized version
        -b or --bitcode to emit llvm Code
        -w or --worker for worker development'

if [ "$#" -lt 2 ]; then
    echo "$USAGE"
    exit 1
fi

while true; do
    case "$1" in
        -d|--dev)
            OPTIM=0
            WORKER=0
            shift
            ;;
        -o|--opt)
            OPTIM=1
            WORKER=0
            shift
            ;;
        -b|--bitcode)
            SAVE_BC=1
            WORKER=0
            shift
            ;;
        -w|--worker)
            OPTIM=0
            SAVE_BC=0
            WORKER=1
            shift
            ;;
        --)
            shift; break ;;
        *) echo "Internal Error"; exit 1 ;;
    esac
done

if [ $SAVE_BC -eq 1 ]
then

    echo "BitCode version, bitecode will be emitted in the file. Optimisation is ON. " mcc2.ll
    em++ \
            --save-bc mcc2.bc \
            -I $BOOST_FOLDER   \
            -I $EIGEN_LIB_FOLDER \
            -O3   \
            --profiling     \
            -s DEMANGLE_SUPPORT=1  \
            -s ASSERTIONS=1  \
            -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
            -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size', '_get_pointset_ptr', '_get_pointset_size', _about' ]" \
            -s NO_EXIT_RUNTIME=1          \
            -Winline         \
            -s TOTAL_MEMORY=30100100    \
            -s ABORTING_MALLOC=0 \
            -DMORE_ABOUT_INFO='"BitCode"'   \
            -Wall  \
            -pedantic -std=c++14  \
        mcc2.cpp  \
            -o  ../build/mcc2.compiled.js
    echo "Producing human readable LLVM ... "
    llvm-dis mcc2.bc -o mcc2.ll && rm mcc2.bc # transform to human-readable form and delete bytecode file.
    exit 0
fi

if [ $WORKER -eq 0 ]
then

if [ $OPTIM -eq 0 ]
then

    echo "** dev version (non-worker): Optimisations are off. Asserts are on. **"

    EMCC_DEBUG=1

    em++ -I $BOOST_FOLDER  \
         -I $EIGEN_LIB_FOLDER \
        -s TOTAL_MEMORY=30100100 \
        -s ABORTING_MALLOC=0 \
        -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size', '_get_pointset_ptr', '_get_pointset_size', '_about' ]" \
        -s NO_EXIT_RUNTIME=1 \
        -s DEMANGLE_SUPPORT=1 \
        -s ASSERTIONS=1 \
        -s BUILD_AS_WORKER=0 -DNOT_BUILD_AS_WORKER  \
        -DMORE_ABOUT_INFO='"DEV"'   \
        -Wall  \
        -pedantic -std=c++14 \
        mcc2.cpp -o ../build/mcc2.compiled.js
#        --profiling \
    exit 0
fi

fi

if [ $WORKER -eq 1 ]
then

if [ $OPTIM -eq 0 ]
then

    echo "** Worker version (dev) **"

    EMCC_DEBUG=1

    em++ -I $BOOST_FOLDER  \
         -I $EIGEN_LIB_FOLDER \
        -s TOTAL_MEMORY=30100100 \
        -s ABORTING_MALLOC=0 \
        -s EXPORTED_FUNCTIONS="['_main', '_build_geometry' , '_about' , '_worker_api__started', '_worker_api__verts', '_workerapi_make_geometry_part1' ]" \
        -s NO_EXIT_RUNTIME=1 \
        -s DEMANGLE_SUPPORT=1 \
        -s ASSERTIONS=1 \
        -s BUILD_AS_WORKER=1 -DBUILD_AS_WORKER  \
        -DMORE_ABOUT_INFO='"Worker+dev"'   \
        -Wall  \
        -pedantic -std=c++14 \
        mcc2.cpp -o ../build/mcc2.compiled.js
#        --profiling \
    exit 0
fi

fi

if [ $OPTIM -eq 1 ]
then
    echo "** optimised version. Asserts are off. **"

    em++ \
        -I $BOOST_FOLDER   \
        -I $EIGEN_LIB_FOLDER \
        -O3   \
        --profiling     \
        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
        -DMORE_ABOUT_INFO='"Optimised -O3"'   \
        -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size', '_get_pointset_ptr', '_get_pointset_size', '_about' ]" \
        -s NO_EXIT_RUNTIME=1          \
        -Winline         \
        -s TOTAL_MEMORY=30100100    \
        -s ABORTING_MALLOC=0 \
        -s ALLOW_MEMORY_GROWTH=1 \
        -s DISABLE_EXCEPTION_CATCHING=0  \
        -s DEMANGLE_SUPPORT=1 \
        -pedantic -std=c++14  \
    mcc2.cpp  \
        -o  ../build/mcc2.compiled.js
    exit 0
fi
#         -Wall  \

