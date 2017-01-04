#!/bin/bash

# Deployment:
# IMPLISOLID=/src
# Local:
#IMPLISOLID=/home/$USER/mp5-private/implisolid


if [ -z "$IMPLISOLID" ]; then
    echo "deploymwnt mode"
    IMPLISOLID=/src
fi

OPTIM=1
DEV=2
WORKER=3
BITCODE=4

MODE=1


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
            MODE=$DEV
            shift
            ;;
        -o|--opt)
            MODE=$OPTIM
            shift
            ;;
        -b|--bitcode)
            MODE=$BITCODE
            shift
            ;;
        -w|--worker)
            MODE=$WORKER
            shift
            ;;
        --)
            shift; break ;;
        *) echo "Internal Error"; exit 1 ;;
    esac
done



ls /lib/

cd $IMPLISOLID/js_iteration_1

BOOST_FOLDER="/lib/boost_1_61_0"
EIGEN_LIB_FOLDER="/lib/eigen"

if [ $OPTIM -eq $MODE ]
then

echo "** optimised **"

em++ \
    -I $BOOST_FOLDER \
    -I $EIGEN_LIB_FOLDER \
    -O3   \
    -Oz \
    -s OUTLINING_LIMIT=100000 \
    -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
    -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size', '_get_pointset_ptr', '_get_pointset_size', '_build_geometry_u', '_about' ]" \
    -s NO_EXIT_RUNTIME=1          \
    -Winline         \
    -s TOTAL_MEMORY=30100100    \
    -s ABORTING_MALLOC=0 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s DISABLE_EXCEPTION_CATCHING=0  \
    -s DEMANGLE_SUPPORT=1 \
    -Wno-dollar-in-identifier-extension \
    -pedantic -std=c++14  \
mcc2.cpp  \
    -o  ../build/mcc2.compiled.js
exit 0

fi


if [ $MODE -eq $DEV ]
then

    echo "** dev compiled **"

    em++ -I $BOOST_FOLDER  \
         -I $EIGEN_LIB_FOLDER \
        -s TOTAL_MEMORY=30100100 \
        -s ABORTING_MALLOC=0 \
        -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size', '_get_pointset_ptr', '_get_pointset_size', '_build_geometry_u', _about' ]" \
        -s NO_EXIT_RUNTIME=1 \
        -s DEMANGLE_SUPPORT=1 \
        -s ASSERTIONS=1 \
        -s BUILD_AS_WORKER=0 -DNOT_BUILD_AS_WORKER  \
        -DMORE_ABOUT_INFO='"DEV"'   \
        -pedantic -std=c++14 \
        mcc2.cpp -o ../build/mcc2.compiled.js
#        --profiling \

exit 0
fi

echo "** none **"
