#!/bin/bash

# old names:
# recompile-for-demo1.sh
# demo1-build.sh

REPO_ROOT=$(git rev-parse --show-toplevel)
source $REPO_ROOT/demos/base-locations.sh

# Seed for parameter values:
export IMPLISOLID=$BASELOC1/implisolid
set -e

REPO_ROOT=$(git rev-parse --show-toplevel)
source $REPO_ROOT/demos/build_configuration.sh

# Parameters, from the specific configuration (relative locaation of build, lib, etc):
# target:
# DEMO_LOCATION, BUILD_LOCATION, LIB_FOLDER

printf "DEMO_LOCATION:$DEMO_LOCATION, \nBUILD_LOCATION:$BUILD_LOCATION, \nLIB_FOLDER:$LIB_FOLDER\n"

#export LIB_FOLDER=$IMPLISOLID/demos/build/lib
#LIB_DIR
echo "LIB_FOLDER: $LIB_FOLDER"
ls -alt $LIB_FOLDER


#echo $LIB_FOLDER
#ls $LIB_FOLDER
#exit

#sources:

#export IMPLISOLID=$BASELOC1/implisolid

# MacOS-specific
#export IMPLISOLID=$BASELOC1/implisolid
#IMPLISOLID=/home/$USER/mp5-private/implisolid

#export SOURCE_FOLDER=$IMPLISOLID/../js_iteration_1/mcc2.cpp

#export MAIN_SOURCE_FILE=$IMPLISOLID/js_iteration_1/mcc2.cpp
#export MAIN_SOURCE_FOLDER=$IMPLISOLID/js_iteration_1
#unfortunately it has to be one folder higher, becaausee both js_iteration_1 & ../js_iteration_2 are used.
export MAIN_SOURCE_FOLDER=$IMPLISOLID


# errors if doesnt exist
ls $MAIN_SOURCE_FILE >/dev/null

mkdir -p $BUILD_LOCATION
echo $BUILD_LOCATION
ls $BUILD_LOCATION >/dev/null


# BOOST_FOLDER="/lib/boost_1_61_0"
# EIGEN_LIB_FOLDER="/lib/eigen"


#BOOST_FOLDER="/src-lib/boost_1_75_0/boost"
#EIGEN_LIB_FOLDER="/src-lib/eigen/Eigen"
BOOST_FOLDER="/src-lib/boost_1_75_0"
EIGEN_LIB_FOLDER="/src-lib/eigen"



cd $IMPLISOLID/js_iteration_1


[[ $OSTYPE == 'darwin'* ]] || echo "Error: This bash script only tested on MacOS"

export OPTIM=1
export DEV=2
#export WORKER=3
#export BITCODE=4

MODE=$DEV

export WASM=0

if [ $MODE -eq $OPTIM ]
then

    echo "** optimised mode **"

    export CLI_ARGS=" \
        -I $BOOST_FOLDER \
        -I $EIGEN_LIB_FOLDER \
        -O3   \
        -Oz \
        -s OUTLINING_LIMIT=100000 \
        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
        -s NO_EXIT_RUNTIME=1          \
        -Winline         \
        -s TOTAL_MEMORY=30146560    \
        -s ABORTING_MALLOC=0 \
        -s ALLOW_MEMORY_GROWTH=1 \
        -s DISABLE_EXCEPTION_CATCHING=0  \
        -s DEMANGLE_SUPPORT=1 \
        -Wno-dollar-in-identifier-extension \
        -pedantic -std=c++14
        -s WASM=0 \
        "

    #mcc2.cpp  \
    #    -o  ../build/mcc2.compiled.js \

    #em++ $CLI_ARGS
    #exit 0
fi

# Custom flags:
#    BUILD_AS_WORKER is defined by me.

if [ $MODE -eq $DEV ]
then

    echo "** dev compiling mode **"

    export CLI_ARGS=" \
        -I $BOOST_FOLDER  \
        -I $EIGEN_LIB_FOLDER \
        -s TOTAL_MEMORY=30146560 \
        -s ABORTING_MALLOC=0 \
        -s NO_EXIT_RUNTIME=1 \
        -s DEMANGLE_SUPPORT=1 \
        -s ASSERTIONS=1 \
        -s BUILD_AS_WORKER=0 -DNOT_BUILD_AS_WORKER  \
        -DMORE_ABOUT_INFO='\"DEV\"'   \
        -pedantic -std=c++14 \
        -s WASM=0 \
        "

    #mcc2.cpp \
    #    -o ../build/mcc2.compiled.js "

    #   --profiling \

    #em++ $CLI_ARGS
    #exit 0
fi

pwd
echo \$CLI_ARGS:
echo $CLI_ARGS

set -e

#docker run \
#  --rm \
#  -v $MAIN_SOURCE_FOLDER:/src \
#  -v $BUILD_LOCATION:/build \
#  -u $(id -u):$(id -g) \
#  emscripten/emsdk \
#  emsdk list --old
#  #emsdk list

# tested on emscripten/emsdk:2.0.22
# tested on emscripten/emsdk:3.1.8  # detected a flaw

docker run \
  --rm \
  -v $MAIN_SOURCE_FOLDER:/src \
  -v $LIB_FOLDER:/src-lib \
  -v $BUILD_LOCATION:/build \
  -u $(id -u):$(id -g) \
  emscripten/emsdk \
    emcc \
      $CLI_ARGS \
      -s EXPORTED_FUNCTIONS="['_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size', '_get_pointset_ptr', '_get_pointset_size', '_build_geometry_u', '_about' ]" \
      -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' \
      /src/js_iteration_1/mcc2.cpp \
      -o /build/mcc2.compiled.js

# It can be found here:
ls $BUILD_LOCATION/mcc2.compiled.js

  #emcc helloworld.cpp -o helloworld.js

# emcc: warning: EXTRA_EXPORTED_RUNTIME_METHODS is deprecated, please use EXPORTED_RUNTIME_METHODS instead [-Wdeprecated]
