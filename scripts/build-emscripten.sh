#!/bin/bash
set -eux

function assert_env_nonempty() {
  if [ ".$1" = "." ]; then
    echo "shell env is empty"; echo $2
    return 1
  fi
}

# todo: remove $SCRIPTS_DIR
# args:

#assert_env_nonempty $SCRIPTS_DIR "env-argument SCRIPTS_DIR missing. Must contain \$SCRIPTS_DIR/build_configuration.sh"
# build script folder, different to deploy etc

assert_env_nonempty $IMPLISOLID "env-argument IMPLISOLID ..."
#export IMPLISOLID=$BASELOC1/implisolid

assert_env_nonempty $BUILD_LOCATION "env-argument BUILD_LOCATION ..."
assert_env_nonempty $LIB_FOLDER "env-argument LIB_FOLDER ..."

#export BUILD_LOCATION=$IMPLISOLID/demos/build #eliminated
#export LIB_FOLDER=$BUILD_LOCATION/lib


# SCRIPTS_DIR is not really used

# usage: IMPLISOLID=$IMPLISOLID SCRIPTS_DIR=$IMPLISOLID/scripts ./scripts/demo1-build.sh

set -e

# BUILD_LOCATION = where compiled file wil be stored
# LIB_FOLDER = where to find libraries
function old_pattern() {
    IMPLISOLID=$IMPLISOLID source ./scripts/build_configuration.sh
    # output: BUILD_LOCATION,LIB_FOLDER
    assert_env_nonempty $BUILD_LOCATION "env-argument BUILD_LOCATION ..."
    assert_env_nonempty $LIB_FOLDER "env-argument LIB_FOLDER ..."
}

# Does not give you the folder, it tells you the name, that it wuol dbe in: found in, put in, etc.
mkdir -p $BUILD_LOCATION; ls $BUILD_LOCATION >/dev/null

# Parameters, from the specific configuration (relative locaation of build, lib, etc):

printf "BUILD_LOCATION:$BUILD_LOCATION, \nLIB_FOLDER:$LIB_FOLDER\n"

# expects in $LIB_FOLDER the following only : eigen, boost_1_75_0

# todo: tidy up, move up
#sources:
# LIB_FOLDER, MAIN_SOURCE_FOLDER
# targets:
# BUILD_LOCATION

export DEFAULT_MAIN_SOURCE_CPP_FILE="js_iteration_1/mcc2.cpp"
export DEFAULT_TARGET_FILENAME="mcc2.compiled.js"

#unfortunately it has to be one folder higher, becaausee both js_iteration_1 & ../js_iteration_2 are used.
export MAIN_SOURCE_FOLDER=$IMPLISOLID
# Main file is: at $MAIN_SOURCE_FOLDER/js_iteration_1/mcc2.cpp
export MAIN_SOURCE_CPP_FILE="${MAIN_SOURCE_CPP_FILE:-$DEFAULT_MAIN_SOURCE_CPP_FILE}"
export TARGET_FILENAME="${TARGET_FILENAME:-$DEFAULT_TARGET_FILENAME}"

expect_file  "$MAIN_SOURCE_FOLDER/$MAIN_SOURCE_CPP_FILE"

echo "source file: $MAIN_SOURCE_CPP_FILE"
echo "target file: $TARGET_FILENAME"

# old: boost_1_61_0
#BOOST_FOLDER="boost_1_75_0/boost"
#EIGEN_LIB_FOLDER="/eigen/Eigen"
BOOST_FOLDER="boost_1_75_0"
# todo: rename: BOOST_FOLDER -> BOOST_LIB_SUBFOLDER
EIGEN_LIB_FOLDER="eigen"
# todo: rename EIGEN_LIB_FOLDER -> EIGEN_LIB_SUBFOLDER
AUTODIFF_LIB_SUBFOLDER="autodiff"
# $LIB_FOLDER/autodiff/autodiff/forward/dual.hpp

source $IMPLISOLID/scripts/bash-utils.sh

# does pwd matter?
cd $IMPLISOLID/js_iteration_1


# todo: move up
# ./build/lib/boost_1_75_0/boost/array.hpp
expect_file  "$LIB_FOLDER/$BOOST_FOLDER/boost/array.hpp"
# cat ./build/lib/eigen/Eigen/src/Core/MatrixBase.h
expect_file  "$LIB_FOLDER/$EIGEN_LIB_FOLDER/Eigen/src/Core/MatrixBase.h"
expect_file "$LIB_FOLDER/$AUTODIFF_LIB_SUBFOLDER/autodiff/forward/dual.hpp"

export OPTIM=1
export DEV=2
#export WORKER=3
#export BITCODE=4

MODE=$DEV

export WASM=0

[[ $OSTYPE == 'darwin'* ]] || echo "Warning: This bash script only tested on MacOS. MacOS-specific code"
# because of the `if` `then` `fi` conditions

if [ $MODE -eq $OPTIM ]
then

    echo "** optimised mode **"

#         -I /src-lib/$AUTODIFF_LIB_SUBFOLDER \

    export CLI_ARGS=" \
        -I /src-lib/$BOOST_FOLDER \
        -I /src-lib/$EIGEN_LIB_FOLDER \
        -O3   \
        -Oz \
        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
        -s NO_EXIT_RUNTIME=1          \
        -Winline         \
        -s TOTAL_MEMORY=30146560    \
        -s ABORTING_MALLOC=0 \
        -s ALLOW_MEMORY_GROWTH=1 \
        -s DISABLE_EXCEPTION_CATCHING=0  \
        -s DEMANGLE_SUPPORT=1 \
        -Wno-dollar-in-identifier-extension \
        -pedantic -std=c++14 \
        -s WASM=0 \
        "

        # -s AGGRESSIVE_VARIABLE_ELIMINATION=1
        # -s OUTLINING_LIMIT=100000 \
fi

# Custom flags:
#    BUILD_AS_WORKER is defined by me.

if [ $MODE -eq $DEV ]
then

    echo "** dev compiling mode **"

#         -I /src-lib/$AUTODIFF_LIB_SUBFOLDER \

    export CLI_ARGS=" \
        -I /src-lib/$BOOST_FOLDER  \
        -I /src-lib/$EIGEN_LIB_FOLDER \
        -s TOTAL_MEMORY=30146560 \
        -s ABORTING_MALLOC=0 \
        -s NO_EXIT_RUNTIME=1 \
        -s DEMANGLE_SUPPORT=1 \
        -s ASSERTIONS=1 \
        -s BUILD_AS_WORKER=0 -DNOT_BUILD_AS_WORKER  \
        -DMORE_ABOUT_INFO='\"DEV\"'   \
        -Wno-dollar-in-identifier-extension \
        -pedantic -std=c++14 \
        -s WASM=0 \
        -s NO_DISABLE_EXCEPTION_CATCHING \
        "
        # -fexceptions \
fi

set -e

# tested on emscripten/emsdk:2.0.22
# tested on emscripten/emsdk:3.1.8  # detected a flaw

EXPORTED_FUNCTIONS="['_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size', '_get_pointset_ptr', '_get_pointset_size', '_build_geometry_u', '_about' ]"
#EXPORTED_FUNCTIONS="['_main', '_about', '_about2'  ]"
#EXPORTED_FUNCTIONS="['_main', '_about2'  ]"

docker run \
  --rm \
  -v $MAIN_SOURCE_FOLDER:/src \
  -v $LIB_FOLDER:/src-lib \
  -v $BUILD_LOCATION:/build \
  -u $(id -u):$(id -g) \
  emscripten/emsdk \
    emcc \
      $CLI_ARGS \
      -s EXPORTED_FUNCTIONS="$EXPORTED_FUNCTIONS" \
      -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' \
      /src/$MAIN_SOURCE_CPP_FILE \
      -o /build/$TARGET_FILENAME

# The compiled file can be found here:
ls "$BUILD_LOCATION/$TARGET_FILENAME"

expect_file  "$BUILD_LOCATION/$TARGET_FILENAME"

# todo:
# emcc: warning: EXTRA_EXPORTED_RUNTIME_METHODS is deprecated, please use EXPORTED_RUNTIME_METHODS instead [-Wdeprecated]

#
# More notes on usage of compiler
#
# docker run \
#  --rm \
#  -v $MAIN_SOURCE_FOLDER:/src \
#  -v $BUILD_LOCATION:/build \
#  -u $(id -u):$(id -g) \
#  emscripten/emsdk \
#
# emcc helloworld.cpp -o helloworld.js
# emcc mcc2.cpp  -o  ../build/mcc2.compiled.js

# emcc?
# em++?
#emsdk list
#emsdk list --old
#   --profiling
