#!/usr/bin/env bash
set -eux

# why error?
export REPO_ROOT=$(git rev-parse --show-toplevel)
#cd  $REPO_ROOT/scripts/testing

source $REPO_ROOT/scripts/bash-utils.sh
export REPO_ROOT="$HOME/cs/implisolid"

#export BASEPATH="$HOME/cs"
#export IMPLISOLID="$BASEPATH/implisolid"
export IMPLISOLID="$REPO_ROOT"

[[ $OSTYPE == 'darwin'* ]] || echo "Warning: This bash script only tested on MacOS. MacOS-specific code"
# because of the `if` `then` `fi` conditions

export BUILD_LOCATION="$IMPLISOLID/build"
export LIB_FOLDER="$BUILD_LOCATION/lib"
# LIB_FOLDER= /Users/9858770/cs/implisolid/build/lib

# rm -rf /Users/9858770/cs/implisolid/build/lib/autodiff/

echo ">>>>>>>>>"

# ls -1  "$LIB_FOLDER/autodiff/autodiff/forward/dual.hpp"
MAKE_HAPPEN  "$LIB_FOLDER/autodiff/autodiff/forward/dual.hpp" || {
   echo "CLON"
    git clone https://github.com/autodiff/autodiff "$LIB_FOLDER/autodiff"
}
   echo "NCLON"
expect_file "$LIB_FOLDER/autodiff/autodiff/forward/dual.hpp"

export TARGET_FILENAME="autodiff.compiled.js"

# COMPILED, COMPILED_FILE, TARGET_FILE
export COMPILED_FILE="$BUILD_LOCATION/$TARGET_FILENAME"

rm -f $COMPILED_FILE

expect_file "$IMPLISOLID/examples/implicit-functions/cpp/autodiff-sample1.cpp.1"

# echo "Skipping build" || \
MAKE_HAPPEN  "$COMPILED_FILE" || {
                  # The MAKE_HAPPEN pattern
    time \
       LIB_FOLDER="$IMPLISOLID/build/lib" \
        BUILD_LOCATION=$BUILD_LOCATION   \
        MAIN_SOURCE_CPP_FILE="sandbox/autodiff/implicit-functions/cpp/autodiff-sample1.cpp" \
        TARGET_FILENAME=$TARGET_FILENAME \
          bash "$IMPLISOLID/scripts/build-emscripten.sh"
}
expect_file  "$HOME/cs/implisolid/build/$TARGET_FILENAME"
expect_file  "$COMPILED_FILE"
# export COMPILED_FILE=$HOME/cs/implisolid/build/mcc2.compiled.js

node --version
# tested on v12.22.12

# node $COMPILED_FILE

node --trace-uncaught $REPO_ROOT/sandbox/autodiff/autodiff-sanity1.js \
    $COMPILED_FILE
