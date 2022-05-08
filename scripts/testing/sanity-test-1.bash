#!/usr/bin/env bash
set -eux

export REPO_ROOT=$(git rev-parse --show-toplevel)
#cd  $REPO_ROOT/scripts/testing

source $REPO_ROOT/scripts/bash-utils.sh

#export BASEPATH="$HOME/cs"
#export IMPLISOLID="$BASEPATH/implisolid"
export IMPLISOLID="$REPO_ROOT"

[[ $OSTYPE == 'darwin'* ]] || echo "Warning: This bash script only tested on MacOS. MacOS-specific code"
# because of the `if` `then` `fi` conditions

export BUILD_LOCATION="$IMPLISOLID/build"
export LIB_FOLDER="$BUILD_LOCATION/lib"


export TARGET_FILENAME="mcc2.compiled.js"

# COMPILED, COMPILED_FILE, TARGET_FILE
export COMPILED_FILE="$BUILD_LOCATION/$TARGET_FILENAME"

# echo "Skipping build" || \
MAKE_HAPPEN  "$COMPILED_FILE" || {
                  # The MAKE_HAPPEN pattern
    time \
       LIB_FOLDER="$IMPLISOLID/build/lib" \
        BUILD_LOCATION=$BUILD_LOCATION   \
          bash "$IMPLISOLID/scripts/build-emscripten.sh"
}
expect_file  "$HOME/cs/implisolid/build/$TARGET_FILENAME"
expect_file  "$COMPILED_FILE"
# export COMPILED_FILE=$HOME/cs/implisolid/build/mcc2.compiled.js

node --version
# tested on v12.22.12

# node $COMPILED_FILE

node --trace-uncaught $REPO_ROOT/scripts/testing/sanity1.js  $COMPILED_FILE
