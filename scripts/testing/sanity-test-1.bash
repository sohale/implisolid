#!/usr/bin/env bash
set -eux

export REPO_ROOT=$(git rev-parse --show-toplevel)
#cd  $REPO_ROOT/scripts/testing

source $REPO_ROOT/scripts/bash-utils.sh

#export BASEPATH="$HOME/cs"
#export IMPLISOLID="$BASEPATH/implisolid"
export IMPLISOLID="$REPO_ROOT"

# echo "Skipping build" || \
expect_file  $HOME/cs/implisolid/build/mcc2.compiled.js || {
                  # The MAKE_HAPPEN pattern
time \
   LIB_FOLDER="$IMPLISOLID/build/lib" \
    BUILD_LOCATION="$IMPLISOLID/build"   \
      bash "$IMPLISOLID/scripts/build-emscripten.sh"
}

expect_file  $HOME/cs/implisolid/build/mcc2.compiled.js

export COMPILED=$HOME/cs/implisolid/build/mcc2.compiled.js
node --version
# tested on v12.22.12

# node $COMPILED

node --trace-uncaught $REPO_ROOT/scripts/testing/sanity1.js $COMPILED
