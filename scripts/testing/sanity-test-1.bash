#!/usr/bin/env bash
set -eux

source ../bash-utils.sh

function expect_file() {
    export FILE="$1"
    assert_env_nonempty $FILE "specify a filepath/name"
    if test -f "$FILE"; then
        # file exists, fine
        return 0
    else
        echo "$FILE does not exist. breaking"
        return -1
    fi
}

export BASEPATH="$HOME/cs"
export IMPLISOLID="$BASEPATH/implisolid"

echo || time \
   LIB_FOLDER="$IMPLISOLID/demos/build/lib" \
    BUILD_LOCATION="$IMPLISOLID/demos/build"   \
      bash "$IMPLISOLID/scripts/build-emscripten.sh"

expect_file  $HOME/cs/implisolid/demos/build/mcc2.compiled.js

export COMPILED=$HOME/cs/implisolid/demos/build/mcc2.compiled.js
node --version
# tested on v12.22.12

node $COMPILED

node sanity1.js $COMPILED
