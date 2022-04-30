#!/bin/bash
# test all builds
# Fresh test
set -ex

# args:
# * pwd (via $ORIG_REPO_ROOT)
# mmust be executed fromo within a folder in the main implisolid repo
# todo : another version for install-dev

#docker pull frolvlad/alpine-gxx
#docker pull groovy
#docker pull bash
#docker pull emscripten/emsdk

# todo:
# Linux: Using docker bash
# Mac

# Python
# all targets: npm

# move out of ./demos

# e2e
# for latest?

# also can be used for installing (installing dev-environment)

# latest: Pulling from frolvlad/alpine-gxx Digest: sha256:d31e4365fd9e0ec840565259d618bff9d4b8387a4d8c1fdea0eddb6d42f8c95e
# latest: Pulling from library/groovy Digest: sha256:49156d340390d288e2c94403dc78d535a679020f30485ceadafb1b59b8191eca
# latest: Pulling from library/bash Digest: sha256:b3abe4255706618c550e8db5ec0875328333a14dbf663e6f1e2b6875f45521e5
# docker.io/library/bash:latest
# tested on emscripten/emsdk:2.0.22
# tested on emscripten/emsdk:3.1.8  # detected a flaw


function __current_script_dir_func0 () {
  cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd
}
source $(__current_script_dir_func0)/bash-utils.sh

set -ex
export ORIG_REPO_ROOT=$(git rev-parse --show-toplevel)
echo "ORIG_REPO_ROOT :::: $ORIG_REPO_ROOT"
assert_env_nonempty $ORIG_REPO_ROOT "ORIG_REPO_ROOT=$ORIG_REPO_ROOT   implisolid repo not found in current directory $(pwd)"

#cd $ORIG_REPO_ROOT; mkdir e2e-sandbox-temp
E2E=$ORIG_REPO_ROOT/e2e-sandbox-temp
rm -rf $E2E
# #### hard reset done #####

mkdir $E2E
cd $E2E
# only updated after actually pushing => requires branch name! active branch name: hot branch: one neing processed. stil hot
git clone --recurse-submodules git@github.com:sohale/implisolid.git
# todo: from local:
#rsync -r $ORIG_REPO_ROOT $E2E
# recursive

cd implisolid
pwd
find $E2E -maxdepth 2

NEWREPO_BASE="$E2E"
NEW_REPO="$E2E/implisolid"
# NEW_REPO=NEWREPO_ROOT

# CACHE_TEMP=$ORIG_IMPLISOLID/demos/build
CACHE_TEMP="$ORIG_REPO_ROOT/temp"
mkdir -p "$CACHE_TEMP"



# export SOHALE_IO_REPO=$BASELOC3/sohale.github.io/
# export BASE_MP5_PRIVATE=$BASELOC2/mp5/mp5-private

# Targets:
#1. for sohale.github.io (uses own repo, aka "remote" or static) -> ./docs/implisolid-build/demo1
#       remote/public deploy (to github-pages)
#2. for separate (liquid) (runtime-created) folder (aka local! deploy) -> ./demos/demo1/
#3. for mp5_private  --> /home/$USER/mp5-private/implisolid

# Two alternatives: pre-buillt, and the new-built:
# simplified:
export DEPLOY_LOCATION=$NEW_REPO/demos/demo1
export DEPLOY_LOCATION=$NEW_REPO/docs/implisolid-build/demo1
# was: export DEPLOY_LOCATION=$NEW_REPO/docs/implisolid-build
export DEPLOY_LOCATION="$ORIG_REPO_ROOT/newapp"


pwd
# must run the internal one!
CACHE_TEMP="$CACHE_TEMP" \

    IMPLISOLID="$NEW_REPO"  \
    LIB_FOLDER="$IMPLISOLID/demos/build/lib"  BUILD_LOCATION="$IMPLISOLID/demos/build" \
    ./scripts/build-clonepull.sh

pwd

IMPLISOLID="$NEW_REPO"  LIB_FOLDER="$IMPLISOLID/demos/build/lib" \
    BUILD_LOCATION="$IMPLISOLID/demos/build"   \
    bash ./scripts/build-emscripten.sh

pwd
IMPLISOLID_REPO="$NEW_REPO"  BUILD_LOCATION="$NEW_REPO/demos/build" \
    DEPLOY_LOCATION="$DEPLOY_LOCATION" \
    bash ./scripts/demos/demo1/demo1-deploy.sh
# also runs demos/demo1-run-local.sh

pwd
APP_RUN_LOCATION="$DEPLOY_LOCATION" \
    bash ./scripts/demos/demo1/demo1-run-local.sh

echo "All successful."
return 0

<< ////

Scripts:
        # based on how they can move: 1. e2e/overall  2.

    e2e-test-builds.bash
    //    * is not for demo => not in ./demo
    // also do: install-dev
    // also a script foor docker? for test under linux
    demo1-clonepull.sh

    #use the configuration:
        * is not for demo => not in ./demo
        * demo1-build.sh // no! does not use the configuration of demo

        * demo1-deploy.sh
        * demo1-run-local.sh
        * build_configuration.sh

    * bash-utils.sh

    # deprecated
       * base-locations.sh
       * macos-specific.bash

folders:
    build
    demo1

conceptual places (folders):
    demo1 (where deploy is copied. no scripts.)
    demo2 (potenially other demos1)
    e2e
    scripts
    build (where build is located)
    mainrepo
    other repos?

    how about examples?

    build-binary-publish?
    two others

    python demo
    npm
    ...

export BUILD_LOCATION=$IMPLISOLID_REPO/demos/build

ideal names:
    $SCRIPTS_DIR/bash-utils.sh

////


<< ////
base-locations.sh:
# use with `source` only

echo "invalid"
export BASELOC1=
export BASELOC2=
export BASELOC3=
return

# export USER=a9858770
export USER_HOME=/Users/$USER

# Parameters
#export BASELOC1=$USER_HOME/cs/mp5
#export BASELOC2=$USER_HOME/cs/mp5
export BASELOC1=$USER_HOME/cs
export BASELOC2=$USER_HOME/cs
export BASELOC3=$USER_HOME/cs

# $BASELOC1 -> implisolid
# $BASELOC2 -> mp5-private
# $BASELOC3 -> sohale.github.io

# /Users/$USER/cs/mp5/implisolid
# mycomputer-specific

[[ $OSTYPE == 'darwin'* ]] || "Error: MacOS-specific code"

# Base locations for implisolic, mp5-private, sohale.github.io
# Other scrips use this as the seed parameters

////
