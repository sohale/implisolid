#!/bin/bash
# test all builds
# Fresh test


# args:
# * pwd (via $REPO_ROOT)

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

set -ex
REPO_ROOT=$(git rev-parse --show-toplevel)
#cd $REPO_ROOT; mkdir e2e-sandbox-temp
E2E=$REPO_ROOT/e2e-sandbox-temp
rm -rf $E2E
# #### hard reset done #####

mkdir $E2E
cd $E2E
# only updated after actually pushing => requires branch name! active branch name: hot branch: one neing processed. stil hot
git clone git@github.com:sohale/implisolid.git
# todo: from local:
#rsync -r $REPO_ROOT $E2E
# recursive

cd implisolid
pwd
find $E2E -depth 2

NEWREPO_ROOT=$E2E/implisolid
NEWREPO_BASE=$E2E

# todo: remove, to make it explicit.
export BASELOC1=$NEWREPO_BASE
export BASELOC2=$NEWREPO_BASE
export BASELOC3=$NEWREPO_BASE

pwd
# must run the internal one!
BASELOC1=$BASELOC1 ./demos/demo1-clonepull.sh
pwd
BASELOC1=$BASELOC1 ./demos/demo1-deploy.sh

<< ////

Scripts:
        # based on how they can move: 1. e2e/overall  2.

    e2e-test-builds.bash
    demo1-clonepull.sh

    #use the configuration:
        * demo1-build.sh // no! does not use the configuration of demo

        * demo1-deploy.sh
        * demo1-run-local.sh
        * build_configuration.sh

    * utils.sh

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

    build-binary-publish?
    two others

    python demo
    npm
    ...


ideal names:
    $SCRIPTS_DIR/utils.sh

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
