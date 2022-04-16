#!/bin/bash
# test all builds
# Fresh test



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
cd implisolid
pwd
find $E2E -depth 2

pwd
# must run the internal one!
./demos/demo1-clonepull.sh
pwd
./demos/demo1-deploy.sh
