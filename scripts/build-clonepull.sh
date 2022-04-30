#!/bin/bash
set -ex
function assert_env_nonempty() {
  if [ ".$1" = "." ]; then
    echo "shell env is empty"; echo $2
    return 1
  fi
}

# run:  bash deploy-demo-1.sh

# This script downloads the prequisites

echo "11IMPLISOLID=$IMPLISOLID"
# args:
assert_env_nonempty $IMPLISOLID "env-argument IMPLISOLID= not specified"
#assert_env_nonempty $BASELOC1
#export IMPLISOLID=$BASELOC1/implisolid
# arg: pwd (unused)
assert_env_nonempty $CACHE_TEMP "env-argument CACHE_TEMP= not specified"

# not used?!: BASELOC2, BASELOC3

# 1. Clone the repo itself (including this script)
# 2. Run this to load libraries from github, etc

#REPO_ROOT=$(git rev-parse --show-toplevel)
#source $REPO_ROOT/demos/base-locations.sh

# # Parameters
# #export BASELOC1=/Users/$USER/cs/mp5
# #export BASELOC2=/Users/$USER/cs/mp5
# export BASELOC1=/Users/$USER/cs
# export BASELOC2=/Users/$USER/cs


#export IMPLISOLID=$BASELOC1/implisolid
# BASE_MP5_PRIVATE
#export MP5_PRIVATE=$BASELOC2/mp5-private

# manually:
# reset
#
# rm -rf $IMPLISOLID
# rm -rf $MP5_PRIVATE
# first time
#cd $BASELOC1
#git clone git@github.com:sohale/implisolid.git
#cd $BASELOC2
#git clone git@github.com:sohale/mp5-private.git # not used here

set -e

# target:
export BUILD_LOCATION=$IMPLISOLID/demos/build
mkdir -p $IMPLISOLID/demos/demo1  # why?  .DEMO_LOCATION
mkdir -p $BUILD_LOCATION

#sources:
export DEMO0=$IMPLISOLID/js_iteration_2/examples/mp5interactive
export JSI2=$IMPLISOLID/js_iteration_2
export JS_EX1=$IMPLISOLID/js_iteration_2/examples/js
export BUILT=$IMPLISOLID/docs/implisolid-build

# clones the repo and submodules
# pulls latest version
# pull latest dsocker
# ...

prime_docker() {
  echo "docker pull emscripten/emsdk"
  docker pull emscripten/emsdk
}
printf "\n\n\n"

cd $BUILD_LOCATION
pwd
mkdir -p $BUILD_LOCATION/lib
export LIB_FOLDER=$IMPLISOLID/demos/build/lib

get_boost() {
    echo "downloading library: Boost (1.75.0)"

    # check https://www.boost.org/users/download/
    # old: boost_1_61_0

    # two alternatives
    #ORIG_IMPLISOLID=$IMPLISOLID/..
    #CACHE_TEMP=$ORIG_IMPLISOLID/demos/build

    mkdir -p $CACHE_TEMP
    ls -1 $CACHE_TEMP/boost.tar.gz >/dev/null || \
    wget https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz -O $CACHE_TEMP/boost.tar.gz

    #wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz -O boost.tar.gz
    #echo rem || wget https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz -O boost.tar.gz
    cp $CACHE_TEMP/boost.tar.gz ./boost.tar.gz

    # boost_1_75_0.tar.gz
    echo Unzipping boost\'s .tar.gz
    ## gunzip -c boost_1_75_0.tar.gz | tar xopf -
    ## gunzip -c boost.tar.gz | tar xopf -
    gunzip -c boost.tar.gz | tar xopf -
    # # boost_1_75_0
    mv -vn ./boost_1_75_0 ./lib/
    #mv -vn ./boost_1_75_0/boost ./lib/boost

    # todo: wget filename
    # wget --server-response -q -O - "https://very.long/url/here" 2>&1 |   grep "Content-Disposition:" | tail -1 |   awk 'match($0, /filename=(.+)/, f){ print f[1] }' )

}
get_boost
export BOOST=$LIB_FOLDER/boost_1_75_0
echo "Boost downloaded in: $BOOST"


get_eigen() {
  echo "downloading library: Eigen  (3.?.?)"
  # see https://eigen.tuxfamily.org/index.php?title=Main_Page
  # seems to be 3.3.9

  # todo: 831133cc =  3.4.0-rc1
  # master latest was: 853a5c4b843a3f1de5de2a25429eefd62dbd153a

  pushd .
  mkdir -p $CACHE_TEMP; cd $CACHE_TEMP

  # two alternatives
  #git clone https://gitlab.com/libeigen/eigen.git  --depth 1
  #git clone https://github.com/libigl/eigen  --depth 1
  ls -1 $CACHE_TEMP/eigen >/dev/null || \
  git clone https://github.com/libigl/eigen  --depth 1

  popd

  mkdir -p ./eigen
  cp -R $CACHE_TEMP/eigen .

  #mv -vn ./eigen/Eigen ./lib/eigen
  mv -vn ./eigen ./lib/
}

get_eigen
export EIGEN=$LIB_FOLDER/eigen
echo "Eigen downloaded in: $EIGEN"

prime_docker

# prepare

# todo:
#  submodules

# first time: (`clone`s sub-modules)
# git submodule update --init --recursive
# ?
# git submodule update --recursive --remote

# not first time:
# git pull --recurse-submodules

# Look at the branches !
# * [new branch]      assert_fixing     -> origin/assert_fixing
# * [new branch]      optimization_task -> origin/optimization_task
# * [new branch]      researchOnSDF     -> origin/researchOnSDF
