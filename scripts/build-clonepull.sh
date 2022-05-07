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

# target:
#export BUILD_LOCATION=$IMPLISOLID/build
assert_env_nonempty $BUILD_LOCATION "env-argument BUILD_LOCATION= not specified"
assert_env_nonempty $LIB_FOLDER "env-argument LIB_FOLDER= not specified"

source $IMPLISOLID/scripts/bash-utils.sh

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
#export BUILD_LOCATION=$IMPLISOLID/build
#mkdir -p $IMPLISOLID/demos/demo1  # why?  .DEMO_LOCATION
mkdir -p $BUILD_LOCATION

#sources:
export DEMO0=$IMPLISOLID/examples/mp5interactive
export JSI2=$IMPLISOLID/js_iteration_2
export EX_JSLIB=$IMPLISOLID/examples/js-lib
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

export QQ="$BUILD_LOCATION"
cd $BUILD_LOCATION
cd $QQ # BUILD_LOCATION

pwd
#mkdir -p $BUILD_LOCATION/lib
#export LIB_FOLDER=$IMPLISOLID/build/lib
#export LIB_FOLDER=$BUILD_LOCATION/lib
mkdir -p $LIB_FOLDER

get_boost() {
    echo "downloading library: Boost (1.75.0)"

    # check https://www.boost.org/users/download/
    # old: boost_1_61_0

    # two alternatives
    #ORIG_IMPLISOLID=$IMPLISOLID/..
    #CACHE_TEMP=$ORIG_IMPLISOLID/build

    mkdir -p $CACHE_TEMP
    ls -1 $CACHE_TEMP/boost.tar.gz >/dev/null || \
    wget https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz -O $CACHE_TEMP/boost.tar.gz

    #wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz -O boost.tar.gz
    #echo rem || wget https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz -O boost.tar.gz
    cp $CACHE_TEMP/boost.tar.gz  $QQ/boost.tar.gz

    # boost_1_75_0.tar.gz
    echo "Unzipping boost\'s .tar.gz"
    ## gunzip -c boost_1_75_0.tar.gz | tar xopf -
    ## gunzip -c boost.tar.gz | tar xopf -
    #gunzip -c boost.tar.gz | tar xopf -
    #pushd .
    #cd $QQ
    gunzip -c $QQ/boost.tar.gz | tar xopf - -C $QQ
    #popd
    # # boost_1_75_0
    mv -vn $QQ/boost_1_75_0 $LIB_FOLDER/
    #mv -vn $QQ/boost_1_75_0/boost $LIB_FOLDER/boost

    # todo: wget filename
    # wget --server-response -q -O - "https://very.long/url/here" 2>&1 |   grep "Content-Disposition:" | tail -1 |   awk 'match($0, /filename=(.+)/, f){ print f[1] }' )

}
export BOOST_FOLDER="boost_1_75_0"
MAKE_HAPPEN "$LIB_FOLDER/$BOOST_FOLDER/boost/array.hpp" \
  || {
      get_boost
     }
export BOOST="$LIB_FOLDER/$BOOST_FOLDER"
echo "Boost downloaded in: $BOOST"

# ./build/lib/boost_1_75_0/boost/array.hpp
expect_file  "$LIB_FOLDER/$BOOST_FOLDER/boost/array.hpp"


get_eigen() {
  echo "downloading library: Eigen  (3.?.?)"
  # see https://eigen.tuxfamily.org/index.php?title=Main_Page
  # seems to be 3.3.9

  # todo: 831133cc =  3.4.0-rc1
  # master latest was: 853a5c4b843a3f1de5de2a25429eefd62dbd153a

  # A benefit of this approach is that the newest docker will be pulled, and the respective updates will be compelled.

  #pushd .
  mkdir -p $CACHE_TEMP #; cd $CACHE_TEMP

  # two alternatives
  #git clone https://gitlab.com/libeigen/eigen.git  --depth 1
  #git clone https://github.com/libigl/eigen  --depth 1
  ls -1 $CACHE_TEMP/eigen >/dev/null || \
  git clone https://github.com/libigl/eigen  --depth 1  $CACHE_TEMP/eigen

  #popd

  mkdir -p $QQ/eigen
  cp -R $CACHE_TEMP/eigen $QQ

  #mv -vn $QQ/eigen/Eigen $LIB_FOLDER/eigen
  mv -vn $QQ/eigen $LIB_FOLDER/
}

export EIGEN_LIB_FOLDER="eigen"
MAKE_HAPPEN "$LIB_FOLDER/$EIGEN_LIB_FOLDER/Eigen/src/Core/MatrixBase.h" \
  || {
    get_eigen
  }
export EIGEN="$LIB_FOLDER/$EIGEN_LIB_FOLDER"
echo "Eigen downloaded in: $EIGEN"

# cat ./build/lib/eigen/Eigen/src/Core/MatrixBase.h
expect_file  "$LIB_FOLDER/$EIGEN_LIB_FOLDER/Eigen/src/Core/MatrixBase.h"

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
