
#!/bin/bash

# run:  bash deploy-demo-1.sh
set -e

# target:
export DEMO_LOCATION=/Users/$USER/cs/mp5/implisolid/demos/demo1
export BUILD_LOCATION=/Users/$USER/cs/mp5/implisolid/demos/build

#sources:
export MP5_PRIVATE=/Users/$USER/cs/mp5/mp5-private
export IMPLISOLID=/Users/$USER/cs/mp5/implisolid
export DEMO0=$IMPLISOLID/js_iteration_2/examples/mp5interactive
export JSI2=$IMPLISOLID/js_iteration_2
export JS_EX1=$IMPLISOLID/js_iteration_2/examples/js
export BUILT=$IMPLISOLID/docs/implisolid-build

# clones the repo and submodules
# pulls latest version
# pull latest dsocker
# ...

prime_docker() {
  docker pull emscripten/emsdk
}
echo "\n\n\n"

cd $BUILD_LOCATION
pwd
mkdir -p lib

get_boost() {
    echo "downloading library: Boost (1.75.0)"

    # check https://www.boost.org/users/download/
    # old: boost_1_61_0

    wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz -O boost.tar.gz

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

get_eigen() {
  echo "downloading library: Eigen  (3.?.?)"
  # see https://eigen.tuxfamily.org/index.php?title=Main_Page
  # seems to be 3.3.9

  # todo: 831133cc =  3.4.0-rc1
  # master latest was: 853a5c4b843a3f1de5de2a25429eefd62dbd153a
  git clone https://gitlab.com/libeigen/eigen.git  --depth 1
  #mv -vn ./eigen/Eigen ./lib/eigen
  mv -vn ./eigen ./lib/
}

get_eigen

# prepare

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
