echo "Don't run this script. It is instructions for manual installation"
exit

This file contains instructions for installing Google Test fraework for Emscripten.
Synopsis:
1- choose a folder
2- put files (git clone)
3- set em_prepare.sh
4- mkdir build
5- set EMSCRIPTEN_ROOT
6- make a symbolic link
7- cmake and make
8- check if libgtest.a is created
9- compile some tests

# go to an installation folder: for example:
export G1=~/mp5  # change according to this installation

cd $G1
git clone https://github.com/google/googletest.git
export GOOGLETEST_ROOT=${G1}/googletest
echo Add this line to em_prepare.sh
echo export GOOGLETEST_ROOT=~/mp5/googletest

cd $GOOGLETEST_ROOT
# cd googletest/
mkdir build
cd build/
emcmake cmake -Dgtest_disable_pthreads=ON ..
EMSCRIPTEN_ROOT=$EMSCRIPTEN
ln -s ${EMSCRIPTEN_ROOT}/system/lib/libcxxabi/include/cxxabi.h ${EMSCRIPTEN_ROOT}/system/local/include
emmake make

echo This should exist:
ls ${GOOGLETEST_ROOT}/build/googlemock/gtest/libgtest.a
# for me: /home/sohail/mp5/googletest/build/googlemock/gtest/libgtest.a

# compile:
# em++ -I ${GOOGLETEST_ROOT}/googletest/include/ [...]
#                           your_test.cpp ${GOOGLETEST_ROOT}/build/googlemock/gtest/libgtest.a -o your_test.js

#Then go to the following folder
#      /mp5-private/solidmodeler/js_iteration_2/tests$
# and run:
#      bash build_tests.sh timer_test.cpp
