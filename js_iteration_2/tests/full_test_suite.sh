GTEST_ROOT=$GOOGLETEST_ROOT
fullfilename=$1 # get first argument as input file
BOOST_FOLDER=/usr/local/include/boost_1_61_0/


if [[ -n $EM_PREPARE ]]
then
    echo ""
else
    echo "Error: You need to first run \"source em_prepare.sh\""
    echo
    exit 1;
fi

if [ ! -d "$GTEST_ROOT" ]; then
    echo "Error: googletest directory not found at $GTEST_ROOT"
    echo
    exit 1;
fi

if [ ! -d "$BOOST_FOLDER" ]; then
    echo "Error: Boost 1.61.0+ not found at $BOOST_FOLDER"
    echo
    exit 1;
fi


mkdir -p build
em++    -I $BOOST_FOLDER    \
        -I $EIGEN_LIB_FOLDER \
        -I $GTEST_ROOT/googletest/include \
        -s ASSERTIONS=1 \
        -pedantic -std=c++14 \
        main_all_tests.cpp  test_*.cpp\
        ${GTEST_ROOT}/build/googlemock/gtest/libgtest.a -o \
        build/testsuit.compiled.js \
        && echo "compile success." \
        && node build/"$filename".compiled.js \
        && echo "test executed." \
        || echo "something went wrong. tests were not executed."




exit
bash build_tests.sh  make_random_pm1.cpp

bash build_tests.sh    timer_test.cpp

#bash finitediff_test0.cpp
bash build_tests.sh    primitives_test.cpp
#bash finitediff_test.cpp
bash build_tests.sh    test_crisp_subtract.cpp
bash build_tests.sh    gradient_test.cpp
bash build_tests.sh    implicit_function_output.cpp
bash build_tests.sh    test_vectorised_algorithms.cpp
bash build_tests.sh    make_random_pm1.cpp
bash build_tests.sh    perf_crisp_subtract.cpp

bash build_tests.sh    test_svd.cpp
