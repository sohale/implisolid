#!/bin/bash

# echo "@rem em++ -I /usr/local/inlcude/boost_1_61_0/ -s EXPORTED_FUNCTIONS="['_main']" -pedantic -std=c++14  -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  -O3     vect_test_2.cpp    -o vect_test_2.js"


#BOOST_FOLDER=/usr/local/include/boost_1_61_0
#BOOST_FOLDER=./boost_1_61_0


#if [ -z ${EM_PREPARE+x} ]
if [[ -n "$EM_PREPARE" ]]
then
    echo ""
else
    echo "First run \"source em_prepare.sh\""
    exit
fi

case "$1" in
    0)
    OPTIM=0
    ;;
    1)
    OPTIM=1
    ;;
    *)
    echo "Select how to compile [0 / 1] or press Enter for default(Development)"
    echo

    echo Development Version: 0
    echo Optimized Version: 1
    read OPTIM
    ;;
esac

case "$OPTIM" in
    0)
    echo "Will now compile in Development Version, optimizations are turned off"
    OPTIM=0
    ;;
    1)
    echo "Will now compile in Optimized Version, optimizations are turned on"
    OPTIM=1
    ;;
    *)
    OPTIM=0
    echo "Running default, Development Version"
    ;;
esac

if [ $OPTIM -eq 0 ]
then

    echo "** dev version **"

    em++ -I $BOOST_FOLDER  \
        -s TOTAL_MEMORY=30100100 \
        -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size' ]" \
        -s NO_EXIT_RUNTIME=1 \
        -s DEMANGLE_SUPPORT=1 \
         -s ASSERTIONS=1 \
        -pedantic -std=c++14 \
        --profiling \
        mcc2.cpp -o mcc2.compiled.js

fi


if [ $OPTIM -eq 1 ]
then
    echo "** optimised version **"

    em++ \
            -I $BOOST_FOLDER   \
            -O3   \
            --profiling     \
            -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  \
            -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr',   '_set_object', '_unset_object', '_set_x', '_unset_x', '_calculate_implicit_values', '_get_values_ptr', '_get_values_size', '_calculate_implicit_gradients', '_get_gradients_ptr', '_get_gradients_size' ]" \
            -s NO_EXIT_RUNTIME=1          \
            -Winline         \
            -s TOTAL_MEMORY=30100100    \
            -pedantic -std=c++14  \
        mcc2.cpp  \
            -o  mcc2.compiled.js
fi

# copy the path for em++
# echo "@rem em++ -I /usr/local/inlcude/boost_1_61_0/ -s EXPORTED_FUNCTIONS="['_main']" -pedantic -std=c++14  -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  -O3     vect_test_2.cpp    -o vect_test_2.js"
#   ~/emsdk_portable/emscripten/master/em++ -I /usr/local/include/boost_1_61_0/  \
#   -s TOTAL_MEMORY=30100100 -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr']" -O3 -s NO_EXIT_RUNTIME=1 -s ASSERTIONS=1 -s DEMANGLE_SUPPORT=1 -pedantic -std=c++14 mcc2.cpp -o mcc2.cpp.js

# ::     -s EXPORTED_FUNCTIONS="['_main', '_make_object']" ^
# ::     -s EXPORTED_FUNCTIONS="['_main']" ^1
# :: @rem     -O3   ^
# :: @rem     -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
#
#
# ::     -O3   ^
# ::     -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
#
#
# ::        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
# ::        -O3                           ^
#

# em++ ^
#         -I C:\sohail\March\emscripten\boost_1_61_0\   ^
#         -s ASSERTIONS=1               ^
#         -s DEMANGLE_SUPPORT=1   ^
#         -s TOTAL_MEMORY=30100100    ^
#                                     ^
#         -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr']" ^
#         -s NO_EXIT_RUNTIME=1          ^
#                                           ^
#         -pedantic -std=c++14  ^
#     mcc2.cpp  ^
#         -o  mcc2.cpp.js
