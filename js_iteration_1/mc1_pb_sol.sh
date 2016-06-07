


# em++ -I ./boost_1_61_0/ -s EXPORTED_FUNCTIONS="['_main']" -pedantic -std=c++14 -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS -O3 vect_test_2.cpp -o vect_test_2.js

  # ~/Desktop/emsdk_portable/emscripten/master/em++ -I /usr/local/include/boost_1_61_0/  -s EXPORTED_FUNCTIONS="['_make_object']"  -pedantic -std=c++14  mc1_pb.cpp  -o mc1_pb.js && echo DONE
   # -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS \
  # ~/Desktop/emsdk_portable/emscripten/master/em++ -I /usr/local/include/Eigen/  -s EXPORTED_FUNCTIONS="['_make_object']"  -pedantic -std=c++14  mc1_eigen.cpp  -o mc1_eigen.js && echo DONE
    # -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS \


    # echo "@rem em++ -I /usr/local/inlcude/boost_1_61_0/ -s EXPORTED_FUNCTIONS="['_main']" -pedantic -std=c++14  -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  -O3     vect_test_2.cpp    -o vect_test_2.js"
    #   ~/Desktop/emsdk_portable/emscripten/master/em++ -I /usr/local/include/boost_1_61_0/ -s EXPORTED_FUNCTIONS="['_make_object', '_main']" -s NO_EXIT_RUNTIME=1 -s ASSERTIONS=1 -pedantic -std=c++14 mc2_pb.cpp -o mc2_pb.cpp.js
  #  echo "@rem em++ -I /usr/local/inlcude/boost_1_61_0/ -s EXPORTED_FUNCTIONS="['_main']" -pedantic -std=c++14  -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  -O3     vect_test_2.cpp    -o vect_test_2.js"
#      ~/Desktop/emsdk_portable/emscripten/master/em++ -I /usr/local/include/boost_1_61_0/ -s EXPORTED_FUNCTIONS="['_make_object', '_main']" -s NO_EXIT_RUNTIME=1 -s ASSERTIONS=1 -pedantic -std=c++14 mc1_eigen.cpp -o mc1_eigen.cpp.js
      ~/Desktop/emsdk_portable/emscripten/master/em++ -I /usr/local/include/Eigen/ -s EXPORTED_FUNCTIONS="['_make_object', '_main']" -s NO_EXIT_RUNTIME=1 -s ASSERTIONS=1 -pedantic -std=c++14 mc2_armadillo.cpp -o mc2_armadillo.cpp.js
    # ::     -s EXPORTED_FUNCTIONS="['_main', '_make_object']" ^
    # ::     -s EXPORTED_FUNCTIONS="['_main']" ^
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
