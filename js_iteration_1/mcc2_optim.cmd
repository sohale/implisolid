@echo
@echo "--------------------------------------------------- optim ------------------------------------------------------"

em++ ^
        -I C:\sohail\March\emscripten\boost_1_61_0\   ^
        -O3   ^
        -DNDEBUG -DBOOST_UBLAS_NDEBUG -DBOOST_DISABLE_ASSERTS  ^
 ^
        -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry', '_get_f_ptr', '_get_v_ptr' ]" ^
        -s NO_EXIT_RUNTIME=1          ^
        --llvm-lto 1     ^
        -Winline         ^
        -s TOTAL_MEMORY=30100100    ^
 ^
        -pedantic -std=c++14  ^
    mcc2.cpp  ^
        -o  mcc2.cpp.js

@rem ******* DISABLED OPTIONS: **********
@rem         -s ASSERTIONS=1   ^
@rem      -s ALLOW_MEMORY_GROWTH=1  ^   # This makes it 3X slower!
@rem     -s TOTAL_MEMORY=16777216
