@echo "---------------------------------------------------"
@rem set EMCC_DEBUG=1

@rem python "C:\Program Files\Emscripten\emscripten\1.35.0\tools\webidl_binder.py" mcc1.idl  mcc1-glue
@rem echo Web-IDL done


@rem compiler option:           --post-js mcc1-glue.js  ^




em++ ^
        -I C:\sohail\March\emscripten\boost_1_61_0\   ^
        -s ASSERTIONS=1               ^
        -s DEMANGLE_SUPPORT=1   ^
        -s TOTAL_MEMORY=30100100    ^
                                    ^
        -s EXPORTED_FUNCTIONS="['_produce_object_old2', '_main', '_build_geometry', '_get_v_size', '_get_f_size', '_get_f', '_get_v', '_finish_geometry']" ^
        -s NO_EXIT_RUNTIME=1          ^
                                          ^
        -pedantic -std=c++14  ^
    mcc2.cpp  ^
        -o  mcc2.cpp.js




@rem        -s EXPORTED_FUNCTIONS="['_main']" ^
@rem        -s NO_EXIT_RUNTIME=0          ^
@rem      --profiling   ^

@rem      -s ALLOW_MEMORY_GROWTH=1  ^   # This makes it 3X slower!
@rem     -s TOTAL_MEMORY=16777216
