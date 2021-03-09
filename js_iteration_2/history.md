
## 9 March 2021

* produce_object_old2() is no longer used. IT was exported but the function wass already dremoved from mcc2.cpp

* removing `memoryInitializerPrefixURL` due to the following error:
`Assertion failed: Module.memoryInitializerPrefixURL option was removed, use Module.locateFile instead`

* mcc2.compiled.js:215903 'cwrap' was not exported. add it to EXTRA_EXPORTED_RUNTIME_METHODS (see the FAQ)


### intents:
* todo: mcc2.cpp should be moved to iteration2
* Remaining hpp in js_iteration1 are: ...
* Note about js_iteration1: Why some hpp files still remain there, should be move js, or cpp there?

* What is the preferred direcrtory structure for cpp versus js files?
  * We need to have the C++ parts standalone, if there is no js / Emscriptedn involved.
  * The interface parts need to be in a separate cpp file. mcc2.cpp, which will contain that, shoul dbe strcitly minimal.

* Make it use progressive/non-progressive/web-worker
* Use CMake
* separate 3 levels of API into separate files
* remove mcc1.cpp if no mention of it is used.
* split mcc2 into three parts: 1. residue of intrface to the MarchingCubes class, 2. very surface interface : only extern "C", 3. All the rest? (grand_algorithm - specific)
* Have two versions ready: with webworker and non-webworker.
* test suit for interface?
* Keep `polygonizer_algorithm_ob02` and other algorithms in the tests loop
* Remove all of those pointers
* RAND_MAX, make it long, or double, float, etc. cast inline.
* refactor long elseif of object_factory in (separation and conglumeration of concerns).

### Questions
* What was the point about `MORE_ABOUT_INFO`?
    * It is shown on console
* A different way of pulling the libraries `/lib`
* What is `mcc2_MS.cpp` ?

### reminder:
* `cwrap` is used for calling a function inside?
* `BUILD_AS_WORKER` is defined by me.
