
Also see [demos/notes-on-folders.md](demos/notes-on-folders.md)
## 6 May 2022

  TyH's changes identified:
    * exploded produce_mesh into its usage. Why?
    * extrusion.hpp (by Blobixx) is removed by TyH.
        but some dependencies are not removed
    * Why did he rename convex_polygon to concave_polygon?
    * Made changes to implicit_function/screw.hpp , but why those arn't there? implicit_function/screw.hpp
    * Why removed `isFinal` from implisolid_wroker.js ? (default: true)
    * Why removed `is_progress_update` also. (Default: false. Code killed)
    * Why removed `wwapi.send_progress_update` ???
    * Why so much code removed from above `class MarchingCubes {`?

  I am going through this
    `git diff 69ce747c7bbf4959d79a16298187c0af2ecf8c5e..99dc78cf6242556383db31c54a9a576836849c56`
    Reached page 5 (reached 16) (out of 41 pages)
    TBC

  Random questions:
      what was js_iteration_2/examples/js/example_objects.js ?
      What is the make_fg() ? Why some  asm.js for random object is there?
      What is `_build_geometry_u`?

  I found out:
      * -Wno-dollar-in-identifier-extension  Seems too be explained here https://github.com/emscripten-core/emscripten/issues/7113
  Currently, the only error that remains today is:
      `Uncaught (in promise)`
## 9 March 2021

* produce_object_old2() is no longer used. IT was exported but the function wass already dremoved from mcc2.cpp

* Updated to latest Emscripten, Eigen,m and Boost:
   * removing `memoryInitializerPrefixURL` due to the following error:
      `Assertion failed: Module.memoryInitializerPrefixURL option was removed, use Module.locateFile instead`

   * 'cwrap' was not exported. add it to EXTRA_EXPORTED_RUNTIME_METHODS (see the FAQ)

* moved to branch: revival-sohale from commit
* reviving some code destoryed by others.
* will have to check commits; master, my last change ( 69ce747c7bbf4959d79a16298187c0af2ecf8c5e ) and Tiger's  (proobably) last change (  99dc78cf6242556383db31c54a9a576836849c56 ):
  `git diff 69ce747c7bbf4959d79a16298187c0af2ecf8c5e..99dc78cf6242556383db31c54a9a576836849c56`

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

* Check who added `_build_geometry_u`
* Some good links:
    * My good version: https://github.com/sohale/implisolid/blob/69ce747c7bbf4959d79a16298187c0af2ecf8c5e/js_iteration_2/implicit_function/javascript_implicit_function.hpp
* Possible plan:
  1. Compare all changes and revert all, except for ftures added.
        * until all changes  (diff) from that old version of mine are minimal? yes.

  2. Or alternatively, continue all from an old commit? NO.
    * starting from there (my old chabges), then I can apply all of these later changes that I made today, as a separate diff.
### Questions
* What was the point about `MORE_ABOUT_INFO`?
    * It is shown on console
* A different way of pulling the libraries `/lib`
* What is `mcc2_MS.cpp` ?

### reminder:
* `cwrap` is used for calling a function inside?
* `BUILD_AS_WORKER` is defined by me.
