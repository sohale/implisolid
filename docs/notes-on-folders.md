Newer notes are above, followed by older notes.
Also see [js_iteration_2/history.md](js_iteration_2/history.md)

## More technical debt


void flush_geometry_queue(std::ostream& cout, int& normals_start, std::vector<REAL> &normals,  std::vector<REAL> &verts3, std::vector<vertexindex_type> &faces3, MarchingCubes<false,false>::e3map_t &e3map, int& next_unique_vect_counter)
void flush_geometry_queue(std::ostream& cout, int& normals_start, std::vector<REAL> &normals,  std::vector<REAL> &verts3, std::vector<vertexindex_type> &faces3, MarchingCubes<false,false>::e3map_t &e3map, int& next_unique_vect_counter)
.

  this.update_normals__()  deprecated in  js_iteration_2/geometry79.js
      remove the function.
  simply kill file
    src/cpp/implisolid_js_service/old_experimental.hpp

  js_iteration_2/implisolid_main_old.js
    versus ?

  .update_geometry_()
  versus
  .update_geometry1()

  What is `_build_geometry_u`? Should I include it?

  Just to know:
    * __set_range_of_used_faces()
    * __set_needsUpdate_flag()


    comparing _MS to their counterparts:
      by checkout to history at the time of their latest commits

        ./js_iteration_1/debug_methods_MS.hpp
        ./js_iteration_1/legacy/build_mcc2_MS.sh
        ./js_iteration_1/legacy/geometry77_MS.js
        ./js_iteration_1/legacy/mcc2_MS_3js.html
        ./js_iteration_1/legacy/mcc2_MS.cpp
        ./js_iteration_1/legacy/mcc2_marching_cubes_MS.hpp

    gemomemetry.js:
      Keep greometry77.js though
      Keep greometry73.js though
      threejs_r71 is already included in geometry79.js
      geometry79.js: the main one.

      Full list:
          ./js_iteration_1/legacy/geometry77_MS.js
          ./js_iteration_1/legacy/geometry77.js
          ./docs/implisolid-build/demo1/js/geometry79.js
          ./js_iteration_2/geometry79.js
          ./js_iteration_2/js/geometry73.js


  Technicala debt/understanding:
  This note:
       /* deprecaed. Use ImpliSolid.update_geometry() instead.
       i.e. swap geometry & implicit_service:
       geom.update_geometry(IMPLISOLID, true)  ->  IMPLISOLID.update_geometry(geom, true)
       */
      this.update_geometry_ = function(implicit_service, ignoreDefaultNormals) {

  Move all scripts like js_iteration_1/build_mcc2.sh to one folder at least. (older ones to its ./legacy?)
Technical debt:
        ```
        bool check_state() {
            return true;  // fixme. NOT WHEN UPDATING
        ```
Hierarchies/levels:
      *   .update_geometry_()   versus  .update_geometry1()
      * IMPLICIT: layers
      * OIMPLICIT_WORKER: layers
      * ...
Aspirations:
      * globalbox multiple boxes
      * multiplle strides


#### Automated tests:
    Write some unit tests. for js-compiled code. Tesds are in Javasfcript.
      (separate from C++ tests)
    Using nodejs! yes!
    Call those baasic functions first.
    The (implicit) objects are built on top of the basic operations.

    Write automated tests for the *compiled* file.
        Compiled versions:
          1. non-error, positive only: module of cpp: without .. objects
          2. objects
          3. queries: set notrmal, etc (also mesh operations?)
          4. mesh function. (subdivide, setnorm, getnorm, etc)

      Based on the compilation, a subset of them is activted?
            (unfortunately we dont use linkers)
            separate tests for each block?
            on-the-fly compilation (flag! the flag deletes the cache (js file), and re-compiles)
      Browser-level tests?
            * (Slenium?)
            * (How does threejs tests?) -- interestingly, threejs uses node. -- https://github.com/mrdoob/three.js/tree/dev/test/unit
              (also https://github.com/mrdoob/three.js/commit/1c685e6eb56e54c4609b7b3af764e2e5110abfbf )

      test multiple combinaation of subsets: using #define flags and recompiling.

      Only write tests that are ...
      Re-run scripts: anoother configuraiton with the same build script. (but without lengthy cloning and downloding)
        how to have different options?
        do the clone-... once, nd then repeat other e2e tests:using another e2e script. It first checks the clone & library download is done probperly. Then rune build. Then runs multiple test sets. (a test script?)


#### Ok I know
  Disabale worker
  * check createMyModule
  * check these exceptions throw FS.genericErrors[44];
  * it breaks inside the following trace (stack):
          init_scene
              geometry = IMPLICIT.getLiveGeometry_from_json(shape_json, polygonization_json, true, geom_callback__1, ignore_root_matrix);
          return service2.make_geometry(mp5_str, polygonization_setttings_json_str,
          impli2.make_geometry = function (mp5_str, polygonization_params_str, geometry_callback, allocate_buffer) {
          impli1.build_geometry(mp5_str, polygonization_params_str);
          which is in C++

  *

## Technical debt: Due refactor

Todo:

Tiger's contributions are in his branch `tiger-master`.
Those need to be moved to the main branch `revival-sohale`.
However, the main branch is moving away.
Hence, the contributions will not be automtically mergable.
It is ineavitable since keeping the old structure is doing more harm and slowing down ("delaying") the progress.
I made this decision on "Mon 2 May 2022 16:54 BST".
Hence, those are of less priority.
Let those be taken care of at a later time.

Soon `revival-sohale` branch will be heavily refactored and overhauled.

Tiger had some contributions which are left behind (burried for now). I leave this note to keep the potential of dreviving them, and mark this as a technical debt,
But Tiger's contributions disabled some important features such as incremental and web Workers (service), and some refinement to the main algorithm (and the `asmjscb` primitive).

So, in some future time,
to revive some of Tiger's contributions, use `diff` between relevant commits.
Then manually add them.

Todo:
1. extract the relevants commits' sha-numbers.
2. git diff sha1 sha2
3. Compare to sha3 (beginning of divergence)
4. Identify and list those functions. Document them.
5. Add each that is considered useful and of priority, one at a time.

### Unrelated note
Note: Also Tiger has done some work on strictly SDF objects. See branch `researchOnSDF`.

I am not sure about two other branches `assert_fixing` and `optimization_task`.

## C++ demos, tests

* `js_iteration_1/mcc2_3js_r79.html`
* worker+incremental demo: ? (threejs)
* `js_iteration_2/examples/mp5interactive/mp5_json_code.html` -- demo1: interactive demo (ace + python3 http.server + threejs)

* reusable module: IMPLICIT_WORKER,

* npm: (not implemented)

* main files:
* `js_iteration_1/build_mcc2.sh`
* python binding (C++)

### Multiple fine-tunings:
*
### build scripts:
* `demo1-*.sh`: e2e (fresh clone), build, run mp5 interactive demo on python
* js_iteration_1/build_g++.sh
* js_iteration_1/build_mcc2.sh
* js_iteration_1/build_mcc2_docker.sh

Tiger:
* `js_iteration_2/examples/function_editor_demos.js`

## Branches:
*  `revival-sohale` -- main one: worker+incremental demo,...
*  `tiger-master` --
*  `assert_fixing` --
*  `optimization_task` --
*  `researchOnSDF` --

## Python
### Main executables
* python_clean_code/ohtake_belyaev_demo_subdivision_projection_qem.py: demo_combination_plus_qem()
* python_implicit/vtk_mc.py
* ...
* experiementation/...

* python_implicit/ohtake_belyaev_2.py
* python_implicit/ohtake_belyaev_3.py
* python_implicit/ohtake_belyaev_4.py
* python_implicit/ohtake_belyaev_5.py

* python_implicit/ohtake_surface_projection.py
* python_implicit/ohtake_surface_projection_2.py
* python_implicit/ohtake_surface_projection_v2_5.py
Small ones:
* python_implicit/cylinder_example.py

* python_implicit/*.py

# Older notes:
## Folders and their purposes:
 ImpliSolid folders: (relative to `$IMPLISOLID_REPO/` )
   demos/demo1-*.sh
   demos/demo1       Target for local deployment
   demos/demo1/build ?
   /docs The github-pages for implisolid itself. (has its own implisolid-build)
   docs/implisolid-build
   docs/implisolid-build/demo1 - Target for github-pages deployment of implisolid (not for sohale.github.io)

   js_iteration_2/examples/     -  Location for original source file for demos (demo1 for now, and demo2, etc in future)
   js_iteration_2/examples/mp5interactive - (i.e. $DEMO0 above) the demo1
       For other source locations, see $JSI2 and $EX_JSLIB above.

   build/ --- ?


### Folders in sohale.github.io
   demos/implisolid-build Also used for publishing (from here to github.com)
   demos/implisolid-build/demo1 - Target for github-pages deployment of sohale.github.io (not for implisolid's own page, https://sohale.github.io/implisolid/)
   demos/implisolid-build/demo1/mp5_json_code.html - todo: Note that this is not yet activated. If this is activated, there won't be need for updating sohale.github.io each time.
               Your site is published at https://sohale.github.io/implisolid/
                        https://sohale.github.io/implisolid/implisolid-build
                                      mmaps to:      ./docs/implisolid-build
                        https://sohale.github.io/implisolid/implisolid-build/demo1/latest-commit-log.txt

### The implisolid-build repo:
      Folder structure of the implisolid-build repo (as submoodule):
       ./demo1
       ./opt    -- old usage for wedesign.live, the location for optimised build (CDN) using C++/Emmscripten



## Useful absolute folders:
   ~/cs/sohale.github.io/demos/implisolid-build/demo1
   ~/cs/mp5/implisolid/demos/



## Steps for actual publish: (not local)
 1. Change ln to cp (See `CP()` vs `CP_not()` )  (On this script)
 2. Change DEPLOY_LOCATION to $IMPLISOLID_REPO/ `docs/implisolid-build/demo1` (On this script)
 3. Change directory: `cd` $IMPLISOLID_REPO/ `demos`
 4. Run `bash demo1-deploy.sh`

 Go to the implisolid-build repo (as a submodule)
 `cd` to  $IMPLISOLID_REPO/  `docs/implisolid-build`
 5. git commit and push `implisolid-build` into repo (Make sure you add files if thery are not added)
    5.1 You may want to squash or remove content from history to save space, and do a `git push -f`, to avoiud bloadting of files.

 Back to this repo:
 6. Revert changes to this file (`CP` and `DEPLOY_LOCATION`).

 Go to repo sohale.github.io:
 7. cd sohale.github.io/demos/implisolid-build
 8. pull changes from github to there
 9. cd .. to sohale.github.io itself
 8. git commit and push the main github.io page (i.e. sohale.github.com ) for github.com to pick up the changes from `implisolid-build`
    Check: https://github.com/sohale/sohale.github.io/actions for build pipeline (uit may publish despite errors here)
    Check https://github.com/sohale/sohale.github.io/settings for errors that fail the actual publishing.
 Test on browser:
    Visit: http://sohale.github.io/demos/implisolid-build/demo1/mp5_json_code.html
    visit http://sohale.github.io/demos/implisolid-build/demo1/mp5_json_code.html
    Check http://sohale.github.io/demos/implisolid-build/demo1/latest-commit.txt
    Check http://sohale.github.io/demos/implisolid-build/demo1/latest-commit-log.txt
 compare: https://sohale.github.io/implisolid/implisolid-build/demo1/mp5_json_code.html
     and: https://sohale.github.io/implisolid/implisolid-build/demo1/latest-commit-log.txt
     and: https://sohale.github.io/implisolid/


 9. (Not part of deployment:) Don't forget to git commit and push current (implisolid) repo for your other changes.
 10. Implisolid has its own github pages. Update that too (instructions to be added here)

 This did not include rebuilding the implisolid iself. See `demo1-build.sh` for that.




## Steps for local publish for test: (on MacOS)
 1. Change ln to ln (See `CP()` vs `CP_not()` )  (On this script)
 2. Change DEPLOY_LOCATION to $IMPLISOLID_REPO/ `docs/implisolid-build/demo1` (On this script)
 3. Change directory: cd $IMPLISOLID_REPO/ `demos`
 4. Run `bash demo1-deploy.sh`
 5. if python http server shows erro, you maay want to kill the previous instance (it may be not required).
 6. ...
    ...  Check http://localhost:8000/latest-commit-log.txt
    ...  Check http://localhost:8000/latest-commit.txt
