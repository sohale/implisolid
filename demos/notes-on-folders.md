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
* clean_code/ohtake_belyaev_demo_subdivision_projection_qem.py: demo_combination_plus_qem()
* implicit/vtk_mc.py
* ...
* experiementation/...

* implicit/ohtake_belyaev_2.py
* implicit/ohtake_belyaev_3.py
* implicit/ohtake_belyaev_4.py
* implicit/ohtake_belyaev_5.py

* implicit/ohtake_surface_projection.py
* implicit/ohtake_surface_projection_2.py
* implicit/ohtake_surface_projection_v2_5.py
Small ones:
* implicit/cylinder_example.py

* implicit/*.py

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
       For other source locations, see $JSI2 and $JS_EX1 above.

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
