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
