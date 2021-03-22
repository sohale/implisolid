
#!/bin/bash

# runs demo1 locally for MacOS

# deploy: to put it "there"
# run:  bash deploy-demo-1.sh
# todo: put there in any given (parametrised) location
set -e

# target:
# local deploy
export DEMO_LOCATION=/Users/a9858770/cs/mp5/implisolid/demos/demo1
# public deploy (to github-pages)
export DEMO_LOCATION=/Users/a9858770/cs/mp5/implisolid/docs/implisolid-build/demo1

#sources:
export MP5_PRIVATE=/Users/a9858770/cs/mp5/mp5-private
export IMPLISOLID=/Users/a9858770/cs/mp5/implisolid
export DEMO0=$IMPLISOLID/js_iteration_2/examples/mp5interactive
export JSI2=$IMPLISOLID/js_iteration_2
export JS_EX1=$IMPLISOLID/js_iteration_2/examples/js

# Folders and their purposes:
# ImpliSolid folders: (relative to `$IMPLISOLID/` )
#   demos/demo1-*.sh
#   demos/demo1       Target for local deployment
#   demos/demo1/build ?
#   /docs The github-paages for implisolid itself. (has its own implisolid-build)
#   docs/implisolid-build
#   docs/implisolid-build/demo1 - Target for github-pages deployment of implisolid (not for sohale.github.io)
#
#   js_iteration_2/examples/     -  Location for original source file for demos (demo1 for now, and demo2, etc in future)
#   js_iteration_2/examples/mp5interactive - (i.e. $DEMO0 above) the demo1
#       For other source locations, see $JSI2 and $JS_EX1 above.
#
# Folders in sohale.github.io
#   demos/implisolid-build Also used for publishing (from here to github.com)
#   demos/implisolid-build/demo1 - Target for github-pages deployment of sohale.github.io (not for implisolid's own page, https://sohale.github.io/implisolid/)
#   demos/implisolid-build/demo1/mp5_json_code.html - todo: Note that this is not yet activated. If this is activated, there won't be need for updating sohale.github.io each time.
#
# The implisolid-build repo:
#      Folder structure of the implisolid-build repo (as submoodule):
#       ./demo1
#       ./opt    -- old usage for wedesign.live, the location for optimised build (CDN) using C++/Emmscripten


# Steps for actual publish: (not local)
# 1. Change ln to cp (See `CP()` vs `CP_not()` )  (On this script)
# 2. Change DEMO_LOCATION to $IMPLISOLID/ `docs/implisolid-build/demo1` (On this script)
# 3. Change directory: `cd` $IMPLISOLID/ `demos`
# 4. Run `bash demo1-deploy.sh`
#
# Go to the implisolid-build repo (as a submodule)
# `cd` to  $IMPLISOLID/  `docs/implisolid-build`
# 5. git commit and push `implisolid-build` into repo (Make sure you add files if thery are not added)
#    5.1 You may want to squash or remove content from history to save space, and do a `git push -f`, to avoiud bloadting of files.
#
# Back to this repo:
# 6. Revert changes to this file (`CP` and `DEMO_LOCATION`).
#
# Go to repo sohale.github.io:
# 7. cd sohale.github.io/demos/implisolid-build
# 8. pull changes from github to there
# 9. cd .. to sohale.github.io itself
# 8. git commit and push the main github.io page (i.e. sohale.github.com ) for github.com to pick up the changes from `implisolid-build`
#    Check: https://github.com/sohale/sohale.github.io/actions for build pipeline (uit may publish despite errors here)
#    Check https://github.com/sohale/sohale.github.io/settings for errors that fail the actual publishing.
# Test on browser:
#    Visit: http://sohale.github.io/demos/implisolid-build/demo1/mp5_json_code.html
#    visit http://sohale.github.io/demos/implisolid-build/demo1/mp5_json_code.html
#    Check http://sohale.github.io/demos/implisolid-build/demo1/latest-commit.txt ?
#    Check http://sohale.github.io/demos/implisolid-build/demo1/latest-commit-log.txt ?
#
# 9. (Not part of deployment:) Don't forget to git commit and push current (implisolid) repo for your other changes.
# 10. Implisolid has its own github pages. Update that too (instructions to be added here)
#
# This did not include rebuilding the implisolid iself. See `demo1-build.sh` for that.


# Steps for local publish for test: (on MacOS)
# 1. Change ln to ln (See `CP()` vs `CP_not()` )  (On this script)
# 2. Change DEMO_LOCATION to $IMPLISOLID/ `docs/implisolid-build/demo1` (On this script)
# 3. Change directory: cd $IMPLISOLID/ `demos`
# 4. Run `bash demo1-deploy.sh`
# 5. if python http server shows erro, you maay want to kill the previous instance (it may be not required).
# 6. ...
#    ...  Check http://localhost:8000/latest-commit-log.txt
#    ...  Check http://localhost:8000/latest-commit.txt


# haldbuild: as built on github
#   $BUILT/opt/mcc2.compiled.js
#export BUILT=$IMPLISOLID/docs/implisolid-build
# $BUILT/opt/mcc2.compiled.js
#export compiled_file=$BUILT/opt/mcc2.compiled.js
# softbuild: locally just built
export BUILD_LOCATION=/Users/$USER/cs/mp5/implisolid/demos/build
# $BUILD_LOCATION/mcc2.compiled.js
export compiled_file=$BUILD_LOCATION/mcc2.compiled.js

# prepare

# first time: (`clone`s sub-modules)
# git submodule update --init --recursive
# ?
# git submodule update --recursive --remote

# not first time:
# git pull --recurse-submodules

# Look at the branches !
# * [new branch]      assert_fixing     -> origin/assert_fixing
# * [new branch]      optimization_task -> origin/optimization_task
# * [new branch]      researchOnSDF     -> origin/researchOnSDF

# fail fast
#set -x; mkdir demo1
mkdir -p demo1
cd demo1
pwd | grep -qE "/demo1$"
touch deleteme.js
# pwd | grep -qE "/demo1$"; echo $?
set -x; pwd | grep -qE "/demo1$"
# never `rm` unless you are there for sure:
rm -rv ../demo1/*
cd ..
rmdir ./demo1

# Somehow assert empty ?

# fail-fast
mkdir demo1
mkdir demo1/js

cd demo1
pwd | grep -qE "/demo1$"

# local run:
CP_not() {
  ln -s $1 $2
}

# server deploy (then push to the "implicit_build" repo)
CP() {
  cp   $1 $2
}

# /Users/a9858770/cs/mp5/mp5-private/implisolid/js_iteration_1/controls/OrbitControls_r79.js

#CP $MP5_PRIVATE/implisolid/js_iteration_1/controls/OrbitControls_r79.js $DEMO_LOCATION/js/
CP $JS_EX1/OrbitControls_r79-copy.js $DEMO_LOCATION/js/OrbitControls_r79.js

CP $DEMO0/mp5_json_code.html $DEMO_LOCATION/
# CP $DEMO0/2222.html $DEMO_LOCATION/
#ls $JSI2/js
#ls -l $JSI2
CP $JSI2/geometry79.js $DEMO_LOCATION/js/
CP $JSI2/implisolid_main.js $DEMO_LOCATION/js/
CP $JSI2/js/js_utils.js $DEMO_LOCATION/js/
CP $JSI2/js/pointset_utils.js $DEMO_LOCATION/js/
CP $JSI2/js/arrow_utils.js $DEMO_LOCATION/js/
CP $JS_EX1/example_objects.js $DEMO_LOCATION/js/
CP $JS_EX1/example_materials.js $DEMO_LOCATION/js/
CP $JS_EX1/performance_graphs.js $DEMO_LOCATION/js/
CP $JS_EX1/misc_props.js $DEMO_LOCATION/js/
CP $JS_EX1/boundingbox_utils.js $DEMO_LOCATION/js/
CP $compiled_file $DEMO_LOCATION/js/
# CP $BUILT/opt/mcc2.compiled.js.mem $DEMO_LOCATION/js/
CP $JS_EX1/simple_assert.js $DEMO_LOCATION/js/

# self-refelection of deploy (version inspect endpoint!)
# implisolid/js_iteration_2/examples/js/simple_assert.js
git rev-parse HEAD >$DEMO_LOCATION/latest-commit.txt
git log  -n 1 >$DEMO_LOCATION/latest-commit-log.txt
echo "\n\ngit diff\n" >>latest-commit-log.txt
git diff >>latest-commit-log.txt
date >>latest-commit-log.txt
echo "for latest commit info click: try http://localhost:8000/latest-commit-log.txt"

# ls -l $DEMO_LOCATION
# examine and produce all errors (the "ln -s" file links that the file is non-existant )
#ls -1 $DEMO_LOCATION | xargs cat 1>/dev/null
echo "Checking errors *******"
pushd $DEMO_LOCATION/js
ls -1 . | xargs cat 1>/dev/null
cd ..

pwd
echo 'fine'

# for local run:
echo "" >$DEMO_LOCATION/run.sh
echo "echo visit http://localhost:8000/mp5_json_code.html" >>$DEMO_LOCATION/run.sh
echo "python3 -m http.server 8000" >>$DEMO_LOCATION/run.sh

python3 -m http.server 8000 &
export server_pid=$!
echo $server_pid >server_pid-$server_pid.pid

popd
pwd

# warning: MacOS-specific code
echo "click on mp5_json_code.html @"
# open -a "Google Chrome" http://localhost:8000/
open -a "Google Chrome" http://localhost:8000/mp5_json_code.html


echo "The current server PID is:"
#ps aux|grep -ie python
ps aux|grep -ie python|grep http
sleep 1
echo "kill $server_pid"
echo "\n\n\n\n ****************"

echo "python processes to kill $(ps aux|grep -ie python|grep http|cut -c 17-25 | xargs echo)"
