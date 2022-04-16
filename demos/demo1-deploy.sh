#!/bin/bash

# runs demo1 locally for MacOS

# deploy: to put it "there"
# run:  bash deploy-demo-1.sh
# todo: put there in any given (parametrised) location

# Prerequisite:
# * compiles files to be ready here: $IMPLISOLID_REPO/demos/build. Use demos/demo1-build.sh
# * three repos are cloned

set -e

REPO_ROOT=$(git rev-parse --show-toplevel)
# only for this?
export SCRIPTS_DIR=$REPO_ROOT/demos
source $SCRIPTS_DIR/base-locations.sh

# Altenative names: BASE_IMPLISOLID, REPOBASE_IMPLISOLID, IMPLISOLID, REPO_IMPLISOLID, IMPLISOLID_REPO
# repo bases
export IMPLISOLID_REPO=$BASELOC1/implisolid
export SOHALE_IO_REPO=$BASELOC3/sohale.github.io/
# export BASE_MP5_PRIVATE=$BASELOC2/mp5/mp5-private
# relative location of repos
export IMPLISOLID_BUILD_REPO=$IMPLISOLID_REPO/docs/implisolid-build



# target:
# local deploy
#export DEMO_LOCATION=$IMPLISOLID_REPO/demos/demo1
# remote/public deploy (to github-pages)
#export DEMO_LOCATION=$IMPLISOLID_BUILD_REPO/demo1

#export LOCAL_DEPLOY_LOCATION=$IMPLISOLID_REPO/demos/demo1
export REMOTE_DEPLOY_LOCATION=$IMPLISOLID_BUILD_REPO/demo1
# Two alternatives:
# why two options?
#export DEPLOY_LOCATION=$LOCAL_DEPLOY_LOCATION
export DEPLOY_LOCATION=$REMOTE_DEPLOY_LOCATION
# simplified:
export DEPLOY_LOCATION=$IMPLISOLID_BUILD_REPO/demo1


#sources deployed files:
export DEMO0=$IMPLISOLID_REPO/js_iteration_2/examples/mp5interactive
export JSI2=$IMPLISOLID_REPO/js_iteration_2
export JS_EX1=$IMPLISOLID_REPO/js_iteration_2/examples/js
ls  $IMPLISOLID_REPO/js_iteration_2/examples/js/OrbitControls_r79-copy.js
ls  $JS_EX1/OrbitControls_r79-copy.js


export START_DEMOS=$IMPLISOLID_REPO/demos

# Folders and their purposes:
# ImpliSolid folders: (relative to `$IMPLISOLID_REPO/` )
#   demos/demo1-*.sh
#   demos/demo1       Target for local deployment
#   demos/demo1/build ?
#   /docs The github-pages for implisolid itself. (has its own implisolid-build)
#   docs/implisolid-build
#   docs/implisolid-build/demo1 - Target for github-pages deployment of implisolid (not for sohale.github.io)
#
#   js_iteration_2/examples/     -  Location for original source file for demos (demo1 for now, and demo2, etc in future)
#   js_iteration_2/examples/mp5interactive - (i.e. $DEMO0 above) the demo1
#       For other source locations, see $JSI2 and $JS_EX1 above.
#
#   build/ --- ?
#
#
# Folders in sohale.github.io
#   demos/implisolid-build Also used for publishing (from here to github.com)
#   demos/implisolid-build/demo1 - Target for github-pages deployment of sohale.github.io (not for implisolid's own page, https://sohale.github.io/implisolid/)
#   demos/implisolid-build/demo1/mp5_json_code.html - todo: Note that this is not yet activated. If this is activated, there won't be need for updating sohale.github.io each time.
#               Your site is published at https://sohale.github.io/implisolid/
#                        https://sohale.github.io/implisolid/implisolid-build
#                                      mmaps to:      ./docs/implisolid-build
#                        https://sohale.github.io/implisolid/implisolid-build/demo1/latest-commit-log.txt
#
# The implisolid-build repo:
#      Folder structure of the implisolid-build repo (as submoodule):
#       ./demo1
#       ./opt    -- old usage for wedesign.live, the location for optimised build (CDN) using C++/Emmscripten

# Useful absolute folders:
#   ~/cs/sohale.github.io/demos/implisolid-build/demo1
#   ~/cs/mp5/implisolid/demos/

# Steps for actual publish: (not local)
# 1. Change ln to cp (See `CP()` vs `CP_not()` )  (On this script)
# 2. Change DEPLOY_LOCATION to $IMPLISOLID_REPO/ `docs/implisolid-build/demo1` (On this script)
# 3. Change directory: `cd` $IMPLISOLID_REPO/ `demos`
# 4. Run `bash demo1-deploy.sh`
#
# Go to the implisolid-build repo (as a submodule)
# `cd` to  $IMPLISOLID_REPO/  `docs/implisolid-build`
# 5. git commit and push `implisolid-build` into repo (Make sure you add files if thery are not added)
#    5.1 You may want to squash or remove content from history to save space, and do a `git push -f`, to avoiud bloadting of files.
#
# Back to this repo:
# 6. Revert changes to this file (`CP` and `DEPLOY_LOCATION`).
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
#    Check http://sohale.github.io/demos/implisolid-build/demo1/latest-commit.txt
#    Check http://sohale.github.io/demos/implisolid-build/demo1/latest-commit-log.txt
# compare: https://sohale.github.io/implisolid/implisolid-build/demo1/mp5_json_code.html
#     and: https://sohale.github.io/implisolid/implisolid-build/demo1/latest-commit-log.txt
#     and: https://sohale.github.io/implisolid/
#
#
# 9. (Not part of deployment:) Don't forget to git commit and push current (implisolid) repo for your other changes.
# 10. Implisolid has its own github pages. Update that too (instructions to be added here)
#
# This did not include rebuilding the implisolid iself. See `demo1-build.sh` for that.


# Steps for local publish for test: (on MacOS)
# 1. Change ln to ln (See `CP()` vs `CP_not()` )  (On this script)
# 2. Change DEPLOY_LOCATION to $IMPLISOLID_REPO/ `docs/implisolid-build/demo1` (On this script)
# 3. Change directory: cd $IMPLISOLID_REPO/ `demos`
# 4. Run `bash demo1-deploy.sh`
# 5. if python http server shows erro, you maay want to kill the previous instance (it may be not required).
# 6. ...
#    ...  Check http://localhost:8000/latest-commit-log.txt
#    ...  Check http://localhost:8000/latest-commit.txt





# haldbuild: as built on github
#   $BUILT/opt/mcc2.compiled.js
#export BUILT=$IMPLISOLID_BUILD_REPO
# $BUILT/opt/mcc2.compiled.js
#export compiled_file=$BUILT/opt/mcc2.compiled.js
# softbuild: locally just built
export BUILD_LOCATION=$IMPLISOLID_REPO/demos/build
# $BUILD_LOCATION/mcc2.compiled.js
export compiled_file=$BUILD_LOCATION/mcc2.compiled.js

# prepare

# submodules


# fail fast
#set -x; mkdir demo1

# First time:
#mkdir -p $START_DEMOS
#mkdir -p $IMPLISOLID_REPO/docs/implisolid-build/demo1/js


function tricks () {
    # create a clean folder called 'demo1' and goes inside it in a safe way
    # demo1 is created inside $START_DEMOS
    # equivalent to "cd ./demo1"

    # all paths need to be absolute
    cd $START_DEMOS

    mkdir -p demo1
    # assert; make sure we are inside 'demo1', and it exists
    cd demo1
    pwd | grep -qE "/demo1$"
    touch deleteme.js
    # pwd | grep -qE "/demo1$"; echo $?
    set -x; pwd | grep -qE "/demo1$"
    # never `rm` unless you are there for sure:
    rm -rv ../demo1/*
    cd ..
    # safe delete of folder
    rmdir ./demo1

    # Somehow assert empty ?

    # fail-fast
    mkdir demo1
    mkdir demo1/js

    cd demo1
    # assert; make sure we are inside 'demo1', and it exists
    pwd | grep -qE "/demo1$"
}

tricks

function gather_files_for_deploy() {
    # Gathers various files and prepares them for deploying.
    # Why not move them permanently on the (git) repo?
    # sources: JS_EX1,JSI2,DEMO0, compiled_file, destination: DEPLOY_LOCATION

    #compiled_file should be ready ( in $BUILD_LOCATION )

    mkdir -p $DEPLOY_LOCATION/js

    # local run:
    CP_not() {
      ln -s $1 $2
    }

    # server deploy (then push to the "implicit_build" repo)
    CP() {
      cp   $1 $2
    }

    echo "JS_EX1: $JS_EX1"
    echo "DEPLOY_LOCATION: $DEPLOY_LOCATION"

    # $BASE_MP5_PRIVATE/implisolid/js_iteration_1/controls/OrbitControls_r79.js

    echo JS_EX1
    ls -alt $JS_EX1
    echo c1
    #CP $BASE_MP5_PRIVATE/implisolid/js_iteration_1/controls/OrbitControls_r79.js $DEPLOY_LOCATION/js/
    CP $JS_EX1/OrbitControls_r79-copy.js $DEPLOY_LOCATION/js/OrbitControls_r79.js
    echo c2

    CP $DEMO0/mp5_json_code.html $DEPLOY_LOCATION/
    #ls $JSI2/js
    #ls -l $JSI2
    CP $JSI2/geometry79.js $DEPLOY_LOCATION/js/
    CP $JSI2/implisolid_main.js $DEPLOY_LOCATION/js/
    CP $JSI2/js/js_utils.js $DEPLOY_LOCATION/js/
    CP $JSI2/js/pointset_utils.js $DEPLOY_LOCATION/js/
    CP $JSI2/js/arrow_utils.js $DEPLOY_LOCATION/js/
    CP $JS_EX1/example_objects.js $DEPLOY_LOCATION/js/
    CP $JS_EX1/example_materials.js $DEPLOY_LOCATION/js/
    CP $JS_EX1/performance_graphs.js $DEPLOY_LOCATION/js/
    CP $JS_EX1/misc_props.js $DEPLOY_LOCATION/js/
    CP $JS_EX1/boundingbox_utils.js $DEPLOY_LOCATION/js/
    CP $compiled_file $DEPLOY_LOCATION/js/
    # CP $BUILT/opt/mcc2.compiled.js.mem $DEPLOY_LOCATION/js/
    CP $JS_EX1/simple_assert.js $DEPLOY_LOCATION/js/

    # self-refelection of deploy (version inspect endpoint!)
    # implisolid/js_iteration_2/examples/js/simple_assert.js
}

################
#bash copy_files_for_deploy.sh
gather_files_for_deploy
################

function find_noprefix() {
  FOLDER=$1
  find $FOLDER |awk -v prefix="$FOLDER/" '{ match($0, prefix); printf substr($0, RSTART+RLENGTH); printf "\n"}'
}

find_noprefix $DEPLOY_LOCATION

# Info about latest commit to be available on browser
git rev-parse HEAD >$DEPLOY_LOCATION/latest-commit.txt
git log  -n 1 >$DEPLOY_LOCATION/latest-commit-log.txt
echo "\n\ngit diff\n" >>$DEPLOY_LOCATION/latest-commit-log.txt
git diff >>$DEPLOY_LOCATION/latest-commit-log.txt
date >>$DEPLOY_LOCATION/latest-commit-log.txt
echo "for latest commit info click: try http://localhost:8000/latest-commit-log.txt"

# ls -l $DEPLOY_LOCATION
# examine and produce all errors (the "ln -s" file links that the file is non-existant )
#ls -1 $DEPLOY_LOCATION | xargs cat 1>/dev/null
echo "Checking missing symbolic links *******"
# pushd $DEPLOY_LOCATION/js
cd $DEPLOY_LOCATION/js
ls -1 . | xargs cat 1>/dev/null
cd ..

cd $DEPLOY_LOCATION

pwd
echo 'fine'


echo 'for local run: (Also see ./demo1-run-local.sh)'

echo "" >$DEPLOY_LOCATION/run.sh
echo "echo visit http://localhost:8000/mp5_json_code.html" >>$DEPLOY_LOCATION/run.sh
echo "python3 -m http.server 8000" >>$DEPLOY_LOCATION/run.sh

cat $DEPLOY_LOCATION/run.sh

echo "going to actually run"
pwd

# $IMPLISOLID_BUILD_REPO/demo1/..
cd $DEPLOY_LOCATION
pwd
bash $SCRIPTS_DIR/demo1-run-local.sh cd $DEPLOY_LOCATION
