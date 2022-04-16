#!/bin/bash

# runs demo1 locally for MacOS
set -ex
function assert_env_nonempty() {
  if [ ".$1" = "." ]; then
    echo "shell env is empty"; echo $2
    return 1
  fi
}

# args:
assert_env_nonempty $BASELOC1 "env-argument BASELOC1= not specified"
assert_env_nonempty $SCRIPTS_DIR "env-argument SCRIPTS_DIR= not specified"
#export SCRIPTS_DIR=$REPO_ROOT/demos
#assert_env_nonempty $REPO_ROOT
# nothing more
# arg: pwd (used implicitly via REPO_ROOT) # no longer!

# not used: BASELOC2, BASELOC3

# deploy: to put it "there"
# run:  BASELOC1=$(pwd) bash deploy-demo-1.sh
# todo: put there in any given (parametrised) location

# Prerequisite:
# * compiles files to be ready here: $IMPLISOLID_REPO/demos/build. Use demos/demo1-build.sh
# * three repos are cloned

set -e

# implicitly, the arg is pwd

# REPO_ROOT=$(git rev-parse --show-toplevel)
# only for this?
#export SCRIPTS_DIR=$REPO_ROOT/demos
# Incorrect: It is hard coded.
#source $SCRIPTS_DIR/base-locations.sh

# rename: scripts
source $SCRIPTS_DIR/utils.sh


# Altenative names: BASE_IMPLISOLID, REPOBASE_IMPLISOLID, IMPLISOLID, REPO_IMPLISOLID, IMPLISOLID_REPO
# repo bases
export IMPLISOLID_REPO=$BASELOC1/implisolid
# same as $REPO_ROOT

#export SOHALE_IO_REPO=$BASELOC3/sohale.github.io/
# export BASE_MP5_PRIVATE=$BASELOC2/mp5/mp5-private
# relative location of repos
export IMPLISOLID_BUILD_REPO=$IMPLISOLID_REPO/docs/implisolid-build
# todo: specify a separate independent folder  => you wil know the structure BEFORE this.
    # qt the tie oof e2e sscript

# based on configuration


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

# see notes-on-folders.md



# hard-built: as built on github
#export BUILT=$IMPLISOLID_BUILD_REPO

# soft-built: locally just built
export BUILD_LOCATION=$IMPLISOLID_REPO/demos/build
export compiled_file=$BUILD_LOCATION/mcc2.compiled.js
# Pre-built from submodule on github:
#export compiled_file=$BUILT/opt/mcc2.compiled.js

cd $START_DEMOS

rm -rv demo1; mkdir demo1; cd demo1

function gather_files_for_deploy() {
    # Gathers various files and prepares them for deploying.
    # Why not move them permanently on the (git) repo?
    # Copies file from 4 sources:
    #    JS_EX1,JSI2,DEMO0, compiled_file
    # To: destination: DEPLOY_LOCATION

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

    # OrbitControls_r79.js was also in `mp5-private/implisolid/js_iteration_1/controls/OrbitControls_r79.js`
    CP $JS_EX1/OrbitControls_r79-copy.js $DEPLOY_LOCATION/js/OrbitControls_r79.js

    CP $DEMO0/mp5_json_code.html $DEPLOY_LOCATION/
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
    CP $JS_EX1/simple_assert.js $DEPLOY_LOCATION/js/

    echo "compiled_file $compiled_file"
    CP $compiled_file $DEPLOY_LOCATION/js/
    # CP $BUILT/opt/mcc2.compiled.js.mem $DEPLOY_LOCATION/js/

    # self-refelection of deploy (version inspect endpoint!)
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

cat >$DEPLOY_LOCATION/run.sh << LAUNCH_SCRIPT

echo visit http://localhost:8000/mp5_json_code.html
python3 -m http.server 8000
LAUNCH_SCRIPT

cat $DEPLOY_LOCATION/run.sh

echo "going to actually run"
pwd

# $IMPLISOLID_BUILD_REPO/demo1/.. # wrong?
cd $DEPLOY_LOCATION
pwd
DEPLOY_LOCATION="$DEPLOY_LOCATION" bash $SCRIPTS_DIR/demo1-run-local.sh
