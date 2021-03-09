
#!/bin/bash

# deploy: to put it "there"
# run:  bash deploy-demo-1.sh
# todo: put there in any given (parametrised) location
set -e

# target:
export DEMO_LOCATION=/Users/a9858770/cs/mp5/implisolid/demos/demo1

#sources:
export MP5_PRIVATE=/Users/a9858770/cs/mp5/mp5-private
export IMPLISOLID=/Users/a9858770/cs/mp5/implisolid
export DEMO0=$IMPLISOLID/js_iteration_2/examples/mp5interactive
export JSI2=$IMPLISOLID/js_iteration_2
export JS_EX1=$IMPLISOLID/js_iteration_2/examples/js

# haldbuild: as built on github
#   $BUILT/opt/mcc2.compiled.js
#export BUILT=$IMPLISOLID/docs/implisolid-build
# $BUILT/opt/mcc2.compiled.js
#export compiled_file=$BUILT/opt/mcc2.compiled.js
# softbuild: loocally just built
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


# /Users/a9858770/cs/mp5/mp5-private/implisolid/js_iteration_1/controls/OrbitControls_r79.js

#ln -s $MP5_PRIVATE/implisolid/js_iteration_1/controls/OrbitControls_r79.js $DEMO_LOCATION/js/
ln -s $JS_EX1/OrbitControls_r79-copy.js $DEMO_LOCATION/js/OrbitControls_r79.js

ln -s $DEMO0/mp5_json_code.html $DEMO_LOCATION/
# ln -s $DEMO0/2222.html $DEMO_LOCATION/
#ls $JSI2/js
#ls -l $JSI2
ln -s $JSI2/geometry79.js $DEMO_LOCATION/js/
ln -s $JSI2/implisolid_main.js $DEMO_LOCATION/js/
ln -s $JSI2/js/js_utils.js $DEMO_LOCATION/js/
ln -s $JSI2/js/pointset_utils.js $DEMO_LOCATION/js/
ln -s $JSI2/js/arrow_utils.js $DEMO_LOCATION/js/
ln -s $JS_EX1/example_objects.js $DEMO_LOCATION/js/
ln -s $JS_EX1/example_materials.js $DEMO_LOCATION/js/
ln -s $JS_EX1/performance_graphs.js $DEMO_LOCATION/js/
ln -s $JS_EX1/misc_props.js $DEMO_LOCATION/js/
ln -s $JS_EX1/boundingbox_utils.js $DEMO_LOCATION/js/
ln -s $compiled_file $DEMO_LOCATION/js/
# ln -s $BUILT/opt/mcc2.compiled.js.mem $DEMO_LOCATION/js/
ln -s $JS_EX1/simple_assert.js $DEMO_LOCATION/js/

# implisolid/js_iteration_2/examples/js/simple_assert.js


# ls -l $DEMO_LOCATION
# examine and produce all errors (the "ln -s" file links that the file is non-existant )
#ls -1 $DEMO_LOCATION | xargs cat 1>/dev/null
echo "Checking errors *******"
pushd $DEMO_LOCATION/js
ls -1 . | xargs cat 1>/dev/null
popd

pwd
echo 'fine'

python3 -m http.server 8000 &
export server_pid=$!
echo $server_pid >server_pid-$server_pid.pid

# warning: MacOS-specific code
#open http://localhost:8000/
open -a "Google Chrome" http://localhost:8000/

echo "The current server PID is:"
#ps aux|grep -ie python
ps aux|grep -ie python|grep http
sleep 1
echo "kill $server_pid"
echo "\n\n\n\n ****************"
