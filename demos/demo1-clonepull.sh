
#!/bin/bash

# run:  bash deploy-demo-1.sh
set -e

# target:
export DEMO_LOCATION=/Users/a9858770/cs/mp5/implisolid/demos/demo1

#sources:
export MP5_PRIVATE=/Users/a9858770/cs/mp5/mp5-private
export IMPLISOLID=/Users/a9858770/cs/mp5/implisolid
export DEMO0=$IMPLISOLID/js_iteration_2/examples/mp5interactive
export JSI2=$IMPLISOLID/js_iteration_2
export JS_EX1=$IMPLISOLID/js_iteration_2/examples/js
export BUILT=$IMPLISOLID/docs/implisolid-build

# clones the repo and submodules
# pulls latest version
# pull latest dsocker
# ...

docker pull emscripten/emsdk

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
