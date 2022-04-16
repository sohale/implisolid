# use with `source` only
# todo: use groovy

function assert_env_nonempty() {
  if [ ".$1" = "." ]; then
    echo "shell env is empty"; echo $2
    return 1
  fi
}

# args:
assert_env_nonempty $IMPLISOLID "mising IMPLISOLID="
# Seed for parameter values:
#REPO_ROOT=$(git rev-parse --show-toplevel)
#source $REPO_ROOT/demos/base-locations.sh

# Parameters,
# Parameter values from the specific configuration (relative locaation of build, lib, etc):
# target:
export DEMO_LOCATION=$IMPLISOLID/demos/demo1
export BUILD_LOCATION=$IMPLISOLID/demos/build
export LIB_FOLDER=$BUILD_LOCATION/lib

printf "DEMO_LOCATION:$DEMO_LOCATION, \nBUILD_LOCATION:$BUILD_LOCATION, \nLIB_FOLDER:$LIB_FOLDER\n"


# configuration-specific

# Other old files:
# $IMPLISOLID_REPO/docs/implisolid-build/demo1/js
# $BUILT/opt/mcc2.compiled.js
# $BUILT/opt/mcc2.compiled.js.mem

# target:
# local deploy
#export DEMO_LOCATION=$IMPLISOLID_REPO/demos/demo1
# remote/public deploy (to github-pages)
#export DEMO_LOCATION=$IMPLISOLID_BUILD_REPO/demo1


# probably incorrect
#IMPLISOLID=/home/$USER/mp5-private/implisolid
