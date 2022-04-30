# use with `source` only
# todo: use groovy

function assert_env_nonempty() {
  if [ ".$1" = "." ]; then
    echo "shell env is empty"; echo $2
    return 1
  fi
}

# Parameters,
# Parameter values from the specific configuration (relative location of build, lib, etc):

# args:
assert_env_nonempty $IMPLISOLID "mising IMPLISOLID="

# target:
export BUILD_LOCATION=$IMPLISOLID/demos/build
export LIB_FOLDER=$BUILD_LOCATION/lib
printf "BUILD_LOCATION:$BUILD_LOCATION, \nLIB_FOLDER:$LIB_FOLDER\n"


# BUILD_LOCATION is not just the TARGE, but alsoo conotains the two libraries (eigen, boost)

# configuration-specific

# Other old files:
# $IMPLISOLID_REPO/docs/implisolid-build/demo1/js
# $BUILT/opt/mcc2.compiled.js
# $BUILT/opt/mcc2.compiled.js.mem


