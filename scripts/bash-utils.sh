# usage: SCRIPT_DIR=$( __current_script_dir_func )
set -u

function __current_script_dir_func () {
  cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd
}

function assert_env_nonempty() {
  #[[ $OSTYPE == 'darwin'* ]] || "Error: MacOS-specific code"
  #echo ">$1<"
  #echo [[ $1 == '' ]]
  #[[ $1 == '' ]] ; echo $?
  #[ -z "$1" ]; echo $?
  if [ ".$1" = "." ]; then
    # $STRING is empty
    echo "shell env is empty"
    echo $2 # The error message
    #exit 1 # return 1
    return 1 # can be called using `source` too
  fi
  #exit 0 # return 0
  return 0 # can be called using `source` too
}

function expect_file() {
    export FILE="$1"
    assert_env_nonempty $FILE "specify a filepath/name"
    if test -f "$FILE"; then
        # file exists, fine
        return 0
    else
        echo "$FILE does not exist. breaking"
        return -1
    fi
}
export -f expect_file

function MAKE_HAPPEN() {
  expect_file $1
}
export -f MAKE_HAPPEN
