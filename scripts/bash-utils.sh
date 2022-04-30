# usage: SCRIPT_DIR=$( __current_script_dir_func )

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