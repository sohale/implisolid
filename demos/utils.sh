# usage: SCRIPT_DIR=$( script_dir_func )

function script_dir_func () {
  cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd
}
