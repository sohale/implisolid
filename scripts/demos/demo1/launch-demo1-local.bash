#!/bin/bash

set -ex
set -u
function assert_env_nonempty() {
  if [ "$#" -ne 2 ]; then
    echo "assert_env_nonempty: expected 2 arguments, got $#: $*"
    return 1
  fi

  if [ ".$1" = "." ]; then
    echo "shell env is empty"; echo $2
    return 1
  fi
}

export MYPORT8000=8099

# Runs the deployed one
# Currennt folder (pwd) should be where the served files are (ie the app) (root of url resources)
# alt name: demo1-localrun.sh

#args:
assert_env_nonempty "$APP_RUN_LOCATION" "env-argument APP_RUN_LOCATION= missing"
# APP_RUN_LOCATION was DEPLOY_LOCATION
# Template

# Kill the previous server process
ps aux|grep python|grep http.server |cut -c10-17 | xargs kill || :

cd $APP_RUN_LOCATION
echo "Running python server from: $(pwd)"
# Keep `--bind 127.0.0.1` to access locally via 127.0.0.1:8000 . Login should use `-L 8000:127.0.0.1:8000`. It would be "inbound SSH to your Mac". Never bind to 0.0.0.0 on an internet-facing host.
python3 -m http.server $MYPORT8000 --bind 127.0.0.1 &
export server_pid=$!
echo $server_pid >$APP_RUN_LOCATION/server_pid-$server_pid.pid


# cd $APP_RUN_LOCATION/js

public_ip="$(curl https://ipinfo.io/ip)"
echo "public ip: $public_ip"
echo "http://${public_ip}:$MYPORT8000/mp5_json_code.html"
GREEN="\e[1;32m" RESET="\e[0m"
echo -e "Click here: ${GREEN}http://${public_ip}:$MYPORT8000/mp5_json_code.html${RESET}"

echo "OSTYPE= $OSTYPE"

echo "click on mp5_json_code.html @"
[[ $OSTYPE == 'darwin'* ]] || "Warning: MacOS-specific code: for `open`"
[[ $OSTYPE == 'darwin'* ]] || \
open -a "Google Chrome" http://localhost:$MYPORT8000/mp5_json_code.html


echo "The current server PID is:"
ps aux|grep -ie python|grep http
sleep 1
echo "kill $server_pid" | tee -a $APP_RUN_LOCATION/js/processes-to_kill.log
printf "\n\n\n\n ****************"

echo "python processes to kill $(ps aux|grep -ie python|grep http|cut -c 17-25 | xargs echo)"
cat $APP_RUN_LOCATION/js/processes-to_kill.log || :


<< ////
  Three ways to log the pid  of the process to kill
     ./processes-to_kill.log,
     server_pid-*.pid
     ps|cut, and
     direct console out "kill ..."
////

pwd
