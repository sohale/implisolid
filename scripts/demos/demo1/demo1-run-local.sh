#!/bin/bash

set -ex
function assert_env_nonempty() {
  if [ ".$1" = "." ]; then
    echo "shell env is empty"; echo $2
    return 1
  fi
}

# Runs the deployed one
# Currennt folder (pwd) should be where the served files are (ie the app) (root of url resources)
# alt name: demo1-localrun.sh

#args:
assert_env_nonempty $APP_RUN_LOCATION "env-argument APP_RUN_LOCATION= missing"
# APP_RUN_LOCATION was DEPLOY_LOCATION
# Template

cd $APP_RUN_LOCATION
echo "Running python server from: $(pwd)"
python3 -m http.server 8000 &
export server_pid=$!
echo $server_pid >$APP_RUN_LOCATION/server_pid-$server_pid.pid


# cd $APP_RUN_LOCATION/js

public_ip="$(curl https://ipinfo.io/ip)"
echo "http://${public_ip}:8000/mp5_json_code.html"
echo "public ip: $public_ip"

echo "click on mp5_json_code.html @"
[[ $OSTYPE == 'darwin'* ]] || "Warning: MacOS-specific code: for `open`"
open -a "Google Chrome" http://localhost:8000/mp5_json_code.html


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
