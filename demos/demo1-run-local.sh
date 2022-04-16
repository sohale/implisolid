#!/bin/bash
# alt name: demo1-localrun.sh

source utils.sh

echo "SCRIPT_DIR: $( script_dir_func )"

DEPLOY_LOCATION=$1
#cd $DEPLOY_LOCATION
# that is, $REPO/demos

echo "Running python server from: $(pwd)"
python3 -m http.server 8000 &
export server_pid=$!
echo $server_pid >server_pid-$server_pid.pid


#popd
#pwd
cd $DEPLOY_LOCATION/js

[[ $OSTYPE == 'darwin'* ]] || "Warning: MacOS-specific code: for `open`"

echo "click on mp5_json_code.html @"
# open -a "Google Chrome" http://localhost:8000/
open -a "Google Chrome" http://localhost:8000/mp5_json_code.html


echo "The current server PID is:"
#ps aux|grep -ie python
ps aux|grep -ie python|grep http
sleep 1
echo "kill $server_pid" | tee -a processes-to_kill.log
printf "\n\n\n\n ****************"

echo "python processes to kill $(ps aux|grep -ie python|grep http|cut -c 17-25 | xargs echo)"
cat processes-to_kill.log || :


<< ////
  Three ways to log the pid  of the process to kill
     ./processes-to_kill.log,
     server_pid-*.pid
     ps|cut, and
     direct console out "kill ..."
////

pwd
