# demo1-localrun.sh

python3 -m http.server 8000 &
export server_pid=$!
echo $server_pid >server_pid-$server_pid.pid

popd
pwd

# warning: MacOS-specific code
echo "click on mp5_json_code.html @"
# open -a "Google Chrome" http://localhost:8000/
open -a "Google Chrome" http://localhost:8000/mp5_json_code.html


echo "The current server PID is:"
#ps aux|grep -ie python
ps aux|grep -ie python|grep http
sleep 1
echo "kill $server_pid"
echo "\n\n\n\n ****************"

echo "python processes to kill $(ps aux|grep -ie python|grep http|cut -c 17-25 | xargs echo)"
