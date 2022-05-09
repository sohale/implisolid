
U="ssthere"
DM="157.245.40.6"
#ssh $U@$DM -- 'find /'
ssh $U@$DM -- 'ls ~/.ssh/'
# make it one-liner
SCRIPT_TO_RUN=<<< MYSCRIPT
mkdir -p ~/implisolid-on-prod; git clone https://github.com/sohale/implisolid.git  ~/implisolid-on-prod/implisolid ; cd ~/implisolid-on-prod/implisolid; bash ./scripts/e2e-test-builds.bash
MYSCRIPT
#ssh $U@$DM -- "$SCRIPT_TO_RUN"
