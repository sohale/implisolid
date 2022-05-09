
U="****"
DM="****"
#ssh $U@$DM -- 'find /'
ssh $U@$DM -- 'ls ~/.ssh/'
#ssh $U@$DM -- 'mkdir -p ~/implisolid-on-prod; git clone hhttps://github.com/sohale/implisolid.git  ~/implisolid-on-prod/implisolid ; cd ~/implisolid-on-prod/implisolid; bash ./scripts/e2e-test-builds.bash'

