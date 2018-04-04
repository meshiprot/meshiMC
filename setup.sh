#! /bin/bash
currentPath=`pwd -P`
# sed -i 's@ .*\/parameters@ '"${currentPath}"'\/parameters@g' commands
# sed -i 's@ .*\/score_functions@ '"${currentPath}"'\/score_functions@g' commands
sed -i 's@APP_PATH=.*@APP_PATH='"${currentPath}"'@g' ./data/setEnv.sh
sed -i 's@APP_PATH=.*@APP_PATH='"${currentPath}"'@g' ./start.sh
