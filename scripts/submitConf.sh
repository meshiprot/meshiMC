#!/bin/sh

echo $@

if [ $# -ne 2 ]; then
  echo "USAGE: [APP NAME] [pwd] [rnd]"
  exit
fi
currDir=$1
rnd=$2
hostname
# str="path([path ':$APP_PATH/mat' ':$APP_PATH/scripts/mat'] ); cd ${currDir}/Configurations; createConfigs(num2str($rnd)) ; end; quit;"
str="path([path ':$APP_PATH/mat'] ); cd ${currDir}/Configurations; createConfigs(num2str($rnd)) ; end; quit;"

matlab -r "path([path ':$APP_PATH/mat' ':$APP_PATH/scripts/mat'] ); cd ${currDir}/Configurations; createConfigs(num2str($rnd)); quit;"
