#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Error - Missing script name variable."
  echo "Usage [APP NAME] [CURR PARH] [MATLAB SCRIPT NAME]"
  exit
fi

date
hostname
pwd
echo $1
echo $2
 
echo "cd $1; addpath $APP_PATH/scripts/mat; addpath $MESHI_scorePATH; $2; exit;"
matlab -r "cd $1; addpath $APP_PATH/scripts/mat; addpath $MESHI_scorePATH; $2; exit;"

