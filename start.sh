#!/bin/sh

function print_help(){
  echo "USAGE: meshiMC [OPTION]"
  echo "[-ca]  Create configuration array."
  echo "[-vd]  Create validation data."
  echo "[-cv]  Cross-validation."
  echo "[-cvn]  Cross-validation - With existing partition file."
}

APP_PATH=/home/cluster/users/siditom/code/meshiMC

source $APP_PATH/data/setEnv.sh

if [ $# -lt 1 ] || [ "$1" = "-help" ]; then
  print_help
  exit
fi

if [ "$1" = "-ca" ];then
#  $APP_PATH/scripts/submitCA.sh 3 params/pv_test.m test
  $APP_PATH/scripts/submitCA.sh ${@:2}
  exit
fi

if [ "$1" = "-vd" ] && [ $# -eq 3 ];then
  $APP_PATH/scripts/submitVD.sh $2 $3 
  exit
else if [ "$1" = "-vd" ] && [ $# -ne 3 ]; then
  $APP_PATH/scripts/submitVD.sh
  exit
fi
fi

if [ "$1" = "-cv" ] && [ $# -ge 4 ];then
  $APP_PATH/scripts/submitCV.sh $2 $3 ${@:4}
  exit
fi

if [ "$1" = "-cvn" ] && [ $# -ge 4 ];then
  $APP_PATH/scripts/submitCV.sh $2 $3 ${@:4}
  exit
fi

print_help
