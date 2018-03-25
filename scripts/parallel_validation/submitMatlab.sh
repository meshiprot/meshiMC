#!/bin/sh

if [ $# -ne 1 ]; then
  echo "Error - Missing script name variable."
  exit
fi
 

matlab -r "addpath $APP_PATH; addpath $MESHI_scorePATH; $1; exit;"

