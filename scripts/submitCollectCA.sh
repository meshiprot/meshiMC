#!/bin/sh
date
hostname
pwd 
echo $1

echo "cd $1; addpath $APP_PATH/scripts/mat; addpath $APP_PATH/mat; collectAllConfig('Configurations'); exit;"
matlab -r "cd $1; addpath $APP_PATH/scripts/mat; addpath $APP_PATH/mat; collectAllConfig('Configurations'); exit;"

mv Configurations/pv.m .
# rm -r Configurations

