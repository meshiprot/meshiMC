#!/bin/sh

source /home/cluster/users/siditom/code/meshiMC/data/setEnv.sh

$APP_PATH/scripts/submitCA.sh 1 params/CASP13_params.m CASP13_Score
