#!/bin/sh 

date
hostname
pwd

data_path=$1
APP_PATH=/fastspace/users/siditom/MESHI_cloud/scripts/Matlab/parallel_validation
if [ $# -lt 1 ]; then
  echo "create_validation_data.sh"
  echo
  echo "This script is for parralel creation of validation data."
  echo "Assumptions: (1) the Experiment Data File contains the list of targets in the experimentData."
  echo "             (2) There is a configuration array ca.mat in the current path."
  echo
  echo "[AppName] [ExperimentData Folder] "
  exit
fi
hold_jobs=","
if [ $# -eq 2 ];then
  hold_jobs=$2
  echo "UPDATE: hold_jobs=$2"
fi 

dirs=`cat $data_path`

run_folder=`echo $data_path | awk -F "/" '{print $NF}'`_VD
mkdir $run_folder

if [ ! -e logs ]; then
  mkdir logs
fi



echo $dirs

current_jobs=","
#sumiting jobs to calculat validation data per target
for target_name in $dirs
do
    echo $target_name processing...
      #echo "score.sh `pwd` ca.mat $data_path/$target_name VD_$target_name $run_folder" 
    job_name=${run_folder}_${target_name}

      qsub -N $job_name -hold_jid $hold_jobs -q intel_all.q -e `pwd -P`/logs/$target_name.error -o `pwd -P`/logs/$target_name.log $APP_PATH/score.sh `pwd -P` ca.mat $target_name VD_$target_name $run_folder 
    
    current_jobs="${current_jobs},${job_name}"
done

echo ${current_jobs}
