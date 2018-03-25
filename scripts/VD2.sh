#!/bin/sh

if [ $# -ne 2 ]; then
  echo "USAGE: [APP NAME] [Logs Path] [PWD]"
  exit
fi

logs=$1
vdDir=$2
jobName=$vdDir.Break
# qsub -l h_vmem=6G -cwd -S /bin/sh -V -q $servers -o $logs/vd.log -e $logs/vd.error -N $jobName $APP_PATH/scripts/execMatlab.sh breakExData
$APP_PATH/scripts/execMatlab.sh `pwd -P` breakExData

echo $jobName > Hold.validate

$APP_PATH/scripts/createVD.sh ExTargets > tmp.$$


hold_jobs=`cat tmp.$$ | tail -1`
jobName=$vdDir.Collect

qsub -l h_vmem=6G -cwd -S /bin/sh -V -hold_jid $hold_jobs -q $servers -o $logs/vd.log -e $logs/vd.error -N $jobName $APP_PATH/scripts/execMatlab.sh `pwd -P` assembleVData
