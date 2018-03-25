#!/bin/sh
date
hostname
pwd -P

MESHI_scorePATH=/home/cluster/users/siditom/code/meshiMC/mat
APP_PATH=/home/cluster/users/siditom/code/meshiMC/scripts/parallel_validation
export MESHI_scorePATH
export APP_PATH

if [ $# -ne 2 ]; then
  echo "Usage: [App Name] [Experiment Data Path] [MESHI_score Path]"
  echo "Clarifications: [MESHI_score Path] is in configurationArray form (i.e. ca.mat)"
  exit
fi

testSetPath=$1
MScorePath=$2
vdDir="vd_$$"

mkdir $vdDir
touch $vdDir/meta.dat
echo "ExperimentData - $testSetPath" >> $vdDir/meta.dat
echo "MESHI_score/ConfigurationArray - $MScorePath">> $vdDir/meta.dat
echo "ValidationData id - $$" >> $vdDir/meta.dat

cp $testSetPath $vdDir/experimentData.mat
cp $MScorePath $vdDir/ca.mat

cd $vdDir


logs=`pwd -P`/logs
mkdir $logs

jobName=$vdDir.Break
# qsub -cwd -S /bin/sh -V -q intel_all.q -o $logs/vd.log -e $logs/vd.error -N $jobName $APP_PATH/submitMatlab.sh breakExData
$APP_PATH/submitMatlab.sh breakExData

echo $jobName > Hold.validate

$APP_PATH/qHoldJobs.sh Hold.validate 60

$APP_PATH/create_validation_data.sh ExTargets > tmp.$$

hold_jobs=`cat tmp.$$ | tail -1`
jobName=$vdDir.Collect

qsub -cwd -S /bin/sh -V -hold_jid $hold_jobs -q keasar.q -o $logs/vd.log -e $logs/vd.error -N $jobName $APP_PATH/submitMatlab.sh assembleVData
