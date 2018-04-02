#!/bin/sh
date
hostname
pwd -P

if [ $# -lt 2 ]; then
  echo "Usage: [App Name] [Experiment Data Path] [MESHI_score Path] [Job Dependecies file]"
  echo "Clarifications: [MESHI_score Path] is in configurationArray form (i.e. ca.mat)"
  exit
fi

JobDependencies="NO_JOB"
if [ ! -z "${3+x}" ];then
  JobDependencies=`cat $5 | tr '\n' ',' | tr -d ' '`
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

${APP_PATH}/scripts/qWait.sh $MJobs $tIntervalConfig 
qsub -l h_vmem=6G -cwd -S /bin/sh -V -hold_jid $JobDependencies -N submitVD -o $logs/submitVD.log -e $logs/submitVD.error $APP_PATH/scripts/parallelVD.sh $logs $vdDir
