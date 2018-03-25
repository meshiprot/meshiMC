#!/bin/sh
date
hostname
pwd -P

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

qsub -l h_vmem=6G -cwd -S /bin/sh -V -N submitVD -o $logs/submitVD.log -e $logs/submitVD.error $APP_PATH/scripts/VD2.sh $logs $vdDir
