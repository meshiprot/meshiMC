#!/bin/sh

if [ $# -ne 3 ]; then
	echo "USAGE: [APP NAME] [MULTIPLICATIONS] [PATH TO PARAMETERS] [OUTPUT PATH]"
	echo "[MULTIPLICATIONS] - number of times to create a configuration array with parameters [PATH TO PARAMETERS]"
	exit
fi

mul=$1
paramsPath=$2
outDir=$3

if [ -e $outDir ];then
  echo "Output directory already exists."
#  exit
fi

if [ ! -e $paramsPath ];then
  echo "Parameters file doesn't exists."
  exit
fi

if [ -z "${APP_PATH+x}" ];then
  echo "Enviroment is not set."
  exit
fi

mkdir $outDir

mkdir $outDir/Configurations
mkdir $outDir/logs
cp $paramsPath $outDir/Configurations/pv.m
rnd=$$
cd $outDir
currDir=`pwd -P`

logs=${currDir}/logs

params=`basename $paramsPath | sed 's/.m//g'`

for i in $(seq 1 $mul)
do
  rnd=$RANDOM
  job_name=Con.${params}.$rnd
  echo ${job_name} >> tmp.$$
  ${APP_PATH}/scripts/qWait.sh $MJobs $tIntervalConfig 
  qsub -l h_vmem=6G -S /bin/sh -cwd -V -N $job_name -q $servers -o $logs/$job_name.log -e $logs/$job_name.error ${APP_PATH}/scripts/execMatlab.sh $currDir/Configurations "createConfigs(num2str($rnd))" 
done 

${APP_PATH}/scripts/qWait.sh $MJobs $tIntervalConfig 

echo "qsub -l h_vmem=10G -S /bin/sh -V -cwd -N ${params}.collect -hold_jid `cat tmp.$$ | tr '\n' ','` -q keasar.q -V -o $logs/collect.log -e $logs/collect.error $APP_PATH/scripts/submitCollectCA.sh ${currDir}"
qsub -l h_vmem=10G -S /bin/sh -V -cwd -N ${params}.collect -hold_jid `cat tmp.$$ | tr '\n' ','` -q $servers -V -o $logs/collect.log -e $logs/collect.error $APP_PATH/scripts/execMatlab.sh ${currDir} collectAllConfig

rm -f tmp.$$ 
