#!/bin/sh

if [ $# -lt 3 ]; then
  echo "USAGE: [APP NAME] [No. of Folds] [Parameters File - ExperimentData location] [Parameters file path 1] [Parameters file path 2] ... [Parameters file path i]"
  exit
fi

echo "Not Supported Yet."

echo "numOfARGS=$#"
nFolds=$1
pv=`realpath $2`

# create test folder
mkdir "cv_$$"
cvDir=`realpath cv_$$`
mkdir $cvDir/logs
logs=`realpath $cvDir/logs`
cp $pv $cvDir/pv.m

cd $cvDir

# copy params info
mkdir params
i=1
for params in "${@:3}"
do
 cp $params params/${i}.m
 i=$(($i + 1))
done

# create cv_data
  job_name=cv_$$
  echo "qsub -l h_vmem=6G -S /bin/sh -cwd -V -N $job_name -q $servers -o $logs/$job_name.log -e $logs/$job_name.error ${APP_PATH}/scripts/execMatlab.sh $cvDir cross_validation(num2str($nFolds),num2str($RANDOM))"
  qsub -l h_vmem=6G -S /bin/sh -cwd -V -N $job_name -q $servers -o $logs/$job_name.log -e $logs/$job_name.error ${APP_PATH}/scripts/execMatlab.sh $cvDir "cross_validation(num2str($nFolds),num2str($RANDOM))" 
  echo "$job_name" >> dependencies.$$.tmp

# submit train/validate for each parameters file 
hjids=""
for iFold in $(seq 1 $nFolds)
do
  iParam=0
  for params in "${@:3}"
  do
    iParam=$(($iParam + 1))
    caJobs=`cat $params | grep cvJobs | awk '{print $NF}' | tr -d ';'`

    ${APP_PATH}/scripts/submitCA.sh $caJobs $params fold${iFold}.${iParam} ${iFold} `realpath dependencies.$$.tmp` 

    depedencies=`cat fold${iFold}.${iParam}/CAJobs.*.tmp | tr '\n' ','`
    jobName="vd_fold${iFold}.${iParam}"
    hjids="${hjids},${jobName}"
    qsub -l h_vmem=10G -S /bin/sh -V -cwd -N "$jobName" -hold_jid $depedencies -q $servers -V -o $logs/$jobName.log -e $logs/$jobName.error $APP_PATH/scripts/execMatlab.sh ${cvDir}/fold${iFold}.${iParam} "validateCA(${iFold})"
  done
done

# collect validation

jobName="collect_cv.${nFolds}.${iParam}"
qsub -l h_vmem=10G -S /bin/sh -V -cwd -N "$jobName" -hold_jid $hjids -q $servers -V -o $logs/$jobName.log -e $logs/$jobName.error $APP_PATH/scripts/execMatlab.sh ${cvDir} "collectCV(${nFolds},${iParam})"

