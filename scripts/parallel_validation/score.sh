#!/bin/csh
set here = $1
set ca = $2
set targetName = $3
set resVD = $4
set results_path = $5

if ($#argv < 4) then
	echo "score.sh"
	echo "This script is used for matlab activation for getting results on a test set(experimentData)."
	echo "Assumptions: (1) [experimentData] and [ConfigurationArray]  is in the current [pwd] directory."
	echo "             (2) The MESHI package is in path: /fastspace/users/siditom/MESHI_cloud/MY_MATLAB13_tm"
	echo 
	echo "[AppName] [pwd] [ConfigurationArray] [ExperimentData] [result validation Data name]"
	echo
	exit
endif

echo $here
echo $ca
echo $targetName
echo $resVD
echo $results_path
/fastspace/matlab/bin/matlab -r "cd $here;addpath /fastspace/users/siditom/MESHI_cloud/MY_MATLAB13_tm;load $ca;load experimentData.mat; experimentData = experimentData.extractTargets({'$targetName'},[experimentData.name '_' '$targetName']); [stats experimentData] = CoefOptimization.validate(experimentData,ca,'',10);save('$results_path/$resVD','experimentData'); exit;"
