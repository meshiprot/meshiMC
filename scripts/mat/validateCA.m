%********************************************************************************************
% This file specifies the parameters for the CoefOptimization.repeatedCoefOptimization program.
% It will be added to the CopyrightAndReadme field of the resulted experimentData object.
% ********************************************************************************************
function collect_configurations(iFold)
 tic
	pv
	load(['../experimentDataPartitions.mat']);
	load(exPath);
	exData = experimentData.extractTargets(experimentDataPartitions.testSetArray{iFold},['test_set.' num2str(iFold)]);
	exData = exData.duplicate(exData.name);
	load('ca.mat');
	[stats validationData] = CoefOptimization.validate(exData,ca,'',10);
	save('vd.mat','validationData','stats');        
 toc
 end
