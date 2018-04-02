

%********************************************************************************************
% This file specifies the parameters for the CoefOptimization.repeatedCoefOptimization program.
% It will be added to the CopyrightAndReadme field of the resulted experimentData object.
% ********************************************************************************************
function cross_validation(folds,seed)
 tic
	pv
	load(exPath);
        cv_data = ExperimentData.create_XFold_testingData(experimentData,str2num(folds), str2num(seed));
%	save(['cv.' num2str(seed) '.' num2str(folds) '.mat'],'cv_data');
 toc
 end

