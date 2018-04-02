

%********************************************************************************************
% This file specifies the parameters for the CoefOptimization.repeatedCoefOptimization program.
% It will be added to the CopyrightAndReadme field of the resulted experimentData object.
% ********************************************************************************************
function createConfigs(seed,nFold)
 tic
	pv
	load(exPath);
        if ( nargin > 1 )
		load(['experimentDataPartitions.' num2str(partition_seed) '.mat']);
		exData = expeimentData.extractTargets(experimentDataPartitions.trainSetArray{fold_num},['train_set.' num2str(fold_num)]);
		exData = exData.duplicate(exData.name);
	end
        CO = CoefOptimization(experimentData,'myCoefOptimization.m');
	CO.numberOfCoefsRange = numberOfCoefsRange;
	CO.configEnergy = configEnergy;
	CO.objectiveFunction = objectiveFunction;
	CO.deltaObjective = deltaObjective;
	CO.learningFraction = learningFraction;
	CO.featuresToLearn = featuresToLearn;
 
        configurationArray =  CO.repeatedCoefOpt_NoTrainSetPatrtition(normalizers,nonZeroCoefsNum,numOfConfigsInArray,MCMCSteps,initialTemp,str2num(seed));

	configuration = configurationArray;
	save(['configuration.' num2str(seed) '.mat'],'configuration');
 toc
 end

