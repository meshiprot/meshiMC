

%********************************************************************************************
% This file specifies the parameters for the CoefOptimization.repeatedCoefOptimization program.
% It will be added to the CopyrightAndReadme field of the resulted experimentData object.
% ********************************************************************************************
function myCoefOptimization(experimentData)
 tic

    CO = CoefOptimization(experimentData,'myCoefOptimization.m');
	CO.numberOfCoefsRange = 10;
	CO.configEnergy = ConfigEnergyErAndMCorrCMTWeightsGDT_TS;
    configurationArray =  CO.repeatedCoefOpt_NoTrainSetPatrtition({'length' 'contacts8' 'contacts15'}, ... %normalizers
                                                          15,... %number of non-zero coeficients (15-25)
                                                          1,... %number of configurations
                                                          10,... % maxSteps
                                                          0.00001,... %initial temperature
                                                          1);


	configuration = configurationArray{1};
	save(['configuration.mat'],'configuration');
	

 toc
 end

