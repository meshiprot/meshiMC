

classdef OptEngine < handle
    properties(Constant=true)
        default_normalizers = {'length' 'contacts8' 'contacts15'};  %normalizers
        default_NonZeroCoefs = 15;                                  %number of non-zero coeficients (15-25)
        default_NumOfConfigurationsInCArray = 1;                    %number of configurations
        default_maxSteps = 1;                                       % maxSteps
        default_initialTemperature = 0.0001;                        %initial temperature
        default_seed = 0
    end
    properties(SetAccess='public')
        experimentData;
        
        normalizers;
        nonZeroCoefs;
        numOfConfigurationsInCArray;
        maxSteps;
        initialTemperature;
        rnd;
        
        data_10Fold;
    end
    methods
        function obj = OptEngine(experimentData)
            obj.experimentData = experimentData;
            obj.normalizers = OptEngine.default_normalizers;
            %obj.normalizers = {'length' 'contacts8' 'contacts15'};
            obj.nonZeroCoefs = OptEngine.default_NonZeroCoefs;
            %obj.nonZeroCoefs = 15;
            obj.numOfConfigurationsInCArray = OptEngine.default_NumOfConfigurationsInCArray;
            %obj.numOfConfigurationsInCArray = 1;
            obj.maxSteps = OptEngine.default_maxSteps;
            %obj.maxSteps = 1;
            obj.initialTemperature = OptEngine.default_initialTemperature;
            %obj.initialTemperature = 0.0001;
            obj.rnd = RandStream('mcg16807','Seed',OptEngine.default_seed);
            %obj.seed = 0;
        end
      
        function validationDataAVG = flow_10Fold_preExistingPartition(obj,data_10Fold,saveDataFlag)
            
            seed = obj.rnd.randi([0 intmax('int32')]);
            configurationArrays = obj.CoefOptimization_10Fold(data_10Fold);
            [validationDataArray,validationDataAVG] = OptEngine.validateDataArrays(data_10Fold.testSetArray, configurationArrays);
            
            if (saveDataFlag)
                save(['configurationArrays.'    num2str(seed) '.mat'],'configurationArrays');
                save(['ValidationdataArray.' num2str(seed) '.mat'],'validationDataArray');
                save(['ValidationdataAVG.' num2str(seed) '.mat'],'validationDataAVG');
            end
            
            %disp(validationDataAVG);
        end
        
        function validationDataAVG = flow_10Fold_withPartition(obj,saveDataFlag)
            
            
            seed = obj.rnd.randi([0 intmax('int32')]);
            data_10Fold = OptEngine.create_10Fold_testingData(obj.experimentData,obj.rnd);
            configurationArrays = obj.CoefOptimization_10Fold(data_10Fold);
            [validationDataArray,validationDataAVG] = OptEngine.validateDataArrays(data_10Fold.testSetArray, configurationArrays);
            
            if (saveDataFlag)
                save(['data_10Fold.'    num2str(seed) '.mat'],'data_10Fold');
                save(['configurationArrays.'    num2str(seed) '.mat'],'configurationArrays');
                save(['ValidationdataArray.' num2str(seed) '.mat'],'validationDataArray');
                save(['ValidationdataAVG.' num2str(seed) '.mat'],'validationDataAVG');
            end
            
           % disp(validationDataAVG);
            
        end
    
        
        % CoefOptimization_10Fold
        % Function that calculates a validation array avrage on a set of
        % train/test data. 
        %
        % assumptions - there is a 10Fold array with seed number
        % 'exData_index' '.mat' file.
         function configurationArrays = CoefOptimization_10Fold(obj,data_10Fold)

            tic
                %disp(obj);

                %disp(data_10Fold);
                
                trainSetArray = data_10Fold.trainSetArray;
                testSetArray = data_10Fold.testSetArray;
                
                if (length(trainSetArray)~=length(testSetArray))
                    disp(['trainSetArray and testSetArray length is not equal' ]);
                    exit;
                end
                for i = 1:length(trainSetArray)
                    disp(['# 10Fold run number ' num2str(i)]);
                    CO = CoefOptimization(trainSetArray{i},'CoefOptimizationTest.m');
                    CO.numberOfCoefsRange = 10;
                    configurationArrays{i} = CO.repeatedCoefOpt_NoTrainSetPatrtition(obj.normalizers, ... %normalizers
                                                                  obj.nonZeroCoefs,... %number of non-zero coeficients (15-25)
                                                                  obj.numOfConfigurationsInCArray,... %number of configurations
                                                                  obj.maxSteps,... % maxSteps
                                                                  obj.initialTemperature,... %initial temperature
                                                                  obj.rnd.randi([0 intmax('int32')]));
                end
            toc
         end
         
    end
    methods(Static=true)
         % validateDataArrays
         %
         % Assums 
         % that the files configurationArrays and data_10Fold
         % already exists after activating the function
         % CoefOptimization_10Fold. 
         %
         % Returns
         % a validationArray that contains for each configurationArray i in
         % the configurationArrays array, validationData which co-responds
         % to the data_10Fold.testSet{i}. it also creates a
         % validationDataAVG of all the above validationData.
         function [validationDataArray validationDataAVG] = validateDataArrays(testSetArray, configArrays)
             
             validationDataArray = {};
             
             for i = 1:length(testSetArray)
                 validationDataArray{i} = Configuration.validation(testSetArray{i},configArrays{i},'');
             end
             
             validationDataAVG = OptEngine.validationDataAVG(validationDataArray);
             
         end
         
         % validationDataAVG
         % 
         % receives - a validation data array, an array that contains
         % validationData types that already been calculated.
         %
         % result - an avrage accross that array for each category.
         function validationDataAVG = validationDataAVG(validationDataArray)
             
             for i = 1:length(validationDataArray)
                 for iTarget = length(validationDataArray{i}.targetsData):-1:1
                    errorMedian(iTarget)         = validationDataArray{i}.targetsData{iTarget}.statistics.errorMedian;
                    errorWeightedMedian(iTarget) = validationDataArray{i}.targetsData{iTarget}.statistics.errorWeightedMedian;
                    errorMean(iTarget)           = validationDataArray{i}.targetsData{iTarget}.statistics.errorMean;
                    corrMedian(iTarget)          = validationDataArray{i}.targetsData{iTarget}.statistics.corrMedian;
                    corrWeightedMedian(iTarget)  = validationDataArray{i}.targetsData{iTarget}.statistics.corrWeightedMedian;
                    corrMean(iTarget)            = validationDataArray{i}.targetsData{iTarget}.statistics.corrMean;
                 end
                validationData.errorMedian(i) = median(errorMedian);
                validationData.errorWeightedMedian(i) = median(errorWeightedMedian);
                validationData.errorMean(i) = median(errorMean);
                validationData.corrMedian(i) = median(corrMedian);
                validationData.corrWeightedMedian(i) = median(corrWeightedMedian);
                validationData.corrMean(i) = median(corrMean);
             end            

            validationDataAVG.errorMedian = mean(validationData.errorMedian);
            validationDataAVG.errorWeightedMedian = mean(validationData.errorWeightedMedian);
            validationDataAVG.errorMean = mean(validationData.errorMean);
            validationDataAVG.corrMedian = mean(validationData.corrMedian);
            validationDataAVG.corrWeightedMedian = mean(validationData.corrWeightedMedian);
            validationDataAVG.corrMean = mean(validationData.corrMean);
             
         end
         
         % create_10Fold_testingData
         % 
         % args - experimentData , randomStream
         % return - a data array with 10 different partitions of the data
         % in experimentData 
         % to train/test sets using the randomStream variable. it also
         % saves the new data array to disk.
         function data_10Fold = create_10Fold_testingData(experimentData,randomStream)
            x = 10;
            trainingSetFraction = 0.9;
            seed = randomStream.Seed;
            
            for i = 1:x
                [trainSet, testSet] = experimentData.partitionExData(trainingSetFraction,randomStream);
                data_10Fold.trainSetArray{i} = trainSet;
                data_10Fold.testSetArray{i} = testSet;
            end
            
            data_10Fold.seed = seed;
            data_10Fold.seed = seed;
            experimentDataPartitions = data_10Fold;
            save(['experimentDataPartitions.mat'],'experimentDataPartitions');
         end
         
         
         function createNewGDT_TS(ed)
            targetsNum = length(ed.targetsData);

            for i = 1:targetsNum 

                B = sortrows(ed.targetsData{i}.values,9);
                [modelsNum, fieldsNum] = size(B);

                for j = 1: modelsNum

                end

            end
         end
    end
end
