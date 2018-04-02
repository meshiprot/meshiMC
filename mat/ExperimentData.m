%DISCAUMER - This MATLAB class was written by Chen Keasar, BGU (chen.keasar@gmail.com). 
%it is free without any restrictions. Yet, please be advised that it was not tested very carefully. 
%Thus, the author takes NO responsibility to any damage that this class may cause to your computer, research, or health.  
%If you find bugs or add a new functionality please drop me a line (chen.keasar@gmail.com).

classdef ExperimentData < handle
    %A container class for TargetData objects (see the TargetData.m file). 
    %A targetData object contains information regarding multiple models of the same protein (target). 
    %It is assumed that all data was collected in a uniform procedure described in the readme field (an experiment),
    %  and that all targetTata elements use the same fields in the same order. 
    %A typical way to generate a data set is to use the static method getExperimentData(experimentDirectory), 
    %  where experimentDirectory is a directory that includes multiple TargetData files.  
    properties
        name %Name of the data set
        targetsData 
        copyrightAndReadme % Free format description of the data and terms of use. May be printed using the printReadme method.
        flagComparisonTable % 0 - No comparison table data. 1 - there is comparison table data.
        statistics
    end
    methods
        
        %constructor
        function obj = ExperimentData(name)
            obj.name = name;
            obj.flagComparisonTable = 0;
        end
                
        %Duplicate this object
        function dup = duplicate(obj,newName)
            dup = ExperimentData(newName);
            dup.flagComparisonTable = obj.flagComparisonTable;
            for i = length(obj.targetsData()):-1:1
                dup.targetsData{i} = obj.targetsData{i}.duplicate(dup);
            end
        end
        
        %The data fields in this dataset
        function fields = getFields(obj)
            fields = obj.targetsData{1}.fields;
        end

        %A list of all targets
        function targetNames = getTargetNames(obj)
            targetNames = cellfun(@(X)X.targetName,obj.targetsData,'UniformOutput',false);
        end
        
        
        % Generate a new object with a copy of some fields
        % featureNames - a cell array with field names { 'field1' 'field2' ...} 
        function extracted = extract(obj,featureNames,newName)
            extracted = ExperimentData(newName);
            for i = length(obj.targetsData()):-1:1
                if (isempty(obj.targetsData{i}))
                    disp(['Target number ' num2str(i) ' is empty.']);
                else
                    extracted.targetsData{i} = obj.targetsData{i}.extract(featureNames,extracted,obj.flagComparisonTable);
                end
            end
            extracted.flagComparisonTable = obj.flagComparisonTable;
            
        end
        
        
        % Generate a new object with a copy of some targets
        % targetNames - a cell array with field names { 'field1' 'field2' ...} 
        function extractedTargets = extractTargets(obj,targetNames,newName)
            extractedTargets             = ExperimentData(newName);
            extractedTargets.flagComparisonTable = obj.flagComparisonTable;
            allTargetNames               = obj.getTargetNames;
            indices                      = Util.getIndices(allTargetNames,targetNames);
            extractedTargets.targetsData = obj.targetsData(indices);
        end
        
        %In order to break the data set into training and test subsets,
        %assuming that you already know the training set (a cell array with field names { 'field1' 'field2' ...} 
        
        function validationSet = extractValidationSet(obj,trainingSet)
           targetNames          = obj.getTargetNames;
           trainingSetIndices   = Util.getIndices(targetNames,trainingSet);
           validationSetIndices = true(length(targetNames),1);
           validationSetIndices(trainingSetIndices)= false;
           validationSet        = targetNames(validationSetIndices);
           validationSet        = obj.extractTargets(validationSet,[obj.name 'ValidationSet']);
        end
        
        %Generates a new object with a single field (newFieldName) whose
        %values are the subtraction of two fields in the this object
        %(fieldsToSubtract - a cell array).
        function difference = subtract(obj,fieldsToSubtract,newFieldName, newName)
            difference             = ExperimrentData(newName);
            for i = length(obj.targetsData()):-1:1
                difference.targetsData{i} = obj.targetsData{i}.difference(fieldsToSubtract,newFieldName,difference);
            end
        end
            
        % Generates a new object whose fields are a union of this object fields and those of the otherExperimentData argument 
        % Both this and other objects are expected to include same targets in the same order. 
        % Within each target it is expected that the same models appear in the same order.
        function concatenated = concatenate(obj, otherExperimentData,newName)
                concatenated         = ExperimentData(newName);
                if (length(obj.targetsData) ~= length(otherExperimentData.targetsData))
                    error('This is weird');
                end
                for i = length(obj.targetsData):-1:1
                    concatenated.targetsData{i} = obj.targetsData{i}.concatenate(otherExperimentData.targetsData{i},concatenated);
                end
        end
        
        
        function union_ex = union(obj, otherExperimentData,newName)
                union_ex         = ExperimentData(newName);
                
                for i = length(obj.targetsData):-1:1
                    union_ex.targetsData{i} = obj.targetsData{i}.duplicate(union_ex);
                end
                
                for i = length(otherExperimentData.targetsData):-1:1
                    union_ex.targetsData{i+length(obj.targetsData)} = otherExperimentData.targetsData{i}.duplicate(union_ex);
                end
        end
        
        function filter_LeaveOnlyNumberOfModels(obj,numberOfModels)
	    obj.filterByNumberOfModels(numberOfModels);
            cellfun(@(X)X.chooseNModels(numberOfModels),obj.targetsData); 
        end

	function [db_x db_y] = createSimilarityDB(obj,numberOfModels,numberOfRecords)
		db_x = zeros(numberOfRecords,numberOfModels*numberOfModels*3);
		db_y = zeros(numberOfRecords,numberOfModels);
		numberOfTargets=length(obj.targetsData);
		for iRecord = 1:numberOfRecords
			[line tQuality] = obj.targetsData{randperm(numberOfTargets,1)}.createSimilarityPredictionLiner(numberOfModels,'predictionWeightedMedian','gdt_ts','gdt_ts');
			db_x(iRecord,:) = line;
			db_y(iRecord,:) = tQuality;
		end
	end
	function [db_x db_y] = createFeatureDB(obj,numberOfModels)
		obj.filterByNumberOfModels(numberOfModels);
		numberOfRecords=numberOfModels*length(obj.targetsData);
		db_x = zeros(numberOfRecords,length(CoefOptimization.DefaultFeaturesToLearn));
		db_y = zeros(numberOfRecords,1);
		numberOfTargets=length(obj.targetsData);
		for iTarget = 1:numberOfTargets
			[line tQuality] = obj.targetsData{iTarget}.createFeaturesPredictionLiner(numberOfModels,'gdt_ts');
			db_x(((iTarget-1)*numberOfModels+1):(iTarget*numberOfModels),:) = line;
			db_y(((iTarget-1)*numberOfModels+1):(iTarget*numberOfModels)) = tQuality;
		end
	end
        
        function filterByNumberOfModels(obj,minimalNumberOfModels)
            obj.targetsData(cellfun(@(X)X.nModels,obj.targetsData) <minimalNumberOfModels) = [];
        end

	% For example: experimentData.filterByLeaveOnlyModelsWithName('BAKER')
        function filterByLeaveOnlyModelsWithName(obj,str)
		for i=1:length(obj.targetsData)
			obj.targetsData{i} = obj.targetsData{i}.extractModelsByNameList({str});
		end
		obj.filterByNumberOfModels(5);
        end
        
        
        function filterByGdt_ts(obj,minimalGdt_ts,maximalGdt_ts,minimalNumberOfModels)
            for i=1:length(obj.targetsData)
                obj.targetsData{i}.filterByGdt_ts(minimalGdt_ts,maximalGdt_ts);
            end
            obj.targetsData( (cellfun(@(X)X.nModels,obj.targetsData) <minimalNumberOfModels)) = [];
            if (obj.flagComparisonTable == 1)
                obj.calculateCMTWeights();
            end
        end
        
        function filterRandDuplicateModels(obj)
	  for i=1:length(obj.targetsData)
		obj.targetsData{i}.filterRandDuplicateModels();
	  end
	  if (obj.flagComparisonTable == 1)
		obj.calculateCMTWeights();
	  end
        end
        

        function modified = applySigmoidToGdtTs(obj, sigmoid, objectiveFunction, newName)
            modified = ExperimentData(newName);
            
            for i = length(obj.targetsData()):-1:1
                modified.targetsData{i} = obj.targetsData{i}.applySigmoidToGdtTs(sigmoid,objectiveFunction);
            end
            modified.flagComparisonTable = obj.flagComparisonTable;
        end

        function readReadme(obj,fileName,startLine)
            if (nargin == 2)
                startLine = 1;
            end
            readmeFid = fopen(fileName,'r');
            if (readmeFid == -1)
                error(['Cannot find file ' experimentDirectory '/readme']);
            end
            line = 'temp';
            iLine = startLine;
            while (ischar(line))
                line = fgets(readmeFid);
                obj.copyrightAndReadme{iLine} = line;
                iLine = iLine+1;
            end
        end
        
        function printReadme(obj)
            for iLine = 1:length(obj.copyrightAndReadme)-1
                fprintf(1,obj.copyrightAndReadme{iLine}(1:end));
            end
        end
        
        function [x,y] = plotFields(obj,fields,format,xFactor,yFactor,figureHandle)
            figure(figureHandle);
            extracted    = obj.extract(fields,'temp');
            x = []; y = [];
            for iTarget = length(extracted.targetsData):-1:1
                if(isempty(obj.targetsData{iTarget}))
                    disp(['Target number ' num2str(iTarget) ' is empty']);
                else
                    if (length(fields) == 2)
                        plot(extracted.targetsData{iTarget}.values(:,1)*xFactor,extracted.targetsData{iTarget}.values(:,2)*yFactor,format);
                        x = [x;extracted.targetsData{iTarget}.values(:,1)];
                        y = [y;extracted.targetsData{iTarget}.values(:,2)];
                    else
                        plot(extracted.targetsData{iTarget}.values(:,1*100)-extracted.targetsData{iTarget}.values(:,2)*100,...
                             extracted.targetsData{iTarget}.values(:,3),format);
                        x = [x;extracted.targetsData{iTarget}.values(:,1)-extracted.targetsData{iTarget}.values(:,2)];
                        y = [y;extracted.targetsData{iTarget}.values(:,3)];
                    end
                end                
                hold on
            end
        end
        
        function extracted = extractModels(obj, pattern)
            extracted = ExperimentData([obj.name '_' pattern]);
            extracted.copyrightAndReadme = obj.copyrightAndReadme;
            for iTarget = length(obj.targetsData):-1:1
                extracted.targetsData{iTarget} = obj.targetsData{iTarget}.extractModels(pattern,obj);
            end
        end    
        
        function statistic = getStatisticsField(obj,field)
            for iTarget = length(obj.targetsData):-1:1
                if(isempty(obj.targetsData{iTarget}))
                    disp(['Target number ' num2str(iTarget) ' is empty']);
                else
                    statistic(iTarget,:) = obj.targetsData{iTarget}.statistics.(field);
                end
            end
        end
    
        function getStatistics(obj,statistic)
            obj.statistics.(statistic) = getStatisticsField(obj,statistic);
        end
        
        
        function [trainSet, testSet] = partitionExData(obj,trainingSetFraction,seed)
      	    randomStream = RandStream('mcg16807','Seed',seed);    
            partition_id = randomStream.Seed;
            train_targetNames = ExperimentData.trainingSetNamesPartition(obj, trainingSetFraction, randomStream);
            trainSet = obj.extractTargets(train_targetNames,[obj.name '_training_set' num2str(partition_id)]);
            
            [test_targetNames, ia] = setdiff(obj.getTargetNames(),train_targetNames);
            testSet = obj.extractTargets(test_targetNames,[obj.name '_test_set' num2str(partition_id)]);
        end
    
        
	%% This Function calculates the statistics in each of the Targets.
        % objectiveFunction - is the function that is the model true value.
        % scoreList - a list(strings) of score functions used to predict
        % the model quality.
        %
        % stats:
        % Correlation - each line i would be the i target correlation between the
        % `objective_function` and `scoreList`. 
        % Loss - each line contains the true value (`objectiveFunction`) of
        % the best model selected by the `scoreList`
                
        function stats = calculateStatistics(obj, objective_function, scoreList,delta_flag)
           for iT = 1:length(obj.targetsData) 
		if (nargin > 3)
				obj.targetsData{iT}.calculateStatistics(objective_function,scoreList,delta_flag);
		else
   			obj.targetsData{iT}.calculateStatistics(objective_function,scoreList);
		end 
           end 
		   
           
           stats_fields = {'target_corr','MaxPredictedModelScore','iMaxPredictedModelScore','enrichment10','enrichment5','error','median'};
		   
           obj.statistics.ref_fields = [objective_function scoreList];
           for iStat = 1:length(stats_fields) 
				obj.getStatistics(stats_fields{iStat});
           end 
		   
        end 
        
        function calculateCMTWeights(obj)
            if (obj.flagComparisonTable == 0)
                return;
            end
            
            for iTarget = 1:length(obj.targetsData)
                obj.targetsData{iTarget}.calculateCMTWeights();
            end
            
        end
    
    
	    function fold = partitionExDataToFolds(obj,foldFractionSize,seed)
        	  randomStream = RandStream('mcg16807','Seed',seed);    
          	  partition_id = randomStream.Seed;
            
          	  exData = obj.duplicate([obj.name]);
          	  num_of_folds = ceil(1/foldFractionSize);
             	  num_of_targets_in_fold = floor(length(exData.targetsData)/num_of_folds);
            
        	  target_names = obj.getTargetNames();
            
       	    	  for i = (1:num_of_folds-1)
       	         
        	    fold_targetNames = ExperimentData.foldNamesPartition(target_names, num_of_targets_in_fold, randomStream);
                
  	            %fold{i} = exData.extractTargets(fold_targetNames,[exData.name '_fold_no' num2str(i) '_' num2str(partition_id)]);    
	            fold{i} = fold_targetNames;
	            %[remaining_targetNames, ia] = setdiff(exData.getTargetNames(),fold_targetNames);
	            [remaining_targetNames, ia] = setdiff(target_names,fold_targetNames);
          	    %exData = exData.extractTargets(remaining_targetNames,exData.name);
        	    target_names = remaining_targetNames;

       	    	 end
	         %fold{num_of_folds} = exData.duplicate([exData.name '_fold_no' num2str(i) '_' num2str(partition_id)]);
	         fold{num_of_folds} = remaining_targetNames;

	  end

	function exD = conjunctionExData(obj, experimentData)
            exD = ExperimentData([obj.name '_&_' experimentData.name]);
            exD.copyrightAndReadme = experimentData.copyrightAndReadme;
            for iTarget = length(obj.targetsData):-1:1
		modelsNameList = regexprep(obj.targetsData{iTarget}.fileNames,'out.*.pdb.xml|pdb','');
		jTarget =  Util.getIndices(experimentData.getTargetNames,{obj.targetsData{iTarget}.targetName});
                exD.targetsData{iTarget} = experimentData.targetsData{jTarget}.extractModelsByNameList(modelsNameList);
            end
		
	end

	% Receives another experimentData and checks if all of the decoys in obj appears in it. i.e. experimentData contains obj.
	function b = compareExData(obj,experimentData)
		exBase = obj;
		exNext = experimentData;
		b = zeros(length(exBase.targetsData),1)';
		for iTarget = 1:length(exBase.targetsData)
			try
			disp(exBase.targetsData{iTarget}.targetName);
			exTarget = exNext.extractTargets({exBase.targetsData{iTarget}.targetName},exBase.targetsData{iTarget}.targetName);
			changedFileNames = regexprep(regexprep(exTarget.targetsData{1}.fileNames,'out.*.pdb.xml','out.1.pdb.xml'),'.scwrl','');
			%setdiff(exBase.targetsData{iTarget}.fileNames, changedFileNames)	
     			b(iTarget) = length(setdiff(exBase.targetsData{iTarget}.fileNames, changedFileNames));
			if (b(iTarget) > 0)
				setdiff(exBase.targetsData{iTarget}.fileNames, changedFileNames)
			end
			catch error
				disp(error)
			end
		end
		

	end
    end
    methods(Static=true)
        function [experimentData, nModels]=getExperimentData(experimentDirectory)
            experimentData = ExperimentData(experimentDirectory);
            experimentData.readReadme([experimentDirectory '/readme']);
            d = dir([experimentDirectory '/*.targetData.mat']);
            nModels = zeros(length(d),1);
            for i = length(d):-1:1
                load([experimentDirectory '/' d(i).name]);
                nModels(i) = targetData.nModels;
%                if (i < length(d))
%                    if (targetData.differentFields(experimentData.targetsData{length(d)}))
%                        error('This is weird');
%                    end
%                end
                targetData.experimentData = experimentData;
                experimentData.targetsData{i} = targetData;
            end
        end
        
        
        function trainingSetTargetNames = trainingSetNamesPartition(experimentData,trainingSetFraction,randomStream)
            targetNames = experimentData.getTargetNames;
            randomIndices = randomStream.randperm(length(targetNames));
            trainingSet = targetNames(randomIndices(1:round(length(targetNames)*trainingSetFraction)));
            for iTarget = 1:length(trainingSet)
                if (~isempty(find(strcmp(trainingSet{iTarget},{}),1)))
                    disp([trainingSet{iTarget} ' removed from training set']);
                    trainingSet{iTarget} = [];
                end
            end
            trainingSetTargetNames = trainingSet(~cellfun('isempty',trainingSet));
        end
        
       function foldTargetNames = foldNamesPartition(targetNames,fold_targets_num,randomStream)
            %targetNames = experimentData.getTargetNames;
            disp('target names = ');
	    disp(targetNames);
            disp(['target names size = ' num2str(length(targetNames))]);
	    disp(['fold_targets_num = ' num2str(fold_targets_num)]);
		

	    randomIndices = randomStream.randperm(length(targetNames));
            trainingSet = targetNames(randomIndices(1:fold_targets_num));
            for iTarget = 1:length(trainingSet)
                if (~isempty(find(strcmp(trainingSet{iTarget},{}),1)))
                    disp([trainingSet{iTarget} ' removed from training set']);
                    trainingSet{iTarget} = [];
                end
            end
            foldTargetNames = trainingSet(~cellfun('isempty',trainingSet));
       end
        

	
         % create_XFold_testingData
         % 
         % args - experimentData , randomStream
         % return - a data array with 10 different partitions of the data
         % in experimentData 
         % to train/test sets using the randomStream variable. it also
         % saves the new data array to disk.
         function xFoldData = create_XFold_testingData(experimentData,num_of_folds,seed)
            disp(experimentData);
	    foldFraction = (1/num_of_folds);
            
            target_names = experimentData.getTargetNames();
            folds = experimentData.partitionExDataToFolds(foldFraction,seed);
            for i = 1:num_of_folds
               
                %xFoldData.trainSetArray{i} = experimentData.extractTargets(remaining_targetNames,['trainSet_' num2str(i)]);
                 xFoldData.testSetArray{i} = folds{i};
                [remaining_targetNames, ia] = setdiff(target_names, xFoldData.testSetArray{i});
                xFoldData.trainSetArray{i} = remaining_targetNames;
            end
            xFoldData.seed = seed;
            experimentDataPartitions = xFoldData;
            save(['experimentDataPartitions.mat'],'experimentDataPartitions');
            
         end
        
              
         function experimentData = validation2experiment(validationData)
             experimentData = validationData.duplicate(validationData.name);
             index = Util.getIndices(experimentData.targetsData{1}.fields,{'predictionWeightedMean'});
             for iTarget = 1:length(validationData.targetsData)
                 experimentData.targetsData{iTarget} = TargetData(experimentData.targetsData{iTarget}.targetName,experimentData.targetsData{iTarget}.fileNames,...
                     experimentData,experimentData.targetsData{iTarget}.fields(1:index-1),...
                     experimentData.targetsData{iTarget}.values(:,1:index-1));
             end
         end
   	 function experimentData = validation2argrExperiment(validationData,newScoreName,currentScoreInVD)
		experimentData = ExperimentData.validation2experiment(validationData);
	 	for iTarget = 1:length(experimentData.targetsData)
			iPred = Util.getIndices(validationData.targetsData{iTarget}.fields,currentScoreInVD);
			experimentData.targetsData{iTarget}.fields = [experimentData.targetsData{iTarget}.fields newScoreName];
			experimentData.targetsData{iTarget}.values(:,(end + 1):(end + length(iPred))) = validationData.targetsData{iTarget}.values(:,iPred);
		end
	 end
        
    end
end
    
    
    
    
