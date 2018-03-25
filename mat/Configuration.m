classdef Configuration <handle
    properties
        fields             %
        allFieldIndices
        nFields
        normalizers
        coefs
        exponentIndices
        normalizerIndices
        scoreFunctions
        changedFields
        learningSet
        configEnergy
        allValues
        allObjectives
        numberOfCoefs
        fieldIndices
        energy
        errors
        exponentValues
        zeroIndex
        nExponentValues
        trajectory
        correlations
        indicesToReduce %features that tend to be over represented (good features but reduce diversity)
        trainingSet;
        objectiveFunction
	deltaObjectiveFunction
        sigmoid
    end
    
    methods
        function obj = Configuration(fields, normalizers,coefs, exponentIndices, exponentValues, sigmoid,...
                                     nTargets, learningFraction,numberOfCoefs,configEnergy,...
                                     coefOptimization,trainingSet,objectiveFunction,randomStream,weightsExp)
            obj.fields            = fields;
            obj.nFields           = length(fields);
            obj.normalizers       = normalizers;
            obj.normalizerIndices = Util.getIndices(fields,normalizers);
            obj.coefs             = coefs;
            obj.exponentIndices   = exponentIndices;
            obj.changedFields     = 1:length(fields);
            obj.exponentValues    = exponentValues;
            obj.indicesToReduce   = Util.getIndices(fields,CoefOptimization.featuresToReduce);
            obj.sigmoid = sigmoid;
            if (nargin > 6)
                obj.trainingSet       = trainingSet;
                obj.exponentValues    = coefOptimization.exponentValues;
                obj.nExponentValues   = length(coefOptimization.exponentValues);
                obj.exponentIndices   = randomStream.randi(obj.nExponentValues,...
                    length(fields),length(normalizers));
                numberOf7             = length(fields) - numberOfCoefs;
                randomPermutation     = randomStream.randperm(length(fields));
                obj.exponentIndices(randomPermutation(1:numberOf7),1) = obj.nExponentValues+1; %Thease features will have a zero weight by the first normalizer, and thus by all of them.
                for iTarget = nTargets:-1:1
                    obj.scoreFunctions{iTarget} = ScoreFunction(obj,coefOptimization.features.targetsData{iTarget},coefOptimization.exponentValues);
                    obj.scoreFunctions{iTarget}.getExponentWeightsMatrix();
                end
                randomArray           = randomStream.randperm(nTargets);
                obj.learningSet       = randomArray < nTargets*learningFraction;
                obj.configEnergy      = configEnergy;
                obj.numberOfCoefs     = numberOfCoefs;
                obj.allValues         = zeros(100,obj.numberOfCoefs);
                obj.allObjectives     = zeros(100,1);
                obj.zeroIndex         = find(coefOptimization.exponentValues == 0,1);
                obj.objectiveFunction = objectiveFunction;
		obj.deltaObjectiveFunction = coefOptimization.deltaObjective;
            end
        end
        
        function newConfiguration = duplicate(obj)
            newConfiguration = Configuration(obj.fields, obj.normalizers, obj.coefs,obj.exponentIndices, obj.exponentValues,obj.sigmoid);
            newConfiguration.energy            = obj.energy;
            newConfiguration.errors            = obj.errors;
            newConfiguration.correlations      = obj.correlations;
            newConfiguration.fieldIndices      = obj.fieldIndices;
            newConfiguration.trajectory        = obj.trajectory;
            newConfiguration.learningSet       = obj.learningSet;
            newConfiguration.trainingSet       = obj.trainingSet;
            newConfiguration.objectiveFunction = obj.objectiveFunction;
            newConfiguration.deltaObjectiveFunction    = obj.deltaObjectiveFunction;
        end
        function [energy, error, correlation] = evaluate(obj,coefOptimization,randomStream)
            obj.getWeights();
            obj.getValuesAndObjectives(coefOptimization,randomStream);
            tempCoefs               = obj.allValues\obj.allObjectives;
            obj.coefs               = zeros(obj.nFields,1);
            obj.coefs(obj.fieldIndices) = tempCoefs;
            [obj.energy,  obj.errors, obj.correlations] = obj.configEnergy.evaluate(obj,coefOptimization.objective);%,  coefOptimization);
            correlation = median(obj.correlations);
            error       = median(obj.errors);
            energy      = obj.energy;
        end
        
        function getWeights(obj)
            if (~isempty(obj.changedFields))
                for iTarget =1:length(obj.scoreFunctions)
			obj.scoreFunctions{iTarget}.getWeightsTarget();
                end
            end
        end
        
        function getValuesAndObjectives(obj,coefOptimization,randomStream)
            obj.fieldIndices = obj.exponentIndices(:,1) <= obj.nExponentValues;
            totalNmodels = 0;
            for iTarget = find(obj.learningSet)
                nModels = obj.scoreFunctions{iTarget}.nModels;
                modelIndices = randomStream.randperm(nModels);
                if (nModels > coefOptimization.maxNumberOfModels)
                    nModels = coefOptimization.maxNumberOfModels;
                end
                modelIndices = modelIndices(1:nModels);

                try %this line throws a very rare exception. This try/catch is here to figure out what has happend
                    obj.allValues(totalNmodels+1:totalNmodels+nModels,...
                        1:obj.numberOfCoefs) = ...
                        obj.scoreFunctions{iTarget}.weightedValues(modelIndices,obj.fieldIndices);
                catch ex
                    disp(ex.message);
                    disp(['iTarget' num2str(iTarget)]);
                    disp(['totalNmodels+1' num2str(totalNmodels+1)]);
                    disp(['nModels' num2str(nModels)]);
                    disp(['obj.numberOfCoefs' num2str(obj.numberOfCoefs)]);
                    disp('size(obj.scoreFunctions{iTarget}.weightedValues(modelIndices,obj.fieldIndices))');
                    size(obj.scoreFunctions{iTarget}.weightedValues(modelIndices,obj.fieldIndices));
                    disp('obj.fieldIndices')
                    obj.fieldIndices
                    throw(ex);
                end
                obj.allObjectives(totalNmodels+1:totalNmodels+nModels) = ...
                    coefOptimization.objective.targetsData{iTarget}.values(modelIndices,1);
                totalNmodels = totalNmodels+nModels;
            end
            obj.allValues(totalNmodels+1:end,:)=[];
            obj.allObjectives(totalNmodels+1:end)=[];
        end
        function mutate(obj,newConfiguration,randomStream,type)
            newConfiguration.exponentIndices = obj.exponentIndices;
            newConfiguration.changedFields   = obj.changedFields;
            newConfiguration.learningSet     = obj.learningSet;
            newConfiguration.trajectory      = obj.trajectory;
            for i = 1:length(obj.scoreFunctions)
                newConfiguration.scoreFunctions{i}.configuration = newConfiguration;
            end
            if (type == 0)
                type = randomStream.randi(8);
            end
            switch type
                case 1
                    newConfiguration.flipFamilies(obj,randomStream);
                case 2
                    newConfiguration.flipFeatures(obj,randomStream);
                case 3
                    newConfiguration.flipFeaturesToZero(obj,randomStream);
                case 4
                    newConfiguration.flipFeaturesToZero(obj,randomStream);
                case 5
                    newConfiguration.flipFeatures(obj,randomStream);
                case 6
                    newConfiguration.changeWeight(obj,randomStream);
                case 7
                    newConfiguration.changeWeightToZero(obj,randomStream);
                case 8
                    newConfiguration.changeWeightSmall(obj,randomStream);
                otherwise
                    error(['weird type ' num2str(type)]);
            end
        end
        
        function flipFamilies(obj,configuration,randomStream)
            ls = find(configuration.learningSet); %indices of learningSet
            es = find(~configuration.learningSet);
            
            li = randomStream.randi(length(ls)); %random pointer to the indices
            ei = randomStream.randi(length(es));
            
            lr = ls(li); %random pointer into the learning set
            er = es(ei);
            
            obj.changedFields = [];
            obj.learningSet(er) = true; % a new member of the learning set
            obj.learningSet(lr) = false; % no longer in the learning set
            obj.changedFields = 1:length(obj.fields);
        end
        
        function flipFeatures(obj,configuration,randomStream)
            [indices1,indices2,r1,r2] = obj.preFlipFeatures(configuration,randomStream);
            obj.exponentIndices                 = configuration.exponentIndices;
            obj.exponentIndices(indices1(r1),:) = configuration.exponentIndices(indices2(r2),:);
            obj.exponentIndices(indices2(r2),:) = configuration.exponentIndices(indices1(r1),:);
        end
        function flipFeaturesToZero(obj,configuration,randomStream)
            [indices1,indices2,r1,r2] = obj.preFlipFeatures(configuration,randomStream);
            obj.exponentIndices(indices1(r1),:) = ones(1,length(obj.exponentIndices(indices1(r1),:)))*...
                                                  obj.zeroIndex;
            obj.exponentIndices(indices2(r2),:) = obj.nExponentValues+1;
        end
        
        function [indices1,indices2,r1,r2] = preFlipFeatures(obj,configuration,randomStream)
            indices1 = find(configuration.exponentIndices(:,1) == obj.nExponentValues+1);
            indices2 = find(configuration.exponentIndices(:,1) <= obj.nExponentValues);
            n1 = length(indices1);
            n2 = length(indices2);
            if ((n1 == 0)|| (n2 == 0))
                disp(['obj.nExponentValues+1 ' num2str(obj.nExponentValues+1)]); 
                disp('configuration.exponentIndices');
                disp(configuration.exponentIndices')
                disp('indices1');
                disp(indices1);
                disp('indices2');
                disp(indices2);                
                error (['weird n1 or n2 ' num2str(n1) ' ' num2str(n2)]);
            end
            r1 = randomStream.randi(n1); %A feature that was not used until now and will be tested
            if (~isempty(obj.indicesToReduce))
                while ((randomStream.rand() <0.5) && (~isempty(find(obj.indicesToReduce==indices1(r1),1))))
                    r1 = randomStream.randi(n1);
                end
            end
            r2 = randomStream.randi(n2);
            
            obj.changedFields = [indices1(r1) indices2(r2)];
        end
        
        function changeWeight(obj,configuration,randomStream)
            if (size(obj.exponentIndices,2) ~= 1)
                [indices, r1, r2, e, l] = obj.preChangeWeight(configuration,randomStream);
                eNew    = e;
                %looking fo a random new value
                while (e == eNew)
                    eNew = randomStream.randi(l);
                end
                obj.exponentIndices(indices(r1),r2) = eNew;
            end
        end
        function changeWeightToZero(obj,configuration,randomStream)
            if (size(obj.exponentIndices,2) ~= 1)
                [indices, r1, r2, e, ~] = obj.preChangeWeight(configuration,randomStream);
                if (e == obj.zeroIndex)
                    obj.exponentIndices(indices(r1),r2) = e+2*randomStream.randi(2)-3;
                else
                    obj.exponentIndices(indices(r1),r2) = obj.zeroIndex;
                end
            end
        end
        
        function changeWeightSmall(obj,configuration,randomStream)
            if (size(obj.exponentIndices,2) ~= 1)
                [indices, r1, r2, e, l] = obj.preChangeWeight(configuration,randomStream);
                if (e == 1)
                    obj.exponentIndices(indices(r1),r2) = 1+randomStream.randi(2);
                else
                    if (e == l)
                        obj.exponentIndices(indices(r1),r2) = l-randomStream.randi(2);
                    else
                        obj.exponentIndices(indices(r1),r2) = e+2*randomStream.randi(2)-3;
                    end
                end
            end
        end
        
        function [indices, r1, r2, e, l] = preChangeWeight(obj,configuration,randomStream)
            obj.exponentIndices   = configuration.exponentIndices;
            indices = find(configuration.exponentIndices(:,1) <= obj.nExponentValues); %indices of columns with weight > 0
            n       = length(indices); % number of such columns
            r1       = randomStream.randi(n); %a random index into the columns list
            r2       = randomStream.randi(length(configuration.exponentIndices(1,:)));
            e       = configuration.exponentIndices(indices(r1),r2); % exponent index at that column
            l       = obj.nExponentValues;
            %looking fo a random new value
            obj.changedFields = indices(r1);
        end
        
        function [energy, weighted1, coefs, configuration] = validateEnergy(obj,targetData)
            obj.changedFields = 1:length(obj.fields);
            scoreFunction = ScoreFunction(obj,targetData,obj.exponentValues);
            predictions  = scoreFunction.calcScoreTarget();
            observations = targetData.extract({'gdt_ts'},'observations');
            diff = predictions-observations.values;
            %[observations.values predictions diff]
            energy = (sqrt(mean(diff.*diff)));
            
            weighted1 = scoreFunction.weightedValues(1,scoreFunction.configuration.fieldIndices);
            coefs = scoreFunction.configuration.coefs(scoreFunction.configuration.fieldIndices,1);
            configuration = obj;
        end
        
%         function addsScores(obj,experimentData,fieldName)
%             scoreFunction = ScoreFunction(obj,obj.exponentValues);
%             for i = 1:length(experimentData.targetsData)
%                 experimentData.targetsData{i}.fields = {experimentData.targetsData{i}.fields{:} fieldName};
%                 experimentData.targetsData{i}.values = [experimentData.targetsData{i}.values ...
%                                                         scoreFunction.calcScoreTarget(experimentData.targetsData{i})];
%             end
%         end
        
        function trainingSet = extractTrainingSet(obj,experimentData)
            trainingSet = experimentData.extractTargets(obj.trainingSet,[experimentData.name 'TrainingSet']);
        end
        
        function validationSet = extractValidationSet(obj,experimentData)
            validationSet = experimentData.extractValidationSet(obj.trainingSet);
        end
   
        
    end
    
    methods(Static=true)
        
        function energies = getEnergies(configurationArray)
            energies = cellfun(@(X)X.energy,configurationArray);
        end
        function correlations = getCorrelations(configurationArray)
            correlations = cellfun(@(X)median(X.correlations),configurationArray);
        end
        
        function errors = getErrors(configurationArray)
            errors = cellfun(@(X)median(X.errors),configurationArray);
        end
                
        function fields = getfields(configurationArray)
            fields = cellfun(@(X)X.fields(X.fieldIndices),configurationArray,'UniformOutput',false);
        end
        
        function selectedConfigurations = selectConfigurationsByEnergy(configurationArray,fraction)
            energies = Util.getEnergies(configurationArray);
            [~, sortedIndices] = sort(energies);
            selectedConfigurations = configurationArray(sortedIndices(1:round(length(sortedIndices)*fraction)));
        end
        
        function selectedConfigurations = selectConfigurationsByCorrelation(configurationArray,fraction)
            corrs = cellfun(@(X)median(X.correlations),configurationArray);
            [~, sortedIndices] = sort(corrs);
            selectedConfigurations = configurationArray(sortedIndices(round(length(sortedIndices)*(1-fraction)):end));
        end
        
        function energies = validateConfigurationArray(experimentData, configurationArray,debugFlag)
            energies = zeros(length(configurationArray),1);
            if (~isempty(experimentData.targetsData{1}.statistics))
                error('Statistics not empty')
            end
            for i = 1:length(experimentData.targetsData)
                experimentData.targetsData{i}.statistics = ones(1,length(configurationArray))*-9999.9999;
            end
            for i = 1:length(configurationArray)
                fprintf(1,'.');
                if (mod(i,50) == 0)
                        fprintf('\n');
                end
                energies(i) = Configuration.validateConfiguration(configurationArray{i},experimentData,i,debugFlag);
            end
            for i = 1:length(experimentData.targetsData)
                experimentData.targetsData{i}.statistics(experimentData.targetsData{i}.statistics == -9999.9999) = [];
            end
        end
 
        function [energy, errors, correlations] = validateConfiguration(configuration,experimentData,iConfig,debugFlag)
            if (debugFlag)
                newExperimentData = configuration.extractTrainingSet(experimentData);
            else
                newExperimentData = configuration.extractValidationSet(experimentData);
            end
            features       = newExperimentData.extract(configuration.fields,[experimentData.name 'Features']);
            objective      = newExperimentData.extract({configuration.objectiveFunction},[experimentData.name 'Objective']);
            configuration.changedFields = 1:length(configuration.fields);
            configCoefs = configuration.coefs(configuration.fieldIndices,1);
            allTargetNames = experimentData.getTargetNames;
            indices        = Util.getIndices(allTargetNames,configuration.trainingSet);
            for iTarget=length(features.targetsData):-1:1
                scoreFunction  = ScoreFunction(configuration,features.targetsData{iTarget},configuration.exponentValues);
                scoreFunction.calcScoreTarget();
                [errors(iTarget), correlations(iTarget)] = ConfigEnergy.evaluateTarget(configuration,...
                                                                                      scoreFunction,...
                                                                                      configCoefs,...
                                                                                      objective.targetsData{iTarget});
                experimentData.targetsData{indices(iTarget)}.statistics(iConfig) = errors(iTarget);
            end    
            energy = 0.99*sqrt(mean(errors.*errors))+0.01*(1-median(correlations));
        end
        
        
        function trainingConfigurations = extractTrainingConfiguration(targetData, configurationArray) 
            trainingConfigurations = {};
            for iConfig = length(configurationArray):-1:1
                if (isempty(find(strcmp(targetData.targetName,configurationArray{iConfig}.trainingSet),1)))
                    trainingConfigurations{iConfig} = configurationArray{iConfig};
                end
            end
            trainingConfigurations = trainingConfigurations(~cellfun('isempty',trainingConfigurations));
        end
        
        
        function targetDataValidation = predictionTarget(experimentData,  iTarget, configurationArray,weightsExp)
            targetData           = experimentData.targetsData{iTarget};
            targetDataValidation = TargetData(targetData.targetName, targetData.fileNames,...
                                              [], ... % experimentData
                                              {targetData.fields{:} 'predictionWeightedMean' 'predictionMedian' 'predictionWeightedMedian' ...
									'weightedMdianIDR' 'weightedMdianENT' 'predictionExpWeightedMedian' ...
									'predictionMaxWeightedMedian' 'predictionClusterWeightedMedian' 'predictionWeightedMedianCluster' ...
									'predictionRandomWeightedMedian'},...
                                              targetData.values,...
                                              targetData.comparisonTables,...
                                              targetData.fieldsComparisonTables);
            trainingConfigurations = Configuration.extractTrainingConfiguration(targetData,configurationArray);
            disp(['# of training configurations ' int2str(length(trainingConfigurations))]);
                if (isempty(trainingConfigurations))
                    targetDataValidation = [];
                else
                    predictionIndices = Util.getIndices(targetDataValidation.fields,{'predictionWeightedMean' 'predictionMedian' 'predictionWeightedMedian'});
                    predictions       = zeros(targetData.nModels,length(trainingConfigurations));
                    weights           = zeros(targetData.nModels,length(trainingConfigurations));
                    for iConfig = length(trainingConfigurations):-1:1
                        sf                                   = ScoreFunction(trainingConfigurations{iConfig},...
                                                                             targetData,trainingConfigurations{iConfig}.exponentValues);                                    
                        tempPredictions                      = sf.calcScoreTarget();
                        tempPredictions(tempPredictions<0)   = 0;
                        tempPredictions                      = trainingConfigurations{iConfig}.sigmoid.backSigmoid(tempPredictions);
                        predictions(:,iConfig)               = tempPredictions; 
                        %The weight of a prediction is the derivative of the
                        %objectExponenet at that position
                        weights(:,iConfig)                   = trainingConfigurations{iConfig}.sigmoid.derivative(tempPredictions);
                    end
                    targetDataValidation.statistics.predictions = predictions;
                    targetDataValidation.statistics.weights     = weights;
 
                    weightedMean =  sum(predictions.*weights,2)./sum(weights,2);
                    weightedMean(weightedMean>1) = 1;
                    weightedMean(weightedMean<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values weightedMean]; %weighted mean
                    
		    mdian =  median(predictions,2);
                    mdian(mdian>1) = 1;
                    mdian(mdian<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values mdian]; 
                    
                    [weightedMdian, weightedMdianIDR, weightedMdianENT] =  Configuration.weightedMedian(predictions,weights);
                    weightedMdian(weightedMdian>1) = 1;
                    weightedMdian(weightedMdian<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values weightedMdian  weightedMdianIDR weightedMdianENT]; 
                    
                    %exponent weights
                   

                    for iWeight = 1:length(weights(:,1))
                       exp_weights(iWeight,:)=exp(weightsExp*weights(iWeight,:));
                       
                       exp_weights(iWeight,:) = exp_weights(iWeight,:)/sum(exp_weights(iWeight,:));
                       
                    end

                    [weightedEMdian, weightedEMdianIDR, weightedEMdianENT] =  Configuration.weightedMedian(predictions,exp_weights);
                    weightedEMdian(weightedEMdian>1) = 1;
                    weightedEMdian(weightedEMdian<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values weightedEMdian ];

                    
                    %max weighted predictor
                    for iModel = 1:targetDataValidation.nModels
                        sorted_weights = sort(targetDataValidation.statistics.weights(iModel,:));
                        top_weights = sorted_weights(length(sorted_weights)-floor(length(sorted_weights)*0.1));
                        m1 = find(targetDataValidation.statistics.weights(iModel,:)>1.5);
						if (length(m1) < 1)
							m1 = 1:length(targetDataValidation.statistics.weights(iModel,:));
						end
						mPredictors(iModel) = mean(targetDataValidation.statistics.predictions(iModel,m1));

                    end

                    for iModel = 1:targetDataValidation.nModels
                        sorted_weights = sort(targetDataValidation.statistics.weights(iModel,:));
                        top_weights = sorted_weights(length(sorted_weights)-floor(length(sorted_weights)*0.1));
                        m1 = find(targetDataValidation.statistics.weights(iModel,:)>top_weights);
                
                        mPredictors(iModel) = mean(targetDataValidation.statistics.predictions(iModel,m1));
                    end

                    targetDataValidation.values = [targetDataValidation.values mPredictors' ];
                                              
                    if (~isempty(targetData.fieldsComparisonTables))
                        %ClusterPrediction1 - avg the predictions of each
                        %predictor as the score for the model. Then, do
                        %weighted medians for the newPredictions
                        for iModel = 1:targetData.nModels
                            clusterIndices = find(targetData.comparisonTables{2}.CMTable(iModel,:)>0.98);
                            [w l] = size(predictions(clusterIndices,:));
                            
                            if (w~=1)
                            
                                newPredictions(iModel,:) = mean(predictions(clusterIndices,:));
                            else
                                newPredictions(iModel,:) = predictions(clusterIndices,:);
                            end
                        end

                        [weightedClusterMdian, weightedClusterMdianIDR, weightedClusterMdianENT] =  Configuration.weightedMedian(newPredictions,weights);
                        weightedEMdian(weightedEMdian>1) = 1;
                        weightedEMdian(weightedEMdian<0.1) = 0.1;
                        targetDataValidation.values = [targetDataValidation.values weightedClusterMdian ];

				mPredictors = mPredictors';
                for iModel = 1:targetData.nModels
                    clusterIndices = find(targetData.comparisonTables{2}.CMTable(iModel,:)>0.95)';
						ent = [];
						for k = 1:length(clusterIndices)
							ent(k) = entropy(predictions(clusterIndices(k),:)');
						end
						model_entropy = entropy(predictions(iModel,:));
						model_score = weightedMdian(iModel);
						model_weight = (1/(model_entropy/10 + 1));
						clusters_weights = 0;
						try
							clusters_weights = (1./(ent' .* iqr(weights(clusterIndices,:)')').*(median(weights(clusterIndices,:)')').^2);
						catch e
							disp(e)
						end
						clusters_weights = (1./(clusters_weights+1));
						clusters_weights = clusters_weights/sum(clusters_weights);
						weightedMdianCluster(iModel) = sum(mPredictors(clusterIndices).*clusters_weights);
						if (isnan(weightedMdianCluster(iModel)))
							weightedMdianCluster(iModel) = model_score;
						else
							weightedMdianCluster(iModel) = (model_weight*model_score)+(1-model_weight)*weightedMdianCluster(iModel);
						end

                end
                targetDataValidation.values = [targetDataValidation.values weightedMdianCluster' ];
			


			[rows cols] = size(predictions);
			predictions2 = predictions(randperm(rows),randperm(cols));
			weights2 = weights(randperm(rows),randperm(cols));
                   [weightedMdian, weightedMdianIDR, weightedMdianENT] =  Configuration.weightedMedian(predictions2,weights2);
                    weightedMdian(weightedMdian>1) = 1;
                    weightedMdian(weightedMdian<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values weightedMdian ]; 
		         	
			
                    end
               end
        end
        function targetDataValidation = validationTarget(experimentData,  iTarget, configurationArray,weightsExp)
            targetData           = experimentData.targetsData{iTarget};
            targetDataValidation = TargetData(targetData.targetName, targetData.fileNames,...
                                              [], ... % experimentData
                                              {targetData.fields{:} 'predictionWeightedMean' 'predictionMedian' 'predictionWeightedMedian' ...
									'weightedMdianIDR' 'weightedMdianENT' 'predictionExpWeightedMedian' ...
									'predictionMaxWeightedMedian' 'predictionClusterWeightedMedian' 'predictionWeightedMedianCluster' ...
									'predictionRandomWeightedMedian'},...
                                              targetData.values,...
                                              targetData.comparisonTables,...
                                              targetData.fieldsComparisonTables);
            trainingConfigurations = Configuration.extractTrainingConfiguration(targetData,configurationArray);
            disp(['# of training configurations ' int2str(length(trainingConfigurations))]);
                if (isempty(trainingConfigurations))
                    targetDataValidation = [];
                else
                    predictionIndices = Util.getIndices(targetDataValidation.fields,{'predictionWeightedMean' 'predictionMedian' 'predictionWeightedMedian'});
                    objectiveIndex    = Util.getIndices(targetDataValidation.fields,{configurationArray{1}.objectiveFunction});
                    % objectiveDeltaIndex = Util.getIndices(targetDataValidation.fields,{'delta_gdt_ts'});%todo fix this
                    objectiveDeltaIndex = Util.getIndices(targetDataValidation.fields,{configurationArray{1}.deltaObjectiveFunction});%todo fix this
                    predictions       = zeros(targetData.nModels,length(trainingConfigurations));
                    weights           = zeros(targetData.nModels,length(trainingConfigurations));
                    for iConfig = length(trainingConfigurations):-1:1
                        sf                                   = ScoreFunction(trainingConfigurations{iConfig},...
                                                                             targetData,trainingConfigurations{iConfig}.exponentValues);                                    
                        tempPredictions                      = sf.calcScoreTarget();
                        tempPredictions(tempPredictions<0)   = 0;
                        tempPredictions                      = trainingConfigurations{iConfig}.sigmoid.backSigmoid(tempPredictions);
                        predictions(:,iConfig)               = tempPredictions; 
                        %The weight of a prediction is the derivative of the
                        %objectExponenet at that position
                        weights(:,iConfig)                   = trainingConfigurations{iConfig}.sigmoid.derivative(tempPredictions);
                    end
                    targetDataValidation.statistics.predictions = predictions;
                    targetDataValidation.statistics.weights     = weights;
 
                    weightedMean =  sum(predictions.*weights,2)./sum(weights,2);
                    weightedMean(weightedMean>1) = 1;
                    weightedMean(weightedMean<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values weightedMean]; %weighted mean
                    [targetDataValidation.statistics.errorMean, ...
                     targetDataValidation.statistics.corrMean, ...
                     targetDataValidation.statistics.enrichmentMean]   = Util.getErrorAndCorr(weightedMean,...
                                                                                              targetDataValidation.values(:,objectiveIndex));
                    [targetDataValidation.statistics.origErrorMean, ...
                     targetDataValidation.statistics.origCorrMean, ...
                     targetDataValidation.statistics.origEnrichmentMean]   = Util.getErrorAndCorr(weightedMean,...
                                                                                              targetDataValidation.values(:,objectiveIndex)-...
                                                                                              targetDataValidation.values(:,objectiveDeltaIndex));
                    mdian =  median(predictions,2);
                    mdian(mdian>1) = 1;
                    mdian(mdian<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values mdian]; 
                    [targetDataValidation.statistics.errorMedian, ...
                     targetDataValidation.statistics.corrMedian, ...
                     targetDataValidation.statistics.enrichmentMedian]   = Util.getErrorAndCorr(mdian,...
                                                                                              targetDataValidation.values(:,objectiveIndex));
                    [targetDataValidation.statistics.origErrorMedian, ...
                     targetDataValidation.statistics.origCorrMedian, ...
                     targetDataValidation.statistics.origEnrichmentMedian]   = Util.getErrorAndCorr(mdian,...
                                                                                              targetDataValidation.values(:,objectiveIndex)-...
                                                                                              targetDataValidation.values(:,objectiveDeltaIndex));
                    
                    [weightedMdian, weightedMdianIDR, weightedMdianENT] =  Configuration.weightedMedian(predictions,weights);
                    weightedMdian(weightedMdian>1) = 1;
                    weightedMdian(weightedMdian<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values weightedMdian  weightedMdianIDR weightedMdianENT]; 
                    [targetDataValidation.statistics.errorWeightedMedian, ...
                     targetDataValidation.statistics.corrWeightedMedian, ...
                     targetDataValidation.statistics.enrichmentWeightedMedian]   = Util.getErrorAndCorr(weightedMdian,...
                                                                                              targetDataValidation.values(:,objectiveIndex));
                    [targetDataValidation.statistics.origErrorWeightedMedian, ...
                     targetDataValidation.statistics.origCorrWeightedMedian, ...
                     targetDataValidation.statistics.origEnrichmentWeightedMedian]   = Util.getErrorAndCorr(weightedMdian,...
                                                                                              targetDataValidation.values(:,objectiveIndex)-...
                                                                                              targetDataValidation.values(:,objectiveDeltaIndex));
                    
                    %exponent weights
                   

                    for iWeight = 1:length(weights(:,1))
                       exp_weights(iWeight,:)=exp(weightsExp*weights(iWeight,:));
                       
                       exp_weights(iWeight,:) = exp_weights(iWeight,:)/sum(exp_weights(iWeight,:));
                       
                    end

                    [weightedEMdian, weightedEMdianIDR, weightedEMdianENT] =  Configuration.weightedMedian(predictions,exp_weights);
                    weightedEMdian(weightedEMdian>1) = 1;
                    weightedEMdian(weightedEMdian<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values weightedEMdian ];
                    [targetDataValidation.statistics.errorExpWeightedMedian, ...
                     targetDataValidation.statistics.corrExpWeightedMedian, ...
                     targetDataValidation.statistics.enrichmentExpWeightedMedian]   = Util.getErrorAndCorr(weightedEMdian,...
                                                                                                                  targetDataValidation.values(:,objectiveIndex));

                     %max weighted predictor
                    %for iModel = 1:targetDataValidation.nModels
                    %    sorted_weights = sort(targetDataValidation.statistics.weights(iModel,:));
                    %    top_weights = sorted_weights(length(sorted_weights)-floor(length(sorted_weights)*0.1));
                    %    m1 = find(targetDataValidation.statistics.weights(iModel,:)>top_weights);

                    %    mPredictions(iModel,:) = targetDataValidation.statistics.predictions(iModel,m1);
                     %   mWeights(iModel,:) = targetDataValidation.statistics.weights(iModel,m1);

                    %end
                    %[weightedEMdian, weightedEMdianIDR, weightedEMdianENT] =  Configuration.weightedMedian(mPredictions,mWeights);
                    %weightedEMdian(weightedEMdian>1) = 1;
                    %weightedEMdian(weightedEMdian<0.1) = 0.1;

                    %[targetDataValidation.statistics.errorExpWeightedMedian, ...
                    % targetDataValidation.statistics.corrExpWeightedMedian, ...
                    % targetDataValidation.statistics.enrichmentExpWeightedMedian]   = Util.getErrorAndCorr(weightedEMdian,...
                    %                                                                          targetDataValidation.values(:,objectiveIndex));
                    
                    %max weighted predictor
                    for iModel = 1:targetDataValidation.nModels
                        sorted_weights = sort(targetDataValidation.statistics.weights(iModel,:));
                        top_weights = sorted_weights(length(sorted_weights)-floor(length(sorted_weights)*0.1));
                        m1 = find(targetDataValidation.statistics.weights(iModel,:)>1.5);
						if (length(m1) < 1)
							m1 = 1:length(targetDataValidation.statistics.weights(iModel,:));
						end
						mPredictors(iModel) = mean(targetDataValidation.statistics.predictions(iModel,m1));

                    end

                    for iModel = 1:targetDataValidation.nModels
                        sorted_weights = sort(targetDataValidation.statistics.weights(iModel,:));
                        top_weights = sorted_weights(length(sorted_weights)-floor(length(sorted_weights)*0.1));
                        m1 = find(targetDataValidation.statistics.weights(iModel,:)>top_weights);
                
                        mPredictors(iModel) = mean(targetDataValidation.statistics.predictions(iModel,m1));
                    end

                    targetDataValidation.values = [targetDataValidation.values mPredictors' ];
                    targetDataValidation.statistics.errorMaxWeightedMedian =  Util.getErrorAndCorr(mPredictors',...
                                                                                              targetDataValidation.values(:,objectiveIndex));
                                              
                    if (~isempty(targetData.fieldsComparisonTables))
                        %ClusterPrediction1 - avg the predictions of each
                        %predictor as the score for the model. Then, do
                        %weighted medians for the newPredictions
                        for iModel = 1:targetData.nModels
                            clusterIndices = find(targetData.comparisonTables{2}.CMTable(iModel,:)>0.98);
                            [w l] = size(predictions(clusterIndices,:));
                            
                            if (w~=1)
                            
                                newPredictions(iModel,:) = mean(predictions(clusterIndices,:));
                            else
                                newPredictions(iModel,:) = predictions(clusterIndices,:);
                            end
                        end

                        [weightedClusterMdian, weightedClusterMdianIDR, weightedClusterMdianENT] =  Configuration.weightedMedian(newPredictions,weights);
                        weightedEMdian(weightedEMdian>1) = 1;
                        weightedEMdian(weightedEMdian<0.1) = 0.1;
                        targetDataValidation.values = [targetDataValidation.values weightedClusterMdian ];
                        [targetDataValidation.statistics.errorClusterWeightedMedian, ...
                         targetDataValidation.statistics.corrClusterWeightedMedian, ...
                         targetDataValidation.statistics.enrichmentClusterWeightedMedian]   = Util.getErrorAndCorr(weightedClusterMdian,...
                                                                                                                      targetDataValidation.values(:,objectiveIndex));
                        %ClusterPrediction2 - avg the weighted median predictions 


%                        for iModel = 1:targetData.nModels
%                            clusterIndices = find(targetData.comparisonTables{2}.CMTable(iModel,:)>0.98);
%                            weightedMdianCluster(iModel) = mean(weightedMdian(clusterIndices));
%                        end

%                        targetDataValidation.values = [targetDataValidation.values weightedMdianCluster' ];
%                        [targetDataValidation.statistics.errorWeightedMedianCluster, ...
%                         targetDataValidation.statistics.corrWeightedMedianCluster, ...
%                         targetDataValidation.statistics.enrichmentWeightedMedianCluster]   = Util.getErrorAndCorr(weightedMdianCluster',...
%                                                                                                                      targetDataValidation.values(:,objectiveIndex));
				mPredictors = mPredictors';
                for iModel = 1:targetData.nModels
                    clusterIndices = find(targetData.comparisonTables{2}.CMTable(iModel,:)>0.95)';
		%	if (length(clusterIndices)>1)
						ent = [];
						for k = 1:length(clusterIndices)
							ent(k) = entropy(predictions(clusterIndices(k),:)');
						end
						model_entropy = entropy(predictions(iModel,:));
						model_score = weightedMdian(iModel);
						model_weight = (1/(model_entropy/10 + 1));
						clusters_weights = (1./(ent' .* iqr(weights(clusterIndices,:)')').*(median(weights(clusterIndices,:)')').^2);
						clusters_weights = (1./(clusters_weights+1));
						clusters_weights = clusters_weights/sum(clusters_weights);
						weightedMdianCluster(iModel) = sum(mPredictors(clusterIndices).*clusters_weights);
						%weightedMdianCluster(iModel) = sum(weightedMdian(clusterIndices).*clusters_weights);
						if (isnan(weightedMdianCluster(iModel)))
							weightedMdianCluster(iModel) = model_score;
						else
							weightedMdianCluster(iModel) = (model_weight*model_score)+(1-model_weight)*weightedMdianCluster(iModel);
						end
		%	else
		%		weightedMdianCluster(iModel) = weightedMdian(iModel);
		%	end
                end
                targetDataValidation.values = [targetDataValidation.values weightedMdianCluster' ];
                        [targetDataValidation.statistics.errorWeightedMedianCluster, ...
                         targetDataValidation.statistics.corrWeightedMedianCluster, ...
                         targetDataValidation.statistics.enrichmentWeightedMedianCluster]   = Util.getErrorAndCorr(weightedMdianCluster',...
                                                                                                                      targetDataValidation.values(:,objectiveIndex));
			


			%random1
			[rows cols] = size(predictions);
			predictions2 = predictions(randperm(rows),randperm(cols));
			weights2 = weights(randperm(rows),randperm(cols));
                   [weightedMdian, weightedMdianIDR, weightedMdianENT] =  Configuration.weightedMedian(predictions2,weights2);
                    weightedMdian(weightedMdian>1) = 1;
                    weightedMdian(weightedMdian<0.1) = 0.1;
                    targetDataValidation.values = [targetDataValidation.values weightedMdian]; 
                    [targetDataValidation.statistics.errorRandomWeightedMedian, ...
                     targetDataValidation.statistics.corrRandomWeightedMedian, ...
                     targetDataValidation.statistics.enrichmentRandomWeightedMedian]   = Util.getErrorAndCorr(weightedMdian,...
                                                                                              targetDataValidation.values(:,objectiveIndex));
                    [targetDataValidation.statistics.origErrorRandomWeightedMedian, ...
                     targetDataValidation.statistics.origCorrRandomWeightedMedian, ...
                     targetDataValidation.statistics.origEnrichmentRandomWeightedMedian]   = Util.getErrorAndCorr(weightedMdian,...
                                                                                              targetDataValidation.values(:,objectiveIndex)-...
                                                                                              targetDataValidation.values(:,objectiveDeltaIndex));
		         	
			
                    end
               end
        end
        
        
            
        function experimentDataValidation = prediction(experimentData, configurationArray,commentsFile,weightsExp)
            experimentDataValidation                    = ExperimentData([experimentData.name '_prediction']);
            experimentDataValidation.copyrightAndReadme = experimentData.copyrightAndReadme;
            %experimentDataValidation.readReadme(commentsFile,length(experimentDataValidation.copyrightAndReadme)+1);
            experimentDataValidation.flagComparisonTable = experimentData.flagComparisonTable;
            for iTarget = length(experimentData.targetsData):-1:1
                experimentDataValidation.targetsData{iTarget} = Configuration.predictionTarget(experimentData, iTarget, configurationArray,weightsExp);
            end
        end
        
        function experimentDataValidation = validation(experimentData, configurationArray,commentsFile,weightsExp)
            experimentDataValidation                    = ExperimentData([experimentData.name '_validation']);
            experimentDataValidation.copyrightAndReadme = experimentData.copyrightAndReadme;
            %experimentDataValidation.readReadme(commentsFile,length(experimentDataValidation.copyrightAndReadme)+1);
            experimentDataValidation.flagComparisonTable = experimentData.flagComparisonTable;
            for iTarget = length(experimentData.targetsData):-1:1
                experimentDataValidation.targetsData{iTarget} = Configuration.validationTarget(experimentData, iTarget, configurationArray,weightsExp);
            end
        end
           

        function [medians, interdecile, entropies] = weightedMedian(predictions,weights) %matrix rows models; columns predictions or weights theirof.
                [sorted, order] = sort(predictions,2);
                for i = size(predictions,1):-1:1
                    sortedWeights(i,:) = weights(i,order(i,:));
                end
                weightSums              = sum(sortedWeights,2);
                lowDecile               = weightSums/10;
                topDecile               = weightSums*9/10;
                midpoints               = weightSums/2;
                for i = length(weightSums):-1:1;
                    bin = [];
                    for j = 100:-1:1
                        bin(j) = sum(weights(i,(predictions(i,:)>(j-1)/100) & (predictions(i,:)<=(j)/100)))/weightSums(i);
                    end
                    p    = bin(bin>0);
                    entropies(i,1) = -sum(p.*log(p),2);
                    if (entropies(i,1)<-0.001)
                        sum(bin)
                        entropies(i,1)
                        error('this is weird');
                    end
                end
                cumulativeSortedWeights = cumsum(sortedWeights,2);
                for i = size(predictions,1):-1:1
                    lowIndices(i,1) = find(cumulativeSortedWeights(i,:)>lowDecile(i,1),1,'first');
                    topIndices(i,1) = find(cumulativeSortedWeights(i,:)>topDecile(i,1),1,'first');
                    indices(i,1)    = find(cumulativeSortedWeights(i,:)>midpoints(i,1),1,'first');
                end
                for i = size(predictions,1):-1:1
                    top          = sorted(i,indices(i,1));
                    bottomWeight = cumulativeSortedWeights(i,indices(i,1))-midpoints(i,1);
                    if (indices(i,1) == 1)
                        bottom = 0;
                        topWeight = midpoints(i,1);
                    else
                        bottom       = sorted(i,indices(i,1)-1);
                        topWeight    = midpoints(i,1)-cumulativeSortedWeights(i,indices(i,1)-1);
                    end
                    medians(i,1) = (top.*topWeight+bottom.*bottomWeight)./(topWeight+bottomWeight);
                    interdecile(i,1)  = sorted(i,topIndices(i,1)) - sorted(i,lowIndices(i,1));
                end
        end
 
        function allConfigurations = collectConfigurationArrays()
            d =  dir('configurationArray.*');
            fileNames = {d.name};
            allConfigurations = {};
            for i = length(fileNames):-1:1
		disp(fileNames{i})
                load(fileNames{i});
                allConfigurations = {allConfigurations{:} configurationArray{:}};
            end
        end
        
        function fieldHistogram(configurationArray)
            h = zeros(1,configurationArray{1}.nFields);
            for i = 1:length(configurationArray)
                h = h + (configurationArray{i}.coefs' ~= 0)*1;
            end
            hist(h);
            [sh is] = sort(h);
            for i = 1:length(is)
                fprintf(1,'%d %s\n',h(is(i)),configurationArray{1}.fields{is(i)});
            end
        end
        
        function writeXML(configurationArray,fileName)
           rootNode = com.mathworks.xml.XMLUtils.createDocument('configurationArray');
           root     = rootNode.getDocumentElement;
           for iConfig = 1:length(configurationArray)
               iConfig
               config = configurationArray{iConfig};
               configElement = rootNode.createElement('configuration');
               configElement.appendChild(Configuration.getListElementXML('fields',config,rootNode));
               configElement.appendChild(Configuration.getListElementXML('normalizers',config,rootNode));               
               configElement.appendChild(Configuration.getListElementXML('coefs',config,rootNode));  
               configElement.appendChild(Configuration.getListElementXML('exponentIndices',config,rootNode));  
               configElement.appendChild(Configuration.getListElementXML('exponentValues',config,rootNode));  
               configElement.appendChild(Configuration.getListElementXML('objectiveFunction',config,rootNode));
               
               sigmoidElement = rootNode.createElement('sigmoid');
               sigmoidElement.setAttribute('alpha',num2str(config.sigmoid.alpha));
               sigmoidElement.setAttribute('beta',num2str(config.sigmoid.beta));
               sigmoidElement.setAttribute('gamma',num2str(config.sigmoid.gamma));
               sigmoidElement.setAttribute('delta',num2str(config.sigmoid.delta));
               configElement.appendChild(sigmoidElement);
               
               root.appendChild(configElement);
           end
           xmlwrite(fileName,rootNode);
        end
        
        function element = getListElementXML(elementName, config, rootNode)
               element     = rootNode.createElement(elementName);
               elementList = config.(elementName);
               if (iscell(elementList))
                for iItem = length(elementList):-1:1
                    newList{(iItem-1)*2+1} = elementList{iItem};
                    newList{(iItem-1)*2+2} = '  ';
                end
                element.setAttribute(elementName,[newList{:}]);
               elseif (ischar(elementList))
                       element.setAttribute(elementName,elementList);
               else
                [n, m] = size(elementList);
                flattened = zeros(1,n*m);
                for i = 1:n
                    for j = 1:m
                        flattened((i-1)*m+j) = elementList(i,j);
                    end
                end
                element.setAttribute('dim',num2str(size(elementList)));
                element.setAttribute('comment','flattened((i-1)*m+j) = elementList(i,j);');
                element.setAttribute(elementName,num2str(flattened));
               end
        end
                
    end
        
end
