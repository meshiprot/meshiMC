classdef CoefOptimization < handle
    properties(Constant=true)
        defaultObjectiveFunction = 'gdt_ts';%'contacts6_mcc';
	defaultDeltaObjective = 'delta_gdt_ts';  
%        DefaultFeaturesToLearn  = {   'energy' 'distanceConstraints'  'bondEnergy'  'angleEnergy' 'planeEnergy'  'outOfPlaneEnergy'		 ...
%    'ramachandranSidechain'    'ramachSTD'    'ramachandranSidechainCore'    'cooperativeZRamachandranSidechain'                                 ...
%    'cooperativeZstdRamachandranSidechain'    'compositePropensity'    'cooperativeZPropensity'    'cooperativeZstdPropensity'                   ...
%    'atomicPairwisePMFSumma'    'summaStd'    'excludedVolume'    'cooperativeZSumma'                                                            ...
%    'cooperativeSummaPolar'    'cooperativeSummaNonPolar'    'cooperativeSummaNeutral'    'cooperativeSummaPolarNN_OO'                           ...
%    'cooperativeSummaPolarBb'    'solvationEnergy'    'solvationSCpolar'    'solvationSCcarbon'                                                  ...
%    'solvationBBpolar'    'solvationBBcarbon'    'solvationHB'    'solvationBuriedHB'                                                            ...
%    'solvationSTD'    'solvationEntropy'    'bbPolarN'    'bbCarbonN'    'scPolarN'    'scCarbonN'                                               ...
%    'hydrogenBonds'    'hydrogenBondsPairs'    'hydrogenBondsAnglesHOC'    'hydrogenBondsAnglesOHN'    'contactsAndSASA'                         ...
%    'contacts12'    'contacts14'    'contacts14Core'    'SASAratio'    'conservationContacts8'    'contacts8'                                    ...
%    'conservationRgRatio'    'conservation_H_RgRatio'    'conservationContacts11'    'contacts11'    'conservationContacts15'    'contacts15'    ...
%    'conservationContactsHr'    'contactsHr'    'conservationRgRatioHr'    'one'    'tetherEnergy'    'flatRamachEnergy'                         ...
%    'tetherAll'    'samudralaEnergy'    'secondaryStructureFraction' 'length'};
%	 DefaultFeaturesToLearn (Thesis) = {'length' ,'goap','goap_ag','dfire',...
%	    'energy' ...
%           'bondEnergy' 'angleEnergy' 'planeEnergy'  'outOfPlaneEnergy' ...
%            'ramachandranSidechain' 'ramachandranSidechainCore' 'ramachSTD' 'cooperativeZRamachandranSidechain' ... 
%            'cooperativeZstdRamachandranSidechain'...
%            'compositePropensity' 'cooperativeZPropensity' 'cooperativeZstdPropensity' ...
%            'atomicPairwisePMFSumma' 'summaStd' 'excludedVolume' ...
%            'cooperativeZSumma' 'cooperativeSummaPolar' 'cooperativeSummaNonPolar' ...
%            'cooperativeSummaNeutral' 'cooperativeSummaPolarNN_OO' 'cooperativeSummaPolarBb' ...
%            'cooperativeZstdSumma' 'cooperativeStdSummaPolar' 'cooperativeStdSummaNonPolar' ...
%            'cooperativeStdSummaNeutral' 'cooperativeStdSummaPolarNN_OO'  'cooperativeStdSummaPolarBb' ...
%            'solvateEnergy' 'solvationSCpolar' 'solvationSCcarbon' ...%'solvationBBpolar' 
%            'solvationSTD' ...
%            'hydrogenBonds' 'hydrogenBondsPairs' 'hydrogenBondsAnglesHOC' 'hydrogenBondsAnglesOHN'...
%            'rg' 'N_RGhSS' 'hSS' 'bSS' 'cSS' ... 
%            'hSShCoil' 'hSSbSS' 'hSSbCoil' 'hSScSS' 'hSScCoil'...
%            'contacts8'  'contacts11' 'contacts15' 'contactsHr'...
%            'flatRamachEnergy' 'samudralaEnergy' 'one'};
%DefaultFeaturesToLearn = {'optimizationScore_interdecile'...
%				'optimizationScore_weightedMedianScore'  'score1_interdecile'  'score1_weightedMedianScore'...
%				'score2_interdecile'   'score2_weightedMedianScore'  ...
%				'energy'...
%				'bondEnergy'    'angleEnergy'    'planeEnergy'...
%				'ramachandranSidechain'    'ramachSTD'    'cooperativeZRamachandranSidechain'...
%				'cooperativeZstdRamachandranSidechain'    'ramachandranCore'    'compositePropensity'...
%				'cooperativeZPropensity'    'cooperativeZstdPropensity'    'atomicPairwisePMFSumma'...
%				'summaStd'    'excludedVolume'    'cooperativeZSumma'    'cooperativeSummaPolar'...
%				'cooperativeSummaNonPolar'    'cooperativeSummaNeutral'    'cooperativeSummaPolarNN_OO'...
%				'cooperativeSummaPolarBb'    'cooperativeZstdSumma'    'cooperativeStdSummaPolar'...
%				'cooperativeStdSummaNonPolar'    'cooperativeStdSummaNeutral'    'cooperativeStdSummaPolarNN_OO'...
%				'cooperativeStdSummaPolarBb'    'solvationEnergy'    'solvationSCpolar'    'solvationSCcarbon'...
%				'solvationBBpolar'    'solvationBBcarbon'    'solvationHB'    'solvationBuriedHB'...
%				'solvationSTD'    'solvationEntropy'    ...
%				'atomEnvironmentEnergy'    'atomEnvironmentPropensity'    'atomEnvironmentEnergySS'...
%%				'atomEnvironmentPropensitySS'    'solvateEnergy'    'solvateSCpolar'...
%				'solvateSCcarbon'    'solvateBBpolar'    'solvateSTD'    'hydrogenBonds'...
%				'hydrogenBondsPairs'    'hydrogenBondsAnglesHOC'    'hydrogenBondsAnglesOHN'...
%%				'contactsAndSASA'    'contacts12'    'contacts14'    'contacts14Core'    'SASAratio'    'rg'...
%%				'N_RGhSS'  'hSS'...
%				'bSS'    'cSS'    'hSShCoil'    'hSSbSS'   'hSSbCoil'    'hSScSS'    'hSScCoil'    'goap'...
%				'dfire'    'goap_ag'  'contacts8'  ...
%%				'contacts11' ...
%				'contacts15'    'contactsHr' ...
%				'flatRamachEnergy'   'samudralaEnergy'    'secondaryStructureFraction'...
%				'ssCompatibility'    'sasaCompatibility'    'coverage'    'strictCoverage' ...
%				    'scwrl'    'one' 'length' 'nAtoms'}; %CASP12 score function training	
%DefaultFeaturesToLearn = {	'energy'...
%				'bondEnergy'    'angleEnergy'    'planeEnergy'...
%				'ramachandranSidechain'    'ramachSTD'    'cooperativeZRamachandranSidechain'...
%%				'cooperativeZstdRamachandranSidechain'    'ramachandranCore'    'compositePropensity'...
%				'cooperativeZPropensity'    'cooperativeZstdPropensity'    'atomicPairwisePMFSumma'...
%				'summaStd'    'excludedVolume'    'cooperativeZSumma'    'cooperativeSummaPolar'...
%				'cooperativeSummaNonPolar'    'cooperativeSummaNeutral'    'cooperativeSummaPolarNN_OO'...
%				'cooperativeSummaPolarBb'    'cooperativeZstdSumma'    'cooperativeStdSummaPolar'...
%%				'cooperativeStdSummaNonPolar'    'cooperativeStdSummaNeutral'    'cooperativeStdSummaPolarNN_OO'...
%				'cooperativeStdSummaPolarBb'    'solvationEnergy'    'solvationSCpolar'    'solvationSCcarbon'...
%				'solvationBBpolar'    'solvationBBcarbon'    'solvationHB'    'solvationBuriedHB'...
%%				'solvationSTD'    'solvationEntropy'    ...
%				'atomEnvironmentEnergy'    'atomEnvironmentPropensity'    'atomEnvironmentEnergySS'...
%				'atomEnvironmentPropensitySS'    'solvateEnergy'    'solvateSCpolar'...
%				'solvateSCcarbon'    'solvateBBpolar'    'solvateSTD'    'hydrogenBonds'...
%				'hydrogenBondsPairs'    'hydrogenBondsAnglesHOC'    'hydrogenBondsAnglesOHN'...
%				'contactsAndSASA'    'contacts12'    'contacts14'    'contacts14Core'    'SASAratio'    'rg'...
%				'N_RGhSS'  'hSS'...
%				'bSS'    'cSS'    'hSShCoil'    'hSSbSS'   'hSSbCoil'    'hSScSS'    'hSScCoil'    'goap'...
%				'dfire'    'goap_ag'  'contacts8'  ...
%				'contacts11' ...
%				'contacts15'    'contactsHr' ...
%				'flatRamachEnergy'   'samudralaEnergy'    'secondaryStructureFraction'...
%				'ssCompatibility'    'sasaCompatibility'    'coverage'    'strictCoverage' ...
%				    'scwrl'    'one' 'length' 'nAtoms'}; %CASP12 score function training, with no iterative functions
DefaultFeaturesToLearn = {	'energy'...
				'bondEnergy'    'angleEnergy'    'planeEnergy'...
				'ramachandranSidechain'    'ramachSTD'    'cooperativeZRamachandranSidechain'...
				'cooperativeZstdRamachandranSidechain'    'ramachandranCore'    'compositePropensity'...
				'cooperativeZPropensity'    'cooperativeZstdPropensity'    'atomicPairwisePMFSumma'...
				'summaStd'    'excludedVolume'    'cooperativeZSumma'    'cooperativeSummaPolar'...
				'cooperativeSummaNonPolar'    'cooperativeSummaNeutral'    'cooperativeSummaPolarNN_OO'...
				'cooperativeSummaPolarBb'    'cooperativeZstdSumma'    'cooperativeStdSummaPolar'...
				'cooperativeStdSummaNonPolar'    'cooperativeStdSummaNeutral'    'cooperativeStdSummaPolarNN_OO'...
				'cooperativeStdSummaPolarBb'    'solvationEnergy'    'solvationSCpolar'    'solvationSCcarbon'...
				'solvationBBpolar'    'solvationBBcarbon'    'solvationHB'    'solvationBuriedHB'...
				'solvationSTD'    'solvationEntropy'    ...
				'atomEnvironmentEnergy'    'atomEnvironmentPropensity'    'atomEnvironmentEnergySS'...
				'atomEnvironmentPropensitySS'    'solvateEnergy'    'solvateSCpolar'...
				'solvateSCcarbon'    'solvateBBpolar'    'solvateSTD'    'hydrogenBonds'...
				'hydrogenBondsPairs'    'hydrogenBondsAnglesHOC'    'hydrogenBondsAnglesOHN'...
				'contactsAndSASA'    'contacts12'    'contacts14'    'contacts14Core'    'SASAratio'    'rg'...
				'N_RGhSS'  'hSS'...
				'bSS'    'cSS'    'hSShCoil'    'hSSbSS'   'hSSbCoil'    'hSScSS'    'hSScCoil'    'goap'...
				'dfire'    'goap_ag'  'contacts8'  ...
				'contacts11' ...
				'contacts15'    'contactsHr' ...
				'flatRamachEnergy'   'samudralaEnergy'    'secondaryStructureFraction'...
				'ssCompatibility'    'sasaCompatibility'    'coverage'    'strictCoverage' ...
				    'scwrl'    'one' 'length' 'nAtoms' 'score1' 'interdcile1'};	% 2016 paper - ss, sasa, Itrerative model
        defaultLearningFraction = 0.5;
        maxNumberOfModels = 50;
%        configEnergy     = ConfigEnergyErrors;
%        configEnergy     = ConfigEnergyMedians;
        defaultConfigEnergy     = ConfigEnergyMaxScoreWeightedSumErrorsAndMedianCorrelation;
        defaultExponentValues  = [1 0.5000 0 -0.5000 -1 -2];
        reportEvery = 50;
        defaultNumberOfCoefsRange = 5;
        defaultNumberOfNormalizers = 2;
        excludedList = {'T0629-D2'};
	featuresToReduce = {'contacts8'}; %features that tend to be over represented (good features but reduce diversity)
        %featuresToReduce = {'change' 'contacts8'}; %features that tend to be over represented (good features but reduce diversity)
    end
    properties
        origObjective
        objective
	deltaObjective
        features
        energy;
        maxSteps
        exponentValues
        featuresToLearn
        randomStream
        objectiveFunction
        learningFraction
        configEnergy
	    numberOfCoefsRange
        experimentData;
        numberOfNormalizers
        commentsFile
    end
    methods
        function obj = CoefOptimization(experimentData,commentsFile)
            obj.featuresToLearn     = CoefOptimization.DefaultFeaturesToLearn;
	    obj.deltaObjective      = CoefOptimization.defaultDeltaObjective;
            obj.objectiveFunction   = CoefOptimization.defaultObjectiveFunction;
            obj.exponentValues      = CoefOptimization.defaultExponentValues;
            obj.learningFraction    = CoefOptimization.defaultLearningFraction;
            obj.configEnergy        = CoefOptimization.defaultConfigEnergy;
	        obj.numberOfCoefsRange  = CoefOptimization.defaultNumberOfCoefsRange;
            obj.experimentData      = experimentData;
            obj.numberOfNormalizers = CoefOptimization.defaultNumberOfNormalizers;
            obj.commentsFile        = commentsFile;
        end
        
        function setExponentValues(obj,exponentValues)
            obj.exponentValues = exponentValues;
        end
        
        function setFeaturesToLearn(obj,featuresToLearn)
            obj.featuresToLearn = featuresToLearn
            obj.features   = obj.experimentData.extract(obj.featuresToLearn,[obj.experimentData.name '_features']);
        end
        %        function obj = Configuration(fields, normalizers,coefs, exponentIndices, ...
        %                                        nTargets, learningFraction,numberOfCoefs,randomStream)
        
        function [configurationArray, ...
                  validationData] = repeatedCoefOptimization(obj, normalizers, ...
                                                                    numberOfCoefs,...
                                                                    numberOfRepeats,...
                                                                    maxSteps,temperature,trainingSetFraction,seed,weightsExp)

            obj.randomStream  = RandStream('mcg16807','Seed',seed);
            obj.origObjective      = obj.experimentData.extract({obj.objectiveFunction},[obj.experimentData.name '_' obj.objectiveFunction]);
            obj.features           = obj.experimentData.extract(obj.featuresToLearn,[obj.experimentData.name '_features']);
            trainingSet = obj.trainingSet(trainingSetFraction,obj.randomStream);
            testSet     = obj.experimentData.extractValidationSet(trainingSet);
%            testSet     = obj.experimentData.extractTargets(trainingSet,'debug');
            nNormalizers    = length(normalizers);
            origNormalizers = normalizers;
            origNumberOfCoefs = numberOfCoefs;
            for i = numberOfRepeats:-1:1
                numberOfCoefs = origNumberOfCoefs;%%fixed number of Coefs
                %numberOfCoefs = origNumberOfCoefs+obj.randomStream.randi(obj.numberOfCoefsRange);
                randomIndices = obj.randomStream.randperm(nNormalizers-1);
                if (obj.numberOfNormalizers == 1)
                    normalizers = {normalizers{1}};
                else
                    normalizers = {normalizers{1} origNormalizers{randomIndices(obj.numberOfNormalizers-1)+1}};
                end
                %%disp(['round # ' num2str(i)]);
%                try
                    configurationArray{i} = obj.coefOptimization(normalizers, ...
                                                                 numberOfCoefs,...
                                                                 maxSteps,temperature,trainingSet,[],numberOfRepeats,i);
%                catch ex
%                    disp(ex.message)
%                    configurationArray{i} = [];
%                    ex.throw;
%                end
                validationData = CoefOptimization.validate(testSet, configurationArray(i:numberOfRepeats),obj.commentsFile,weightsExp);
            end
        end
        
        
        
        
%       repeatedCoefOpt_NoTrainSetPatrtition 
%
%       This function is a version of the function
%       'repeatedCoefOptimization' only without partitioning the training
%       set in the this object, but, receiving the training set in the
%       args.

        function configurationArray = repeatedCoefOpt_NoTrainSetPatrtition(obj, normalizers, ...
                                                                    numberOfCoefs,...
                                                                    numberOfRepeats,...
                                                                    maxSteps,temperature,seed)

            obj.randomStream  = RandStream('mcg16807','Seed',seed);
            obj.origObjective      = obj.experimentData.extract({obj.objectiveFunction},[obj.experimentData.name '_' obj.objectiveFunction]);
            obj.features           = obj.experimentData.extract(obj.featuresToLearn,[obj.experimentData.name '_features']);

            trainingSet = obj.experimentData.getTargetNames;
            
            
            nNormalizers    = length(normalizers);
            origNormalizers = normalizers;
            origNumberOfCoefs = numberOfCoefs;
            for i = numberOfRepeats:-1:1
                numberOfCoefs = origNumberOfCoefs; %%fixed number of Coefs.
                %numberOfCoefs = origNumberOfCoefs+obj.randomStream.randi(obj.numberOfCoefsRange);
                randomIndices = obj.randomStream.randperm(nNormalizers-1);
                if (obj.numberOfNormalizers == 1)
                    normalizers = {normalizers{1}};
                else
                    normalizers = {normalizers{1} origNormalizers{randomIndices(obj.numberOfNormalizers-1)+1}};
                end
                %%disp(['round # ' num2str(i)]);
 %               try

                    configurationArray{i} = obj.coefOptimization(normalizers, ...
                                                                 numberOfCoefs,...
                                                                 maxSteps,temperature,trainingSet,[],numberOfRepeats,i);

%                catch ex
%                    disp(ex.message)
%                    configurationArray{i} = [];
%                    ex.throw;
%                end

            end
        end
        
        
        function configurationArray = repeatedCoefOpt_ParrallelOptimization(obj, normalizers, ...
                                                                    numberOfCoefs,...
                                                                    numberOfRepeats,...
                                                                    maxSteps,temperature,seed,...
                                                                    numberOfPredictors,predictorNumber)

            obj.randomStream  = RandStream('mcg16807','Seed',seed);
            obj.origObjective      = obj.experimentData.extract({obj.objectiveFunction},[obj.experimentData.name '_' obj.objectiveFunction]);
            obj.features           = obj.experimentData.extract(obj.featuresToLearn,[obj.experimentData.name '_features']);

            trainingSet = obj.experimentData.getTargetNames;
            
            
            nNormalizers    = length(normalizers);
            origNormalizers = normalizers;
            origNumberOfCoefs = numberOfCoefs;
            for i = numberOfRepeats:-1:1
                numberOfCoefs = origNumberOfCoefs; %%fixed number of Coefs.
                %numberOfCoefs = origNumberOfCoefs+obj.randomStream.randi(obj.numberOfCoefsRange);
                randomIndices = obj.randomStream.randperm(nNormalizers-1);
                if (obj.numberOfNormalizers == 1)
                    normalizers = {normalizers{1}};
                else
                    normalizers = {normalizers{1} origNormalizers{randomIndices(obj.numberOfNormalizers-1)+1}};
                end
                configurationArray{i} = obj.coefOptimization(normalizers, ...
                                                             numberOfCoefs,...
                                                             maxSteps,temperature,trainingSet,[],numberOfPredictors,predictorNumber);


            end
        end
            
            
        function [bestConfiguration, lastConfiguration] = coefOptimization(obj, normalizers, ...
                                                               numberOfCoefs,maxSteps,temperature,...
                                                               trainingSet, seed,numberOfRepeats,sigmoidNumber)
            if (~isempty(seed))
                obj.randomStream  = RandStream('mcg16807','Seed',seed);
            end
            
            warning('OFF','MATLAB:rankDeficientMatrix');
            obj.maxSteps  = maxSteps;
            dt            = temperature/(obj.maxSteps+1);
            %sigmoid       = Sigmoid.getSigmoid(numberOfRepeats,sigmoidNumber);
            sigmoid       = Sigmoid.getSigmoid(obj.randomStream);
            obj.objective = obj.origObjective.applySigmoidToGdtTs(sigmoid,obj.objectiveFunction,...
                                                                   [obj.origObjective.name '_sigmoid_' num2str(sigmoid.alpha) '_' num2str(sigmoid.beta)]);
            configuration1 = Configuration(obj.featuresToLearn,normalizers,[],[],[],sigmoid,...
                                           length(obj.features.targetsData),... %#of targets
                                           obj.learningFraction,numberOfCoefs,obj.configEnergy,obj,trainingSet,...
                                           obj.objectiveFunction,obj.randomStream);
                configuration1.trajectory = zeros(round(obj.maxSteps/CoefOptimization.reportEvery)+1,3);
                configuration2 = Configuration(obj.featuresToLearn,normalizers,[],[],[],sigmoid,...
                                               length(obj.features.targetsData),... %#of targets
                                               obj.learningFraction,numberOfCoefs,obj.configEnergy,obj,trainingSet,...
                                                obj.objectiveFunction,obj.randomStream);
             for i = 1:length(configuration1.scoreFunctions)
                configuration2.scoreFunctions{i} = configuration1.scoreFunctions{i};
             end
            [obj.energy, error, correlation] = configuration1.evaluate(obj,obj.randomStream);
            bestConfiguration = configuration1.duplicate();
            configuration1.trajectory(1,:) = [obj.energy error correlation];
            %%disp(obj.energy);
 
            
            for iStep = 1:obj.maxSteps
                if (mod(iStep,CoefOptimization.reportEvery) == 0)
                    fprintf(1,'%7.4f %7.4f %7.4f; ',obj.energy, error, correlation);
                    configuration1.trajectory(iStep/CoefOptimization.reportEvery+1,:) = [obj.energy error correlation];
                end
                if (mod(iStep,200) == 0)
                    fprintf(1,'| %7.4f %7.4f %7.4f \n',bestConfiguration.energy,...
                                                       median(bestConfiguration.errors),...
                                                       median(bestConfiguration.correlations));
                end
                configuration1.mutate(configuration2,obj.randomStream,0);
                [newEnergy,  newError, newCorrelation]    = configuration2.evaluate(obj,obj.randomStream);
                de = newEnergy - obj.energy;
                if (exp(-de/temperature) >= obj.randomStream.rand())
                    temp = configuration2;
                    configuration2 = configuration1;
                    configuration1 = temp;
                    obj.energy = newEnergy;
                    error  = newError;
                    correlation = newCorrelation;
                    if (newEnergy < bestConfiguration.energy)
                       bestConfiguration = configuration1.duplicate();
                    end 
                else
                    configuration1.changedFields = configuration2.changedFields;
                    for i = 1:length(configuration1.scoreFunctions)
                        configuration1.scoreFunctions{i}.configuration = configuration1;
                    end
                    configuration1.getWeights();
%                    configuration1.getValuesAndObjectives(obj)
                end
                temperature = temperature-dt;
            end
            lastConfiguration = configuration1;
        end
        
        
        
%         function featuresToLearn = getFeaturesToLearn(obj,normalizers)
%             featuresToLearn = obj.featuresToLearn;
%             for i = 1:length(normalizers)
%                 featuresToLearn(strcmp(normalizers{i},featuresToLearn)) = [];
%             end
%             r = randperm(length(featuresToLearn));
%             n = round(obj.fractionOfFieldsToLearn*length(featuresToLearn));
%             featuresToLearn = featuresToLearn(r(1:n));
%             featuresToLearn = {featuresToLearn{:} normalizers{:}};
%         end
        
        % good for debug
        function compareConfigs(config1,config2,families)    
            disp('learning set')
            d = config1.learningSet - config2.learningSet;
            if(sum(d.^2)==0)
                disp('equals');
            else
                disp(config1.learningSet);
                disp(config2.learningSet);
            end
            
            disp('exponentIndices')
            d = config1.exponentIndices - config2.exponentIndices;
            if(sum(d.^2)==0)
                disp('equals');
            else
                disp(config1.exponentIndices);
                disp(config2.exponentIndices);
            end
            
            disp('coefs')
            d = config1.coefs - config2.coefs;
            if(sum(d.^2)==0)
                disp('equals');
            else
                disp(config1.coefs);
                disp(config2.coefs);
            end
            
            disp('energy')
            d = config1.energy - config2.energy;
            if(sum(d.^2)==0)
                disp('equals');
            else
                disp(config1.energy);
                disp(config2.energy);
            end
            
            disp('learning set weights')
            ls1 = getWeights(families(config1.learningSet),config1.exponentIndices);
            ls2 = getWeights(families(config2.learningSet),config2.exponentIndices);
            sum(ls1(1).weights)
            sum(ls2(1).weights)
            d = ls1(1).weights-ls2(1).weights;
            if(sum(d.^2)==0)
                disp('equals');
            else
                disp(ls1(1).weights);
                disp(ls1(1).weights);
            end
        end
        
        function trainingSet = trainingSet(obj,trainingSetFraction,randomStream)
            targetNames = obj.features.getTargetNames;
            randomIndices = randomStream.randperm(length(targetNames));
            trainingSet = targetNames(randomIndices(1:round(length(targetNames)*trainingSetFraction)));
            for iTarget = 1:length(trainingSet)
                if (~isempty(find(strcmp(trainingSet{iTarget},obj.excludedList),1)))
                    disp([trainingSet{iTarget} ' removed from training set']);
                    trainingSet{iTarget} = [];
                end
            end
            trainingSet = trainingSet(~cellfun('isempty',trainingSet));
            obj.features = obj.features.extractTargets(trainingSet,[obj.features.name 'TrainingSet']);
            obj.origObjective = obj.origObjective.extractTargets(trainingSet,[obj.origObjective.name 'TrainingSet']);
        end
            
    end
    
    methods(Static=true)
       function  predictionData = predict(testSet, configurationArray, commentsFile,weightsExp)
            predictionData = Configuration.prediction(testSet, configurationArray,commentsFile,weightsExp);
        end 

       function [stats validationData] = validate(testSet, configurationArray, commentsFile,weightsExp)
            disp(testSet)
	    validationData = Configuration.validation(testSet, configurationArray,commentsFile,weightsExp);
            disp(validationData)
            fprintf(1,'\n');
            for iTarget = length(validationData.targetsData):-1:1
                errorMedian(iTarget)         = validationData.targetsData{iTarget}.statistics.errorMedian;
                errorWeightedMedian(iTarget) = validationData.targetsData{iTarget}.statistics.errorWeightedMedian;
                errorMean(iTarget)           = validationData.targetsData{iTarget}.statistics.errorMean;
                corrMedian(iTarget)          = validationData.targetsData{iTarget}.statistics.corrMedian;
                corrWeightedMedian(iTarget)  = validationData.targetsData{iTarget}.statistics.corrWeightedMedian;
		corrMean(iTarget)            = validationData.targetsData{iTarget}.statistics.corrMean
                errorExpWeightedMedian(iTarget) = validationData.targetsData{iTarget}.statistics.errorExpWeightedMedian;
                errorMaxWeightedMedian(iTarget) = validationData.targetsData{iTarget}.statistics.errorMaxWeightedMedian;
                if (validationData.flagComparisonTable~=0)
                	errorRandomWeightedMedian(iTarget)            = validationData.targetsData{iTarget}.statistics.errorRandomWeightedMedian;
	                corrRandomWeightedMedian(iTarget)            = validationData.targetsData{iTarget}.statistics.corrRandomWeightedMedian;
                	errorClusterWeightedMedian(iTarget) = validationData.targetsData{iTarget}.statistics.errorClusterWeightedMedian;
               		errorWeightedMedianCluster(iTarget) = validationData.targetsData{iTarget}.statistics.errorWeightedMedianCluster;        
	                corrClusterWeightedMedian(iTarget) = validationData.targetsData{iTarget}.statistics.corrClusterWeightedMedian;
        	        corrWeightedMedianCluster(iTarget) = validationData.targetsData{iTarget}.statistics.corrWeightedMedianCluster;        
		end
		%%best chosen model stats
        %        bestModelstats = CoefOptimization.bestModelData(validationData.targetsData{iTarget});
        %        bestWeightedMedian(iTarget,1) = bestModelstats.predictionWeightedMedian(1);
        %        bestWeightedMedian(iTarget,2) = bestModelstats.predictionWeightedMedian(2);
        %        
        %        bestMedian(iTarget,1) = bestModelstats.predictionMedian(1);
        %        bestMedian(iTarget,2) = bestModelstats.predictionMedian(2);              
                
        %        bestWeightedMean(iTarget,1) = bestModelstats.predictionWeightedMean(1);
        %        bestWeightedMean(iTarget,2) = bestModelstats.predictionWeightedMean(2);
                 
            end
            stats = {};
	    stats.errorMedian = median(errorMedian);
	    stats.errorMaxWeightedMedian = median(errorMaxWeightedMedian);
	    stats.errorExpWeightedMedian = median(errorExpWeightedMedian);
	    stats.errorWeightedMedian = median(errorWeightedMedian);
	    stats.errorMean = median(errorMean);
	    stats.corrMedian = median(corrMedian);
	    stats.corrWeightedMedian = median(corrWeightedMedian);
	    stats.corrMean = median(corrMean);
                if (validationData.flagComparisonTable~=0)
        		stats.errorClusterWeightedMedian = median(errorClusterWeightedMedian);
		        stats.errorWeightedMedianCluster = median(errorWeightedMedianCluster);
	    		stats.errorRandomWeightedMedian = median(errorRandomWeightedMedian);
			stats.corrRandomWeightedMedian = median(corrRandomWeightedMedian);
        		stats.corrClusterWeightedMedian = median(corrClusterWeightedMedian);
		        stats.corrWeightedMedianCluster = median(corrWeightedMedianCluster);
		end

            disp(['validation by median  and mean ' ...
                  num2str(i) ' ' ...
                  num2str(median(errorMedian)) ' ' ...
                  num2str(median(errorWeightedMedian)) ' ' ...
                  num2str(median(errorMean)) ' ' ...
                  num2str(median(corrMedian)) ' ' ...
                  num2str(median(corrWeightedMedian)) ' ' ...
                  num2str(median(corrMean))]);
        end 

	function stats = printVD(validationData)
            for iTarget = 1:length(validationData.targetsData)
                errorMedian(iTarget)         = validationData.targetsData{iTarget}.statistics.errorMedian;
                errorWeightedMedian(iTarget) = validationData.targetsData{iTarget}.statistics.errorWeightedMedian;
                errorMean(iTarget)           = validationData.targetsData{iTarget}.statistics.errorMean;
                corrMedian(iTarget)          = validationData.targetsData{iTarget}.statistics.corrMedian;
                corrWeightedMedian(iTarget)  = validationData.targetsData{iTarget}.statistics.corrWeightedMedian;
                corrMean(iTarget)            = validationData.targetsData{iTarget}.statistics.corrMean;
                
                
                %%best chosen model stats
                bestModelstats = CoefOptimization.bestModelData(validationData.targetsData{iTarget});
                bestWeightedMedian(iTarget,1) = bestModelstats.predictionWeightedMedian(1);
                bestWeightedMedian(iTarget,2) = bestModelstats.predictionWeightedMedian(2);
                
                bestMedian(iTarget,1) = bestModelstats.predictionMedian(1);
                bestMedian(iTarget,2) = bestModelstats.predictionMedian(2);              
                
                bestWeightedMean(iTarget,1) = bestModelstats.predictionWeightedMean(1);
                bestWeightedMean(iTarget,2) = bestModelstats.predictionWeightedMean(2);
            end
            stats = {};
	    stats.errorMedian = median(errorMedian);
	    stats.errorWeightedMedian = median(errorWeightedMedian);
	    stats.errorMean = median(errorMean);
	    stats.corrMedian = median(corrMedian);
	    stats.corrWeightedMedian = median(corrWeightedMedian);
	    stats.corrMean = median(corrMean);
        stats.bestWeightedMedian = bestWeightedMedian;
        stats.bestMedian = bestMedian;
        stats.bestWeightedMean = bestWeightedMean;
            disp(['validation by median  and mean ' ...
                  num2str(median(errorMedian)) ' ' ...
                  num2str(median(errorWeightedMedian)) ' ' ...
                  num2str(median(errorMean)) ' ' ...
                  num2str(median(corrMedian)) ' ' ...
                  num2str(median(corrWeightedMedian)) ' ' ...
                  num2str(median(corrMean))]);
	
    end
    

        function stats = bestModelData(target)
             %%predictionWeightedMean
                    predictedScorePlace = length(target.fields) - 2;
                    best_model_no = find(target.values(:,1)==max(target.values(:,1)));
                    predicted_model_no = find(target.values(:,predictedScorePlace)==max(target.values(:,predictedScorePlace)));
                    if (length(predicted_model_no)>1)
                            predicted_model_no = predicted_model_no(1);
                    end
                    if (length(best_model_no)>1)
                            best_model_no = best_model_no(1);
                    end
                    stats.predictionWeightedMean(2) = target.values(predicted_model_no,1);
                    stats.predictionWeightedMean(1) = target.values(best_model_no,1);


                    %%predictionMedian
                    predicted_model_no = find(target.values(:,predictedScorePlace+1)==max(target.values(:,predictedScorePlace+1)));
                    if (length(predicted_model_no)>1)
                            predicted_model_no = predicted_model_no(1);
                    end
                    stats.predictionMedian(2) = target.values(predicted_model_no,1);
                    stats.predictionMedian(1) = target.values(best_model_no,1);

                    %%predictionWeightMedian
                    predicted_model_no = find(target.values(:,predictedScorePlace+2)==max(target.values(:,predictedScorePlace+2)));
                    if (length(predicted_model_no)>1)
                            predicted_model_no = predicted_model_no(1);
                    end
                    stats.predictionWeightedMedian(2) = target.values(predicted_model_no,1);
                    stats.predictionWeightedMedian(1) = target.values(best_model_no,1);

        end
       
    end
end






