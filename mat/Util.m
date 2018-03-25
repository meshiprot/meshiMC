classdef Util
    properties
    end
    
    methods(Static=true)
        function indices = getIndices(list,keys)
            listSize = length(keys);
            indices = ones(1,listSize);
            for iKey=1:listSize
                key = keys{iKey};
                index = find(strcmp(list,key));
                if (isempty(index))
%                    disp(list);
                    error(['key not found ' key]);
                end
                if (length(index) > 1)
                    error(['Key found more than once ',key]);
                end
                indices(iKey)= index;
            end
        end
        
        function [error, correlation, enrichment] = getErrorAndCorr(prediction, objective)
            diff        = prediction-objective;
            c           = corr([prediction objective]);
            error       = sqrt(mean(diff.*diff));
            correlation = c(1,2);
            enrichment  = Util.getEnrichment(prediction, objective,0.95);
        end
        
%        function enrichment     = getEnrichment(prediction, objective,x)
%            
%            
%            sortedObjective    = sort(objective);          
%            topObjective        = sortedObjective(round(length(objective)*x))
%            
%            topFraction         = length(find(sortedObjective >= topObjective))/length(prediction);
%            
%            sortedPredictions   = sort(prediction);
%            topPrediction       = sortedPredictions(round(length(objective)*x))
%            topIndices          = find(prediction>=topPrediction)
%            
%            find(objective(topIndices) >= topObjective)
%            both                = length(find(objective(topIndices) >= topObjective))/length(topIndices);
%            enrichment          = both/topFraction;
%        end
        function enrichment     = getEnrichmentOld(prediction, objective,x)
            sortedObjective     = sort(objective);
            topObjective        = sortedObjective(round(length(objective)*x));
            
            topFraction         = length(find(sortedObjective >= topObjective))/length(prediction);
            
            sortedPredictions   = sort(prediction);
            topPrediction       = sortedPredictions(round(length(objective)*x));
            topIndices          = find(prediction>=topPrediction);
            
            find(objective(topIndices) >= topObjective);
            both                = length(find(objective(topIndices) >= topObjective))/length(topIndices);
            enrichment          = both/topFraction;
        end

    function enrichment     = getEnrichment(prediction, objective,x)
            
            prediction          = [prediction, (1:length(prediction))'];
            objective           = [objective, (1:length(objective))'];
            
            
            sortedObjective           = sortrows(objective);
            topObjectiveIndices       = sortedObjective(round(length(objective)*x):end,2);
            
            topFraction         = length(topObjectiveIndices)/length(objective); %normalize x to the number of decoys
            
            sortedPredictions    = sortrows(prediction);
            topPredictionIndices = sortedPredictions(round(length(objective)*x):end,2);
            
            intersection = length(intersect(topObjectiveIndices,topPredictionIndices))/length(topObjectiveIndices);

            enrichment          = intersection/topFraction;
        end
                   
	function combinedConfigurationArray = combine(filesToCombine)
        files = dir(filesToCombine);
        list = {files.name};
        combinedConfigurationArray = {};
		for i = 1:length(list)
			load(list{i});
			combinedConfigurationArray = {combinedConfigurationArray{:} configurationArray{:}}; %#ok<USENS> This is the loaded object.
		end
	end
	
	function energies = validateConfigurationArray(experimentData, configurationArray)
	    energies = zeros(length(configurationArray),1);
	    for i = 1:length(configurationArray)
		fprintf(1,'.');
                if (mod(i,50) == 0)
			fprintf('\n');
	        end
                energies(i) = Configuration.validateConfiguration(configurationArray{i},experimentData,debugFlag);
            end
    end

    function weights = getCMTWeights(comparisonTable)
        [length width] = size(comparisonTable);
        assert(length==width,['The Table size invalid, length=' num2str(length) ' width=' num2str(width)]);
        cm1 = sum(comparisonTable);
        cm2 = 1./cm1;
        cm3 = cm2./sum(cm2);
        weights = cm3';
    end    
    
    
    end
end

