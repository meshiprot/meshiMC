classdef ConfigEnergy
    properties
    end
    
    methods(Static=true)         
        
         function [errors, correlations, enrichments] = evaluateTargets(obj,configuration,objective,...
                                          coefOptimization) %the third argument is only for debugging
            errors       = zeros(length(objective.targetsData),1);
            correlations = zeros(length(objective.targetsData),1);
	    enrichments  = zeros(length(objective.targetsData),1);

            fields = configuration.fieldIndices;
            coefs  = configuration.coefs(fields,1);
            for iTarget = 1:length(configuration.learningSet)
                scoreFunction = configuration.scoreFunctions{iTarget};
                targetObjective = objective.targetsData{iTarget};
                if (nargin > 4) % that is debugging
                    [error, correlation, enrichment] = obj.evaluateTarget(configuration,scoreFunction,coefs,targetObjective,coefOptimization,iTarget);
                else
                    [error, correlation, enrichment] = obj.evaluateTarget(configuration,scoreFunction,coefs,targetObjective);
                end
                errors(iTarget) = error;
                correlations(iTarget) = correlation;
		enrichments(iTarget) = enrichment;
            end
         end
    end
    methods  
         function [error, correlation, enrichment] = evaluateTarget(obj,configuration,scoreFunction,coefs,targetObjective,coefOptimization,iTarget)
             if (std(targetObjective.values) == 0)
	          targetObjective.values = targetObjective.values + (rand(size(targetObjective.values))-0.5)*0.001;
                  targetObjective.values(targetObjective.values>1)=1;
             end	
             predictions = scoreFunction.weightedValues(:,configuration.fieldIndices)*coefs;
            if (std(predictions) == 0)
                  predictions = predictions + (rand(size(predictions))-0.5)*0.001;
                  predictions(predictions>1)=1;
             end
              
            diff        = predictions-targetObjective.values;
            correlation = corr([predictions targetObjective.values],'type','spearman');
            correlation = correlation(1,2);
            error       = sqrt(mean(diff.*diff));
	    enrichment  = Util.getEnrichment(predictions, targetObjective.values,0.9);
            if (nargin > 5)
                    %%%%%%%% for dubugging %%%%%%%%%
                e = sqrt(mean(diff.*diff)); %for debvugging;
                temp = coefOptimization.objective.concatenate(coefOptimization.features,'xxxx');
                changedFields = {configuration.fields{configuration.changedFields}} ; %#ok<CCAT1>
                [e2, validateExponentValues1, validationCoefs, validationConfiguration] = configuration.validateEnergy(temp.targetsData{iTarget}); %for debugging
                disp([e e2]);
                if (e ~= e2)
                    disp('non equal energies:');
                    disp([e e2] );
                    temp = configuration.scoreFunctions{iTarget}.weightedValues(:,configuration.fieldIndices);
                    disp([validateExponentValues1 ; validationCoefs';...
                     temp(1,:); coefs(:,1)']);
                     fields1 = {configuration.fields{configuration.fieldIndices}};
                     fields2 = { validationConfiguration.fields{configuration.fieldIndices}};
                     disp('difference in field');
                     disp(fields1{temp(1,:) ~= validateExponentValues1});
                     temp(1,temp(1,:) ~= validateExponentValues1)
                      fields2{temp(1,:) ~= validateExponentValues1}; %#ok<VUNUS>
                      disp(fields2);
                      validateExponentValues1(1,temp(1,:) ~= validateExponentValues1) %#ok<NOPRT>
                      [configuration.exponentIndices(configuration.fieldIndices,:)';...
                     validationConfiguration.exponentIndices(configuration.fieldIndices,:)'] %#ok<NOPRT>
                     disp('changed fields')
                     disp(changedFields)
                    error('energy problem');
                end
            end
         end
    end
    
end

