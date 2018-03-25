classdef ConfigEnergyWeightedErrorsAndMedianCorrelation < ConfigEnergy
    properties
    end
    
    methods
        function [energy errors correlations] = evaluate(obj,configuration, objective,... 
                                                         coefOptimization) %#ok<MANU> %fourth argument just for debugging
            if(nargin > 4) % that is debug
                 [errors correlations] = ConfigEnergy.evaluateTargets(obj,configuration,objective,coefOptimization);
            else
                 [errors correlations] = ConfigEnergy.evaluateTargets(obj,configuration,objective);
%errors'
%correlations'
            end
%            energy = 0.10*sqrt(mean(errors.*errors))+0.9*(1-median(correlations));
            energy = 0.99*sqrt(mean(errors.*errors))+0.01*(1-median(correlations));
        end
        
        function [error, correlation] = evaluateTarget(obj,configuration,scoreFunction,coefs,targetObjective,coefOptimization,iTarget)
             if (std(targetObjective.values) == 0)
	          targetObjective.values = targetObjective.values + (rand(size(targetObjective.values))-0.5)*0.001;
                  targetObjective.values(targetObjective.values>1)=1;
             end	
             predictions = scoreFunction.weightedValues(:,configuration.fieldIndices)*coefs;
            if (std(predictions) == 0)
                  predictions = predictions + (rand(size(predictions))-0.5)*0.001;
                  predictions(predictions>1)=1;
            end
             
            wScore = configuration.sigmoid.backSigmoid(targetObjective.values);
            wScore = configuration.sigmoid.derivative(wScore);%the weight (derivative of sigmoid) of the objective function 
            wPredScore = configuration.sigmoid.backSigmoid(predictions);
            wPredScore = configuration.sigmoid.derivative(wPredScore);%the weight (derivative of sigmoid) of the objective function 
            weight =  ( (wScore + wPredScore)/sum(wScore + wPredScore) );
            
            diff        = (predictions-targetObjective.values);
            correlation = corr([predictions targetObjective.values]);
            correlation = correlation(1,2);
            error       = sqrt(mean(weight.*diff.*diff));
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
