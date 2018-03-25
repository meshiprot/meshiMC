classdef ConfigEnergyEnrichment < ConfigEnergy
    properties
    end
    
    methods
        function [energy errors correlations] = evaluate(obj,configuration, objective,... 
                                                         coefOptimization) %#ok<MANU> %fourth argument just for debugging
            if(nargin > 3) % that is debug
                 [errors correlations, enrichments] = ConfigEnergy.evaluateTargets(obj,configuration,objective,coefOptimization);
            else
                 [errors correlations, enrichments] = ConfigEnergy.evaluateTargets(obj,configuration,objective);
%errors'
%correlations'
            end
%            energy = 0.10*sqrt(mean(errors.*errors))+0.9*(1-median(correlations));
            energy = 1/(1+sqrt(mean(enrichments.*enrichments)));
        end

    end
    
end
