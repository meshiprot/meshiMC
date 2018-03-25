classdef ConfigEnergyMedianCorrelation < ConfigEnergy
    properties
    end
    
    methods
        function [energy errors correlations] = evaluate(obj,configuration, objective,... 
                                                         coefOptimization) %#ok<MANU> %fourth argument just for debugging
            if(nargin > 4) % that is debug
                 [errors correlations] = ConfigEnergy.evaluateTargets(obj,configuration,objective,coefOptimization);
            else
                 [errors correlations] = ConfigEnergy.evaluateTargets(obj,configuration,objective);
            end
            energy = 0.50*sqrt(mean(errors.*errors))+0.5*(1-median(correlations));
        end

    end
    
end
