classdef ConfigEnergyMedians < ConfigEnergy
    properties
    end
    
    methods
        function [energy errors correlations] = evaluate(obj,configuration, coefOptimization) %#ok<MANU>
            [errors correlations] = ConfigEnergy.evaluateTargets(configuration,coefOptimization);
            energy = 0.9*median(errors)+0.1*(1-median(correlations));
        end

    end
    
end
