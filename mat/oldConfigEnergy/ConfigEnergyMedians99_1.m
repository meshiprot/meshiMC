assdef ConfigEnergyMedians99_1 < ConfigEnergy
    properties
    end
    
    methods
        function [energy errors correlations] = evaluate(obj,configuration, coefOptimization) %#ok<MANU>
            [errors correlations] = ConfigEnergy.evaluateTargets(configuration,coefOptimization);
            energy = 0.99*median(errors)+0.01*(1-median(correlations));
        end

    end
    
end

