classdef ConfigEnergyErrors < ConfigEnergy
    properties
    end
    
    methods
        function [energy errors correlations] = evaluate(obj,configuration, objective,... 
                                                         coefOptimization) %#ok<MANU> %fourth argument just for debugging
            if(nargin > 3) % that is debug
                 [errors correlations] = ConfigEnergy.evaluateTargets(configuration,objective,coefOptimization);
            else
                 [errors correlations] = ConfigEnergy.evaluateTargets(configuration,objective);
            end
            energy = sqrt(mean(errors.*errors));
        end

                  function [energies correlations] = evaluateArray(obj,configurationArray,coefOptimization)
                nTargets     = length(configurationArray);
                energies     = zeros(nTargets,1);
                correlations = zeros(nTargets,1);
                for i = nTargets:-1:1
                    [energy, ~, tempCorrelations] = obj.evaluate(configurationArray{i},coefOptimization);
                    energies(i,1)     = energy;
                    correlations(i,1) = median(tempCorrelations);
                end
         end


    end
    
end
