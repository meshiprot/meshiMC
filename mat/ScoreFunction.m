classdef ScoreFunction < handle
    properties
        exponentValues
        configuration
        exponentWeightMatrix
        values;
        weightedValues
        nModels
    end
    
    methods
        function obj = ScoreFunction(configuration,targetData,exponentValues)
            temp = targetData.extract(configuration.fields,[],1);
            obj.values = temp.values;
            obj.nModels = size(obj.values,1);
            obj.configuration = configuration;
            obj.configuration.changedFields = 1:obj.configuration.nFields;
            obj.exponentValues = exponentValues;
        end
        
        
        function targetDataScores = calcScoreTarget(obj)
                obj.getExponentWeightsMatrix();
                obj.getWeightsTarget();
                targetDataScores = obj.weightedValues(:,obj.configuration.fieldIndices)*...
                                   obj.configuration.coefs(obj.configuration.fieldIndices,1);
        end
        
        function getExponentWeightsMatrix(obj)
                obj.exponentWeightMatrix = zeros(obj.nModels,...
                                       length(obj.exponentValues)+1,... %The last column is only zeros
                                       length(obj.configuration.normalizerIndices)); %The first layer is length based weights
                normalizerIndices = obj.configuration.normalizerIndices;
                for i = 1:obj.nModels
                    for j = 1:length(obj.exponentValues)
                        for k = 1:length(obj.configuration.normalizerIndices)
                            obj.exponentWeightMatrix(i,j,k)= getWeight(obj.exponentValues(j),...
                                                                       obj.values(i,normalizerIndices(k)));
                        end
                    end
                end
                if (sum(sum(sum(isnan(obj.exponentWeightMatrix))))>0)
                    error('weird exponentWeightMatrix')
                end
                
                function weight = getWeight(exponent,l)
                    
                    if ((l <= 0) || (isnan(l)))
                        weight = 0;
                    else
                        switch(exponent)
                            case 1
                                weight = l;
                            case 0.5
                                weight = sqrt(l);
                            case 0
                                weight = 1;
                            case -0.5
                                weight = 1.0/sqrt(l);
                            case -1
                                weight = 1.0/l;
                            case -2
                                weight = 1.0/(l*l);
                            case -3
                                weight = 1.0/(l*l*l);
                            otherwise
                                error(['weird exponent ' num2str(exponent)]);
                        end
                    end
                    if (isnan(weight))
                        error(['weird weight ' num2str(exponent) num2str(l)]);
                    end
                    if (isinf(weight))
                        error(['weird weight 2 ' num2str(exponent) '  ' num2str(l)]);
                    end
                end
        end
        function getWeightsTarget(obj)
                changedFields = obj.configuration.changedFields;
                for i = length(changedFields):-1:1
		    if (changedFields(i) > length(obj.values(1,:)))
                        error('weird changed features');
                    end
                    obj.weightedValues(:,changedFields(i)) = ...
                                           obj.values(:,changedFields(i));
                    for j = 1:length(obj.configuration.exponentIndices(1,:))
                        exI     = obj.configuration.exponentIndices(changedFields(i),j);
                        if (exI <= 0) 
                           % disp('This is weird exI == 0. Please take a look.');
                            exI = 1;
                        end
                        weights = obj.exponentWeightMatrix(:,exI,j);
                        obj.weightedValues(:,changedFields(i)) = obj.weightedValues(:,changedFields(i)).*weights;
                    end
                end
        end

    end

end

