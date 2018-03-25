classdef ScoreParameters
    properties(Constant=true)
        exponentValues = [1 0.5000 0 -0.5000 -1 -2];
    end
    properties
        fields
        normalizers
        coefs
        exponentIndices
    end
    
    methods
        function obj = ScoreParameters(fields, normalizers, coefs, exponentIndices)
            obj.fields           = fields;
            obj.normalizers      = normalizers;
            obj.coefs            = coefs;
            obj.exponentIndices  = exponentIndices;
        end
        
        function targetDataScore = calctargetDataScore(obj,targetData)
            targetData = targetData.extract(obj.fields,targetData.experimentData);
            
        end
            
    end
    
end

