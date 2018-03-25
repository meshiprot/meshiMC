    function parametersArray = fixPA(parametersArray)
    parametersArray = arrayfun(@(X)fix(X),parametersArray);
end

function parameters = fix(parameters)
    parameters.header = {'length' parameters.header{:}};
%    parameters.normalizers = {'length' parameters.normalizers{:}};
    parameters.coefs = [0; parameters.coefs];
    parameters.exponentIndices =  [1 1;parameters.exponentIndices];
end
    