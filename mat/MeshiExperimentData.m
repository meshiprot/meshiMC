classdef MeshiExperimentData < ExperimentData

    properties
    end
    
    methods
        function obj = MeshiExperimentData(filesString)
               obj.targetsDataType = 'MeshiTargetData';
               if (nargin == 1)
                   obj = MeshiExperimentData.readExperiemntData(filesString);
               end
         end
           
        function scores = calcScore(obj,parametersArray,flag)
           scores = MeshiExperimentData(); %new object
           [score interquartile interdecile] = arrayfun(@(X)MeshiExperimentData.calcScoreTarget(X,parametersArray,flag),obj.targetsData);
                scores.fields = {'score' 'interquartile' 'interdecile'};
                for i = length(score):-1:1
                    values    = cell2mat({score{i} interquartile{i} interdecile{i}});
                    fileNames = obj.targetsData{i}.fileNames;
                    scores.targetsData{i} = TargetData(values,fileNames,obj);
                end
           end     
    end
    
    methods(Static = true)
            function [score interquartile interdecile] = calcScoreTarget(target, parametersArray,flag)
                    target = target{1};
                    disp(target.fileNames(1));
                    [score interquartile interdecile] = target.calcScore(parametersArray);
                    if (flag ~= 1)
                        badL = score < 0.11;
                        badH = score > 1;
                        if (flag == 2)
                            score(badL) = 0.11;
                            score(badH) = 1;
                        else
                            bad  = badL | badH;
                            score(bad) = -9999;
                        end
                    end
            end
            
            function obj = readExperiemntData(filesString)
                    files = dir(filesString);
                    obj = MeshiExperimentData();
                    fileNames = {files.name};
                    n = length(fileNames);
                    obj.fields = '';
                    emptyTargets = false(n,1);
                    for i = n:-1:1 %Starting from the last element is a way of pre allocating memory
                        obj.targetsData{i} = TargetData([],[],[]);
                        fileName = fileNames{i};
                        disp(fileName);
                        targetData = load(fileName);
                        targetData = targetData.targetData;
                        if (~isempty(targetData.values))
                            if (isempty(obj.fields))
                                obj.fields = targetData.fields;
                            else
                                if ((length(obj.fields) ~= length(targetData.fields))||...
                                    (length(find(strcmp(obj.fields,targetData.fields))) < length(obj.fields)))
                                    error('This is weird')
                                end
                            end
                            obj.targetsData{i} = TargetData(targetData.values, targetData.fileNames,obj);
                        else 
                            emptyTargets(i) = true;
                        end
                    end
                    obj.targetsData = {obj.targetsData{~emptyTargets}};
               end
    end
end

