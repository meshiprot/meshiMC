classdef ExperimentAnalysis
    methods
    end
    
    methods(Static=true)
        function boxPlotStatisticsFields(experimentData,fields,axes)
            for i = length(fields):-1:1
                values(:,i) = getStatisticsField(experimentData,fields{i});
            end
            boxplot(axes,values,'notch','on','labels',fields);
        end
        
       
        
        function boxPlotExperiements(experiments,field,labels,axes)
            values = [];
            groups = [];
            for i = 1:length(experiments)
                experiment = experiments{i};
                newValues = getStatisticsField(experiment,field);
                values = [values; newValues];
                size(groups)
                size(meshgrid(labels{i},(1:length(newValues))'));
                groups = [groups; meshgrid(labels{i},(1:length(newValues))')];
            end
            boxplot(axes,values,groups,'notch','on');
        end
        
        function boxPlotsAllExperiments(origFlag)
            errorH = subplot(1,3,1);
            corrH  = subplot(1,3,2);
            enricH = subplot(1,3,3);
            experimentNames =   {'original  '  'linear    ' 'unweighted'    'FM+TBM    ' 'FM & TBM  '};
            experimentFolders = {'all_start' 'linear' 'all_noWeights' 'all'    'TBM_FM'};
            for i = length(experimentFolders):-1:1
                cd(experimentFolders{i});
                load validationData.mat
                experiments{i} = validationData;
                cd('..');
            end
            ExperimentAnalysis.boxPlotExperiements(experiments, 'errorWeightedMedian',experimentNames,errorH);
            ExperimentAnalysis.boxPlotExperiements(experiments, 'corrWeightedMedian',experimentNames,corrH);
            ExperimentAnalysis.boxPlotExperiements(experiments, 'enrichmentWeightedMedian',experimentNames,enricH);
        end
        
%         function plotTarget(experimentData, fields, iTarget, axis,format)
%             extracted    = experimentData.extract(fields,'temp');
%             target       = extracted.targetsData{iTarget};
%             values       = extr
    end
end


