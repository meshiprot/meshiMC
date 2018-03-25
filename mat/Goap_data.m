
%% This Object is made for the retrivation of data from the GOAP software.
%% Use: Goap_data(`name`); update_Goap(`Goap_Data_path`); update_QMean(`QMean_Data_path`);
%  `name` - is the name of the database being read. for example: 'CASP10'
%  `Goap_Data_path` - path for the Goap data, where all the targets data
%  summaries are in the same folder with *.RESULTS and *.inp extentions.
%  --important notice-- the .inp file must exist for each .RESULTS file.
%  `QMean_Data_path` - is the path for QMean results for each target in the
%  DB. Its one folder, cotaining files of Targets results. each file
%  contains the models data (of a QMean calculation). NOTE: qmean&observed
%  folder.


%%%%%%%
% TODO
%
% 1. Score_data instead of Goap_data, remove the dependencies in actions (Goap and then QMean)
% 2. update_qmean - do not destroy all target if only one model is missing.
% 3. add the possibility for activating update_qmean first.
% 4. go through createExperimentData, look for bugs: missing files, missing models ... 
%
%%%%%%%

classdef Goap_data < handle
    
   properties
   
        name %Name of the data set
        targetsData 
        copyrightAndReadme % Free format description of the data and terms of use. May be printed using the printReadme method.
        statistics
   
   end
   
   methods
        %%
        %constructor
        function obj = Goap_data(name)
            obj.name = name;
        end
        
        function obj = update_Goap(obj,path)
            
           
            results_files = dir([path filesep '*.RESULTS']);
            inp_files = dir([path filesep '*.inp']);
            
            isub = [results_files(:).isdir]; %# returns logical vector
            results_files = {results_files(~isub).name};
            results_files = sort(results_files);
            
            isub = [inp_files(:).isdir]; %# returns logical vector
            inp_files = {inp_files(~isub).name};
            inp_files = sort(inp_files);
            
            for iTarget = 1:length(results_files)
                n1 = regexprep(results_files{iTarget}, '.RESULTS', '');
                n2 = regexprep(inp_files{iTarget}, '.inp', '');
                if (strcmp(n1, n2) == 0)
                    disp(['problem finding results of target ' n1]);
                    continue;
                end
                
                raw_target_data = importdata([path filesep results_files{iTarget}]);
                raw_target_name = importdata([path filesep inp_files{iTarget}]);
                %arr_str = arrayfun(@num2str, raw_target_data.data, 'UniformOutput', false);
                
                if (isempty(raw_target_data))
                    
                    disp(['problem with results of target ' n1 ' file is empty.']);
                    continue;
                end
                
                
                target_data = {};
                target_data.name = n1;
                
                
                target_data.fields{1} = 'goap';
                target_data.fields{2} = 'dfire';
                target_data.fields{3} = 'goap_ag';
                
                target_data.values = raw_target_data.data;
                target_data.values = target_data.values * (-1);
                
                target_data.fileNames = regexprep(raw_target_name(2:end),'pdb.*.*','pdb');
                
                if (length(target_data.fileNames) ~= length(target_data.values))
                    disp([results_files{iTarget} ' file, doesn`t match ' inp_files{iTarget} ' file.']);
                    continue;
                end
                
                
                obj.addTarget(target_data);
            end
            
            
        end
        
        
        function addTarget(obj,target)
            
            [t index] = obj.find_target(target.name);
            if (index < 0)
                obj.targetsData{length(obj.targetsData)+1} = target;
            end
            
        end
        
        %% Creates a new Goap_data obj with all targets updated with the
        % QMean score values. 
        % new_data - the new Goap_data
        % errorList - contains a List of target names which failed
        % combaining the data in obj and the new QMean Data. NOTE: These Targets
        % would not appear in new_data.
        function [new_data errorList] = update_QMean(obj,path)
            
            errorList = {};
            new_data = obj.duplicate(obj.name);
            
            results_files = dir([ path filesep '*.txt' ]);
            isub = [results_files(:).isdir]; %# returns logical vector
            results_files = {results_files(~isub).name};
            results_files = sort(results_files);
           % results_files = regexprep(results_files,'.txt','');
            for iFile = 1:length(results_files)
                
                target_name = regexprep(results_files{iFile}, '.txt', '');

                %get QMean results data                
                raw_target_QMdata = importdata([path filesep results_files{iFile}]);
                
                %%TESTing the file of QMean Data is not corrupted.
                if (isempty(raw_target_QMdata))
                    disp(['Corrupt Target file in QMean Data Target=' target_name])
                    errorList = [errorList target_name];
                    continue;
                end
                
                qmd_TargetNames = regexprep(raw_target_QMdata.textdata,[target_name '.[R|T][0-9]{4}-*[D]*[0-9]*.'],[target_name '.']);

                %get the current target.
                [target iTarget] = new_data.find_target(target_name);
                
                %%TESTing that number of models in obj.target and QMean
                %%Data target is the same.
                if (length(target.fileNames) ~= length(qmd_TargetNames))
                    disp(['not enough data for target ' target_name]);
                    errorList = [errorList target_name];
                    continue;
                end
                
                indices_in_QMdata = cellfun(@(x) find(strcmp(qmd_TargetNames,x)) ,target.fileNames,'UniformOutput', false);

                missing_QMdata = find(cellfun(@isempty,indices_in_QMdata));

                
                %%TESTing for Missing Model Data
                if (~isempty(missing_QMdata)) %% made for verifing that there every model in obj.target has a QMean score value.
                    disp(['Missing Model in QMean Data Target=' target_name])
                    errorList = [errorList target_name];
                    continue;
                end
                
                values = cell2mat(cellfun(@(iQM) raw_target_QMdata.data(iQM,:),indices_in_QMdata,'UniformOutput', false));
                fields = {'QMean_observed', 'QMean_score'};
                
                
                new_data.targetsData{iTarget}.fields = [new_data.targetsData{iTarget}.fields fields];
                new_data.targetsData{iTarget}.values = [new_data.targetsData{iTarget}.values values];   
                
                
            end
            errorList = [errorList setdiff(obj.getTargetNames(),regexprep(results_files, '.txt', ''))];
            new_data.removeTargetList(errorList);
            
        end
        
        function removeTargetList(obj,errorList)
            
            for iTarget = 1:length(obj.targetsData)
                if (~isempty(find(strcmp(errorList,obj.targetsData{iTarget}.name))))
                    obj.targetsData{iTarget} = [];
                end 
            end
            obj.targetsData = obj.targetsData(~cellfun(@isempty,obj.targetsData));
        end
        
        
        %%Duplicate this object
        function dup = duplicate(obj,newName)
            dup = Goap_data(newName);
            dup.copyrightAndReadme = obj.copyrightAndReadme;
            for i = length(obj.targetsData()):-1:1
                dup.targetsData{i}.name = obj.targetsData{i}.name;
                dup.targetsData{i}.values = obj.targetsData{i}.values;
                dup.targetsData{i}.fields = obj.targetsData{i}.fields;
                dup.targetsData{i}.fileNames = obj.targetsData{i}.fileNames;
            end
        end
        
        
        
        %% search for target `target_name` in obj(Goap_data). if it finds it return the target and its index in the database.
        % else return -1 for index and {} for the target.
        function [target index] = find_target(obj,target_name)
            
            target = {};
            index = -1;
            
            for iTarget = 1:length(obj.targetsData)
                if (strcmp(obj.targetsData{iTarget}.name ,target_name) == 1)
                    target = obj.targetsData{iTarget};
                    index = iTarget;
                end    
            end
            
            
        end
        
        %% search for model `model_name` in target `target_name` in obj(Goap_data). if it finds it return the target and its index in the database.
        % else return -1 for index and {} for the target.
        function [model_data index] = find_model(obj,target_name,model_name)
            model_data = {};
            index = -1;
            
            [target_data iTarget] = obj.find_target(target_name);
            
            index = find(strcmp(obj.targetsData{iTarget}.fileNames, model_name));
            model_data =  obj.targetsData{iTarget}.values(index,:);
            
            
        end
   

        %%this function receives experimentData, and create a new
        %%experimentData containing only stats of gdt_ts and goap.
        function experimentData = createEx_addedGoapData(obj,experimentDataBase)
            
            experimentData = ExperimentData(obj.name);
            
             for iTarget = length(experimentDataBase.targetsData()):-1:1
                 
                 %get goap scores for models from obj
                 [goap_target goap_index] = obj.find_target(experimentDataBase.targetsData{iTarget}.targetName);
                 
                 exDBModelNames = regexprep(experimentDataBase.targetsData{iTarget}.fileNames,'.out.[1-9]|.xml','');
                 if (goap_index > 0) 
                     
                     indices_in_goap_EX = cellfun(@(x) find(strcmp(goap_target.fileNames,x)) ,exDBModelNames,'UniformOutput', false);
                     %indices_in_goap_EX = cell2mat(indices_in_goap_EX);
                     values = cellfun(@(i) goap_target.values(i,:),indices_in_goap_EX,'UniformOutput', false);
                     values = cell2mat(values);
                     
                     gdt_ts_index_in_exDB = Util.getIndices(experimentDataBase.targetsData{iTarget}.fields,{'gdt_ts'});                        
                     experimentDataBase.targetsData{iTarget}.values(:,gdt_ts_index_in_exDB);
                     
                     fields = [ experimentDataBase.targetsData{iTarget}.fields goap_target.fields ];
                     values = [ experimentDataBase.targetsData{iTarget}.values values ];
                     experimentData.targetsData{iTarget} = TargetData(experimentDataBase.targetsData{iTarget}.targetName,experimentDataBase.targetsData{iTarget}.fileNames,experimentData,fields,values);
                     
                     
                 else 
                    disp(['Target ' experimentDataBase.targetsData{iTarget}.targetName ' wasn`t found on this Goat_data obj'])
                 %gdt_ts_index = Util.getIndices(obj.fields,{'gdt_ts'});
                 end
                 
                 
             end
            
             experimentData.targetsData = experimentData.targetsData(~cellfun(@isempty,experimentData.targetsData));
              
        end
        %%this function receives experimentData, and create a new
        %%experimentData containing only stats of gdt_ts and goap.
        function experimentData = createExperimentData(obj,experimentDataBase)
            
            experimentData = ExperimentData(obj.name);
            
             for iTarget = length(experimentDataBase.targetsData()):-1:1
                 
                 %get goap scores for models from obj
                 [goap_target goap_index] = obj.find_target(experimentDataBase.targetsData{iTarget}.targetName);
                 
                 exDBModelNames = regexprep(experimentDataBase.targetsData{iTarget}.fileNames,'.out.[1-9]|.xml','');
                 if (goap_index > 0) 
                     
                     indices_in_goap_EX = cellfun(@(x) find(strcmp(goap_target.fileNames,x)) ,exDBModelNames,'UniformOutput', false);
                     %indices_in_goap_EX = cell2mat(indices_in_goap_EX);
                     values = cellfun(@(i) goap_target.values(i,:),indices_in_goap_EX,'UniformOutput', false);
                     values = cell2mat(values);
                     
                     gdt_ts_index_in_exDB = Util.getIndices(experimentDataBase.targetsData{iTarget}.fields,{'gdt_ts'});                        
                     experimentDataBase.targetsData{iTarget}.values(:,gdt_ts_index_in_exDB);
                     
                     fields = [ 'gdt_ts' goap_target.fields ];
                     values = [ experimentDataBase.targetsData{iTarget}.values(:,gdt_ts_index_in_exDB) values ];
                     experimentData.targetsData{iTarget} = TargetData(goap_target.name,exDBModelNames,experimentData,fields,values);
                     
                     
                 else 
                    disp(['Target ' experimentDataBase.targetsData{iTarget}.targetName ' wasn`t found on this Goat_data obj'])
                 %gdt_ts_index = Util.getIndices(obj.fields,{'gdt_ts'});
                 end
                 
                 
             end
            
             experimentData.targetsData = experimentData.targetsData(~cellfun(@isempty,experimentData.targetsData));
              
        end
        
        %%this function receives experimentData, and create a new
        %%experimentData containing only stats of gdt_ts and goap.
        function experimentData = createExperimentData_EnergyList(obj,experimentDataBase, energy_list)
            
            experimentData = ExperimentData(obj.name);
            
             for iTarget = length(experimentDataBase.targetsData()):-1:1
                 
                 %get goap scores for models from obj
                 [goap_target goap_index] = obj.find_target(experimentDataBase.targetsData{iTarget}.targetName);
                 
                 exDBModelNames = regexprep(experimentDataBase.targetsData{iTarget}.fileNames,'.out.[1-9]|.xml','');
                 if (goap_index > 0) 
                     
                     indices_in_goap_EX = cellfun(@(x) find(strcmp(goap_target.fileNames,x)) ,exDBModelNames,'UniformOutput', false);
                     
                     
                     %indices_in_goap_EX = cell2mat(indices_in_goap_EX);
                     values = cellfun(@(i) goap_target.values(i,:),indices_in_goap_EX,'UniformOutput', false);
                     values = cell2mat(values);
                     
                     fields = goap_target.fields;
                     for iEnergy = 1:length(energy_list)
                         iEnergy_in_exDB = -1;
                         iEnergy_in_exDB = Util.getIndices(experimentDataBase.targetsData{iTarget}.fields,energy_list(iEnergy)); 
                         if (iEnergy_in_exDB > 0)
                            experimentDataBase.targetsData{iTarget}.values(:,iEnergy_in_exDB);
                            fields = [ experimentDataBase.targetsData{iTarget}.fields(iEnergy_in_exDB) fields ];
                            values = [ experimentDataBase.targetsData{iTarget}.values(:,iEnergy_in_exDB) values ];
                         end
                     end
                     
                     
                     experimentData.targetsData{iTarget} = TargetData(goap_target.name,exDBModelNames,experimentData,fields,values);
                     
                 else 
                    disp(['Target ' experimentDataBase.targetsData{iTarget}.targetName ' wasn`t found on this Goat_data obj'])
                 %gdt_ts_index = Util.getIndices(obj.fields,{'gdt_ts'});
                 end
                 
                 
             end
             if (~isempty(experimentData.targetsData))
                experimentData.targetsData = experimentData.targetsData(~cellfun(@isempty,experimentData.targetsData));
             end
        end
        
        
        %%this function receives validationData(experimentData), and create a new
        %%experimentData containing only stats of gdt_ts ,goap, and predictions of gdt_ts from CoefOptimization.
        function experimentData = createValidationData(obj,validationDataBase)
            
            %%should not be - 'PlusPredictions' shouldn`t be added to the
            %%Target name
            for iTarget = 1:length(validationDataBase.targetsData)
                validationDataBase.targetsData{iTarget}.targetName = regexprep(validationDataBase.targetsData{iTarget}.targetName,'PlusPredictions','');
            end
            disp (validationDataBase.getTargetNames());
            field_list = { 'predictionWeightedMean'   ; 'predictionMedian'   ; 'predictionWeightedMedian' ;'gdt_ts'; 'delta_gdt_ts';'weightedMdianIDR';'weightedMdianENT'};
            
            % Remove all fields that the validation data doesn't have.
            test = setdiff(field_list,validationDataBase.getFields());
            if ~isempty(test)
                field_list = setdiff(field_list,test);
            end
            
            
            experimentData = obj.createExperimentData_EnergyList(validationDataBase, field_list) ;
        end
        
        
        function target_names = getTargetNames(obj)
            target_names = {};
            %target_names = cell2mat(cellfun(@(iTarget) obj.targetsData{iTarget}.name, obj ));
            for iTarget = 1:length(obj.targetsData)
                target_names = [target_names obj.targetsData{iTarget}.name];
            end

        end
        
   end
   
   
   
   
   methods(Static=true)
   end


end