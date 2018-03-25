classdef RawData < handle
    %RAWDATAINIT 
    %   
    %
    % use - qmean_data = RawData.getRawQMeanData('QMean_data\Casp10Results_tatiana\qmean&observed')
    % goap_data = RawData.getRawGoapData('Goap_data/CASP10')
    % rawDatas = {goap_data qmean_data }
    % union_data = RawData.union(rawDatas,'Hello')
    % new_experimentData = union_data.createEx_addRawData(experimentData)
    properties
        name %Name of the data set
        targetsData 
        targetsErrorList
        copyrightAndReadme % Free format description of the data and terms of use. May be printed using the printReadme method.
        statistics
        
    end
    properties(Constant=true)
        default_CASP_fields = { 'pred' 'gdt_ts'};
	default_Rosetta_fields = {'rosetta_score','rosetta_fa_atr','rosetta_fa_rep','rosetta_fa_sol','rosetta_fa_intra_rep','rosetta_fa_elec','rosetta_pro_close','rosetta_hbond_sr_bb','rosetta_hbond_lr_bb','rosetta_hbond_bb_sc','rosetta_hbond_sc','rosetta_dslf_fa13','rosetta_rama','rosetta_omega','rosetta_fa_dun','rosetta_p_aa_pp','rosetta_ref'};
	default_Voromqa_fields = {'voromqa_score'};
        default_Goap_fields = { 'goap' 'dfire' 'goap_ag'};
        default_QMean_fields = {'QMean_observed', 'QMean_score'};
	default_ProQ2_fields = {'proQ2'}
        default_ComparisonTable_fields = {'rms','gdt_ts'};
    end
    
    methods
        function obj = RawData(name)
            obj.name = name;
            targetsErrorList = {};
        end
        
        %%Duplicate this object
        function dup = duplicate(obj,newName)
            dup = RawData(newName);
            dup.copyrightAndReadme = obj.copyrightAndReadme;
            for i = length(obj.targetsData()):-1:1
                dup.targetsData{i}.name = obj.targetsData{i}.name;
                dup.targetsData{i}.values = obj.targetsData{i}.values;
                dup.targetsData{i}.fields = obj.targetsData{i}.fields;
                dup.targetsData{i}.fileNames = obj.targetsData{i}.fileNames;
            end
        end
        
	function experimentData = createExperimentData(obj)
%		for i=1:length(obj.targetsData)
%		  assert(length(target_data.fields)==size(target_data.values,2));
%		  assert(length(target_data.fileNames)==size(target_data.values,1));
%		end
		experimentData = ExperimentData(obj.name);
		for iTarget = 1:length(obj.targetsData)
			experimentData.targetsData{iTarget} = TargetData(obj.targetsData{iTarget}.name,obj.targetsData{iTarget}.fileNames,[], obj.targetsData{iTarget}.fields, obj.targetsData{iTarget}.values,{},0);
		end
	end        
        
        % createEx_addedRawData
        % This function creates an experiment data based on
        % experimentDataBase with added features in the RawData object.
        %
        %
        function experimentData = createEx_addRawData(obj,experimentDataBase,allowDataLeak)
            if (nargin < 3)
		allowDataLeak = 1;
	    end 
            experimentData = ExperimentData(obj.name);
            
             for iTarget = length(experimentDataBase.targetsData()):-1:1
% 		for iTarget = 65:65                
                 %get goap scores for models from obj
                 [target index] = obj.find_target(experimentDataBase.targetsData{iTarget}.targetName);
                 exDBModelNames = regexprep(experimentDataBase.targetsData{iTarget}.fileNames,'.scwrl.*|.out.[0-9]+|.xml','.');
		 exDBModelNames = regexprep(exDBModelNames,'.*/','.');
		 indices_in_exBase = cellfun(@(x) find(~cellfun(@isempty,strfind(experimentDataBase.targetsData{iTarget}.fileNames,x))) ,exDBModelNames,'UniformOutput', false);


                 if (index > 0) 
                     
                     indices_in_rawData = cellfun(@(x) find(~cellfun(@isempty,strfind(target.fileNames,x))) ,exDBModelNames,'UniformOutput', false);
		     indices_in_exBase	= cell2mat(indices_in_exBase(~cellfun(@isempty,indices_in_rawData)));
                     %indices_in_goap_EX = cell2mat(indices_in_goap_EX);
                     values = cellfun(@(i) target.values(i,:),indices_in_rawData,'UniformOutput', false);
                     values = cell2mat(values);
                     
                     %gdt_ts_index_in_exDB = Util.getIndices(experimentDataBase.targetsData{iTarget}.fields,{'gdt_ts'});                        
                     %experimentDataBase.targetsData{iTarget}.values(:,gdt_ts_index_in_exDB);
                     
                     fields = [ experimentDataBase.targetsData{iTarget}.fields target.fields ];
 		     %%/TODO
%		     disp(iTarget)
% 		     disp(target.name)
%                     disp(size(values))
% 		     disp(size(experimentDataBase.targetsData{iTarget}.values(indices_in_exBase,:)))
                     values = [ experimentDataBase.targetsData{iTarget}.values(indices_in_exBase,:) values ];
                     experimentData.targetsData{iTarget} = TargetData(experimentDataBase.targetsData{iTarget}.targetName,experimentDataBase.targetsData{iTarget}.fileNames(indices_in_exBase),experimentData,fields,values);
                     
                     assert(size(experimentData.targetsData{iTarget}.values,2)==length(experimentDataBase.targetsData{iTarget}.fields)+length(target.fields), ['Error - RawData.createEx_addRawData - Values of target ' target.name  ' was not added correctly - fields not compatible to values']);
		     assert(length(experimentData.targetsData{iTarget}.fileNames)==size(experimentData.targetsData{iTarget}.values,1),['Error - RawData.createEx_addRawData - Values of target ' target.name  ' - values are not compatible to fileNames']);
		     if (nargin > 2)
			if (~allowDataLeak)
			       assert(length(experimentData.targetsData{iTarget}.fileNames)==length(experimentDataBase.targetsData{iTarget}.fileNames))
			end
		     end
                 else 
		   if (~allowDataLeak)
		       assert(index > 0,['Target ' experimentDataBase.targetsData{iTarget}.targetName ' wasn`t found on this RawData obj'])
  		   else
		       disp(['Target ' experimentDataBase.targetsData{iTarget}.targetName ' wasn`t found on this RawData obj']);
		   end	
                 %gdt_ts_index = Util.getIndices(obj.fields,{'gdt_ts'});
                 end
             end
            
             experimentData.targetsData = experimentData.targetsData(~cellfun(@isempty,experimentData.targetsData));
              
        end
        
        %%this function receives validationData(experimentData), and create a new
        %%experimentData containing only stats of gdt_ts ,goap, and predictions of gdt_ts from CoefOptimization.
        function [errors experimentData] = createValidationData(obj,validationDataBase)
            
            %%should not be - 'PlusPredictions' shouldn`t be added to the
            %%Target name
            for iTarget = 1:length(validationDataBase.targetsData)
                validationDataBase.targetsData{iTarget}.targetName = regexprep(validationDataBase.targetsData{iTarget}.targetName,'PlusPredictions','');
            end
            disp (validationDataBase.getTargetNames());
            field_list = { 'predictionWeightedMean'   ; 'predictionMedian'   ; 'predictionWeightedMedian' ;'gdt_ts'; 'delta_gdt_ts';'weightedMdianIDR';'weightedMdianENT'; 'goap'};
            
            % Remove all fields that the validation data doesn't have.
            test = setdiff(field_list,validationDataBase.getFields());
            if ~isempty(test)
                field_list = setdiff(field_list,test);
            end
            
            
            [errors experimentData] = obj.createExperimentData_EnergyList(validationDataBase, field_list) ;
        end
        
        
        % createExperimentData_EnergyList
        % This function creates a new experimentData with energy_list
        % from experimentDataBase + fields from this RawData object. 
        % 
        % Notice: targets that do not appear in this RawData object are
        % ignored and won't appear in the new experimentData.
        function [errors experimentData] = createExperimentData_EnergyList(obj,experimentDataBase, energy_list)
            
            experimentData = ExperimentData(obj.name);
            
             for iTarget = length(experimentDataBase.targetsData()):-1:1
                 
                 %get goap scores for models from obj
                 [goap_target goap_index] = obj.find_target(experimentDataBase.targetsData{iTarget}.targetName);
                 
                 exDBModelNames = regexprep(experimentDataBase.targetsData{iTarget}.fileNames,'.scwrl.*.xml','');
                 if (goap_index > 0) 
                     
                     indices_in_goap_EX = cellfun(@(x) find(strcmp(goap_target.fileNames,x)) ,exDBModelNames,'UniformOutput', false);
%                     goap_target.fileNames(indices_in_goap_EX); %names in the rawData 
                     dontExistFileNamesInRawData = setdiff(exDBModelNames,goap_target.fileNames(cell2mat(indices_in_goap_EX))); % file names that dont exist in rawData
		     e.targetName = experimentDataBase.targetsData{iTarget}.targetName;
		     e.fileNames = dontExistFileNamesInRawData;

		     foundExDBModelNames = setdiff(exDBModelNames,dontExistFileNamesInRawData); % file names that exists in RawData AND TargetData
                     indices_in_BaseEX = cellfun(@(x) find(strcmp(exDBModelNames,x)) ,foundExDBModelNames,'UniformOutput', false);

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
                            values = [ experimentDataBase.targetsData{iTarget}.values(cell2mat(indices_in_BaseEX),iEnergy_in_exDB) values ];
                         end
                     end
                     
		     if (length(e.fileNames)>0)
			     errors{iTarget} = e;                     
		     end
                     experimentData.targetsData{iTarget} = TargetData(goap_target.name,exDBModelNames,experimentData,fields,values);
                     
                 else 
                    disp(['Target ' experimentDataBase.targetsData{iTarget}.targetName ' wasn`t found on this RawData obj'])
                 %gdt_ts_index = Util.getIndices(obj.fields,{'gdt_ts'});
                 end
                 
                 
             end
             if (~isempty(experimentData.targetsData))
                experimentData.targetsData = experimentData.targetsData(~cellfun(@isempty,experimentData.targetsData));
             end
	     errors = errors(~cellfun(@isempty,errors))
        end
        
        function removeTargetList(obj,errorList)
            
            for iTarget = 1:length(obj.targetsData)
                if (~isempty(find(strcmp(errorList,obj.targetsData{iTarget}.name))))
                    obj.targetsData{iTarget} = [];
                end 
            end
            obj.targetsData = obj.targetsData(~cellfun(@isempty,obj.targetsData));
        end
        
        function addTarget(obj,target)
            
            [t index] = obj.find_target(target.name);
            if (index < 0)
                obj.targetsData{length(obj.targetsData)+1} = target;
            end
            
        end
        function addTarget2ErrorList(obj,targetName)
            if (isempty(obj.targetsErrorList))
                obj.targetsErrorList = {targetName};
                return
            end
            index = getnameidx(obj.targetsErrorList,{targetName});
            if (index == 0)
                obj.targetsErrorList = [obj.targetsErrorList targetName];
            end
        end
        
        % search for target `target_name` in obj. if it finds it return the target and its index in the database.
        % else return -1 for index and {} for the target.
        function [target index] = find_target(obj,target_name)
            
            target = {};
            index = -1;
            
            for iTarget = 1:length(obj.targetsData)
                if (strfind(obj.targetsData{iTarget}.name ,target_name) == 1)
                    target = obj.targetsData{iTarget};
                    index = iTarget;
                end    
            end
            
            
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
                
        
        
        
        % getRawData
        % This function is used for the collection of raw targets data.
        %
        % path - the path of the target data files. The data is organized
        % in files. each file contains the data of one target.
        %
        % format - is the format of the target files. for example:
        %          [model] [field1] [field2]...[fieldk]
        %
        % pattern -  the data files pattern. for example: for
        % T0645_D1.RESULTS the pattern is '*.RESULTS'
        %
        % fields - the fields names that the data files coresponds with.
        %function getRawData(path, name, format, pattern, fields)
            
        %end
        
        function obj = format1_RawDataCollection(path, obj_name, pattern, fields)
            obj = RawData(obj_name);
            results_files = RawData.getFileList(path,pattern);
            targetNameExtractionAgent1 = regexprep(pattern, '*', '');
            for iTarget = 1:length(results_files)
                target_name = regexprep(results_files{iTarget}, targetNameExtractionAgent1, '');
		target_name = regexprep(target_name,'.*/','');
                raw_target_data = importdata(results_files{iTarget});
                %raw_target_data = importdata([path filesep results_files{iTarget}]);
                %%TESTing the file of QMean Data is not corrupted.
                if (isempty(raw_target_data))
                    disp(['Corrupt Target file in QMean Data Target=' target_name]);
                    %errorList = [errorList target_name];
                    continue;
                end
                qmd_TargetModelNames = regexprep(raw_target_data.textdata,[target_name '.[R|T][0-9]{4}-*[D]*[0-9]*.'],[target_name '.']);       
                target_data = {};
                target_data.name = target_name;
                target_data.fields = fields;
                target_data.values = raw_target_data.data ;
                target_data.fileNames = qmd_TargetModelNames;
                obj.addTarget(target_data);
		
		assert(length(target_data.fields)==size(target_data.values,2),'Error -  RawData.format1_RawDataCollection - (target_data.fields~=target_data.values)');
		assert(length(target_data.fileNames)==size(target_data.values,1),'Error -  RawData.format1_RawDataCollection - (target_data.fileNames~=target_data.values)');
		
            end
        end
        
        % format2_RawDataCollection
        % This function is for the collection of Model X Model comparison
        % data.
        function obj = format2_RawDataCollection(path, obj_name, pattern, fields)
            obj = RawData(obj_name);
            results_files = RawData.getFileList(path, pattern);
            targetNameExtractionAgent = regexprep(pattern, '*', '');
            for iTarget = 1:length(results_files)
                target_name = regexprep(results_files{iTarget}, targetNameExtractionAgent, '');
                raw_target_data = importdata([path filesep results_files{iTarget}]);
                %%TESTing the file of Raw Data is not corrupted.
                if (isempty(raw_target_data))
                    disp(['Corrupt Target file in QMean Data Target=' target_name]);
                    obj.addTarget2ErrorList(target_name); 
                    continue;
                end
                [s inx] = unique(strcat(raw_target_data.textdata(:,1)));
                modelNamesAndIndex = raw_target_data.textdata(inx,1);
                
                modelsNumber = length(modelNamesAndIndex);
                comparisonNumber = length(raw_target_data.textdata(:,1));
                if (modelsNumber*modelsNumber ~= comparisonNumber) % check to verify all of the models comparisons have successfully finished.
                    disp(['Corrupt Target file in QMean Data Target=' target_name]);
                    obj.addTarget2ErrorList(target_name); 
                    continue;
                end
                comparisonTables = {};
                for iComparison = 1:comparisonNumber
                    indices = Util.getIndices(modelNamesAndIndex,raw_target_data.textdata(iComparison,:));
                    iModel = indices(1);
                    jModel = indices(2);
                    
                    
                    for iField = 1:length(fields)
                        comparisonTables{iField}.CMTable(iModel,jModel) = raw_target_data.data(iComparison,iField);
                    end                    
                    %rmsComparisonTable(iModel,jModel) = raw_target_data.data(iComparison,1);
                    %gdtComparisonTable(iModel,jModel) = raw_target_data.data(iComparison,2);
                    
                end  
                target_data = {};
                target_data.name = target_name;
                target_data.fields = fields;
                target_data.comparison_tables = raw_target_data.data ;
                target_data.fileNames = modelNamesAndIndex;
                obj.addTarget(target_data);
            end
        end
        

        
        function obj = getRawComparisonTableData(path)
            obj = RawData.format2_RawDataCollection(path, 'ComparisonTable', '*_similarityTable.dat', RawData.default_ComparisonTable_fields);
        end
        
        function obj = getRawQMeanData(path)
            obj = RawData.format1_RawDataCollection(path, 'QMean_Data', '*.txt', RawData.default_QMean_fields);  
        end
        function obj = getRawProQ2Data(path)
            obj = RawData.format1_RawDataCollection(path, 'ProQ2_Data', '*.SCORE', RawData.default_ProQ2_fields);  
        end
        function obj = getRawRosettaData(path)
            obj = RawData.format1_RawDataCollection(path, 'Rosetta_Data', '*.RSCORE', RawData.default_Rosetta_fields); 
        end
        function obj = getRawVoromqaData(path)
            obj = RawData.format1_RawDataCollection(path, 'Voromqa_Data', '*.voromqa', RawData.default_Voromqa_fields);  
        end
        function obj = getRawCASPData(path)
            obj = RawData.format1_RawDataCollection(path, 'CASP_Data', '*.SCORES', RawData.default_CASP_fields);  
        end
                
        
        function obj = getRawGoapData(path)
            % create the RawData object
            fields  = RawData.default_Goap_fields;
            name    = 'Goap_Data';
            pattern1 = '*.RESULTS';
            pattern2 = '*.inp';
            obj = RawData(name);
            
            results_files = RawData.getFileList(path,pattern1);
            inp_files = RawData.getFileList(path,pattern2);
            targetNameExtractionAgent1 = regexprep(pattern1, '*', '');
            targetNameExtractionAgent2 = regexprep(pattern2, '*', '');
            for iTarget = 1:length(results_files)
                n1 = regexprep(results_files{iTarget}, targetNameExtractionAgent1, '');
                n2 = regexprep(inp_files{iTarget}, targetNameExtractionAgent2, '');
                if (strcmp(n1, n2) == 0)                                                    % input file compatability test.
                    disp(['problem finding results of target ' n1]);
                    continue;
                end
                raw_target_data = importdata([path filesep results_files{iTarget}]);
                if (isempty(raw_target_data))                                               % data file size test
                    disp(['problem with results of target ' n1 ' file is empty.']);
                    continue;
                end
                target_data = {};
                target_data.name = n1;
                target_data.fields = fields;
                target_data.values = raw_target_data.data * (-1);
                raw_target_name = importdata([path filesep inp_files{iTarget}]);
                target_data.fileNames = regexprep(raw_target_name(2:end),'pdb.*.*','pdb');
                if (length(target_data.fileNames) ~= length(target_data.values))             % input file - results file, number of results compatability test.
                    disp([results_files{iTarget} ' file, doesn`t match ' inp_files{iTarget} ' file.']);
                    continue;
                end 
                obj.addTarget(target_data);
            end
        end
            
        function obj = union(rawDataList,name)
            assert(iscell(rawDataList),'The rawDataList variable should be a cell-array');
            assert(length(rawDataList)>0,'The rawDataList variable should be at length > 0');
            cellfun(@(x) assert(strcmp(class(x),'RawData'),'Each cell in rawDataList has to be of type RawData'),rawDataList);
            obj = rawDataList{1}.duplicate(name);
            for iRawData = 2:length(rawDataList)
                temp = obj.duplicate(name);
                
                % Go over each target in the next RawData.
                for iTarget = 1:length(rawDataList{iRawData}.targetsData)
                    listTargetData = rawDataList{iRawData}.targetsData{iTarget};
                    [current_target iCurrentTarget] = temp.find_target(listTargetData.name);%get the current target.
                    if (iCurrentTarget<=0) %couldn't find target in obj.targetsData
                        disp(['could not find target ' listTargetData.name ' in rawData']);
                        temp.addTarget2ErrorList(listTargetData.name); 
                        continue;
                    end
                    % Testing - verifing the taregets in
                    % rawDataList{iRawData} are not currupted.
                    fieldsIndecesTarget = find(~getnameidx(current_target.fields, listTargetData.fields));
                    if (length(current_target.fileNames)~= length(listTargetData.fileNames))
                        disp(['RawData ' num2str(iRawData) ' has different number of models. Target ' listTargetData.name ' is deleted.' ]);
                        temp.addTarget2ErrorList(listTargetData.name);
                        continue;
                    end
                    modelsIndicesInlistTargetData = getnameidx(listTargetData.fileNames,current_target.fileNames);
                    if (~isempty(find(~modelsIndicesInlistTargetData,1)))
                        % There is different models in the RawDatas.
                        disp(['RawData ' num2str(iRawData) 'has incompatible models lists. Target ' listTargetData.name 'is deleted.' ]);
                        temp.addTarget2ErrorList(listTargetData.name);
                        continue;
                    end
                    current_target.fields = [current_target.fields listTargetData.fields(fieldsIndecesTarget)];
                    values = listTargetData.values(modelsIndicesInlistTargetData,fieldsIndecesTarget);
                    current_target.values = [current_target.values values];
                    temp.targetsData{iCurrentTarget} = current_target;
                end
                % Remove targets that doesn't have data in the next
                % RawData.
                temp.addTarget2ErrorList(setdiff(temp.getTargetNames(),rawDataList{iRawData}.getTargetNames()));
                temp.removeTargetList(temp.targetsErrorList);
                obj = temp;
            end           
        end
        
       
        % errorCheck1
        %
        % check for raw data size variations
        function error = errorCheck1(rawTargetData,targetName)
            if (isempty(rawTargetData))
                disp(['problem with results of target ' targetName ' ,file is empty.']);
                error = 1;
            else
                error = 0;
            end
        end
        
        
        % errorCheck2
        % 
        % format2 assosiation. check for Model X Model completeness. 
        % verify that for each models mi,mj there is a value.
        function errorCheck2()
        end
        
        
        function fileNamesList = getFileList(path, pattern)
 	    %%changed to enable specific files in many directories. for example: DBTest/*/resetta_score/*.RSCORE
            %fileNamesList = dir([path filesep pattern]);
            %isub = [fileNamesList(:).isdir]; %# returns logical vector
            %fileNamesList = {fileNamesList(~isub).name};
            %fileNamesList = sort(fileNamesList);
	    fileNamesList = strsplit(ls([path filesep pattern]));
	    fileNamesList = fileNamesList(~cellfun(@isempty,fileNamesList));
        end  



        
    end
    
end

