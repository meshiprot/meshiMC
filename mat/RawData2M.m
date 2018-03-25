classdef RawData2M < handle
    %RAWDATAINIT 
    %   
    
    properties
        name %Name of the data set
        targetsData 
        targetsErrorList
        copyrightAndReadme % Free format description of the data and terms of use. May be printed using the printReadme method.
        statistics
        
    end
    properties(Constant=true)
        default_ComparisonTable_fields = {'rms','gdt_ts'};
    end
    
    methods
        function obj = RawData2M(name)
            obj.name = name;
            targetsErrorList = {};
        end
        
        %%Duplicate this object
        function dup = duplicate(obj,newName)
            dup = RawData2M(newName);
            dup.copyrightAndReadme = obj.copyrightAndReadme;
            for i = length(obj.targetsData()):-1:1
                dup.targetsData{i}.name = obj.targetsData{i}.name;
                dup.targetsData{i}.values = obj.targetsData{i}.values;
                dup.targetsData{i}.fields = obj.targetsData{i}.fields;
                dup.targetsData{i}.fileNames = obj.targetsData{i}.fileNames;
            end
        end
        
        function target_list = getTargets(obj,targetNames)
            target_list = RawData2M('target_list');
            for iTarget = 1:length(targetNames)
                target = obj.find_target(targetNames(iTarget));
                if (~isempty(target))
                    target_list.targetsData{length(target_list.targetsData)+1} = target;
                else
                    target_list.addTarget2ErrorList(targetNames(iTarget));
                end
            end
            
        end
        
        % createEx_addedGoapData
        % This function creates an experiment data based on
        % experimentDataBase with comparison Tables from RawData2M object.
        %
        %
        function experimentData = createEx_addRawData(obj,experimentDataBase)
            
            experimentData = ExperimentData(obj.name);
            experimentData.flagComparisonTable = 1;
            for iTarget = length(experimentDataBase.targetsData()):-1:1
                 
                 %get goap scores for models from obj
                 [target index] = obj.find_target(experimentDataBase.targetsData{iTarget}.targetName);
                 
                 exDBModelNames = regexprep(experimentDataBase.targetsData{iTarget}.fileNames,'.out.[1-9]|.xml','');
                 if (index > 0) 
                     
                     newExDataIndices = cellfun(@(x) find(strcmp(target.fileNames,x)) ,exDBModelNames,'UniformOutput', false);
                     
                     missingModelslist = find(cellfun('isempty' ,newExDataIndices));
                     if (length(missingModelslist)>0)
                         disp(['Not enough models in Comparison Table RawData for target' experimentDataBase.targetsData{iTarget}.targetName ]);
                         continue;
                     end
                     newExDataIndices = cell2mat(newExDataIndices);
                    % newExDataIndices = newExDataIndices(~cellfun('isempty' ,newExDataIndices));
                    % values = cellfun(@(i) experimentDataBase.targetsData{iTarget}.values(i,:),newExDataIndices,'UniformOutput', false);
                    % values = cell2mat(values);
                     

                     values = experimentDataBase.targetsData{iTarget}.values;
                     fields = experimentDataBase.targetsData{iTarget}.fields;
                     
                     name = experimentDataBase.targetsData{iTarget}.targetName;
                     fileNames = experimentDataBase.targetsData{iTarget}.fileNames;
                     
                     for iCompare=1:length(target.comparison_tables)
                         comparison_tables{iCompare}.CMTable = target.comparison_tables{iCompare}.CMTable(newExDataIndices,newExDataIndices);
                     end
                     
                     
                     experimentData.targetsData{iTarget} = TargetData(name,fileNames,experimentData,fields,values,comparison_tables,target.fields);
                     
                     
                 else 
                    disp(['Target ' experimentDataBase.targetsData{iTarget}.targetName ' wasn`t found on this RawData2M obj']);
                 %gdt_ts_index = Util.getIndices(obj.fields,{'gdt_ts'});
                 end
                 
                 
             end
            
             experimentData.targetsData = experimentData.targetsData(~cellfun(@isempty,experimentData.targetsData));
              
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
                if (strcmp(obj.targetsData{iTarget}.name ,target_name) == 1)
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
  
        function obj = getRawComparisonTableData(path)
            obj = RawData2M.format2_RawDataCollection(path, 'ComparisonTable', '*_similarityTable.dat', RawData2M.default_ComparisonTable_fields);
        end
        
        % format2_RawDataCollection
        % This function is for the collection of Model X Model comparison
        % data.
        function obj = format2_RawDataCollection(path, obj_name, pattern, fields)
            obj = RawData2M(obj_name);
            results_files = RawData2M.getFileList(path, pattern);
            targetNameExtractionAgent = regexprep(pattern, '*', '');            
            for iTarget = 1:length(results_files)       
                disp(iTarget);
                target_name = regexprep(results_files{iTarget}, targetNameExtractionAgent, '');
                raw_target_data = importdata([path filesep results_files{iTarget}]);
                %%TESTing the file of Raw Data is not corrupted.
                if (isempty(raw_target_data))
                    disp(['Corrupt Target file in Data Target=' target_name]);
                    obj.addTarget2ErrorList(target_name); 
                    continue;
                end
                [s inx] = unique(strcat(raw_target_data.textdata(:,1)));
                modelNamesAndIndex = raw_target_data.textdata(inx,1);
                
                modelsNumber = length(modelNamesAndIndex);
                comparisonNumber = length(raw_target_data.textdata(:,1));

                if (modelsNumber*modelsNumber ~= comparisonNumber) % check to verify all of the models comparisons have successfully finished.
                    disp(['Corrupt Target file in Data Target=' target_name]);
                    obj.addTarget2ErrorList(target_name); 
                    continue;
                end
                %comparisonTables = {};
                for iComparison = 1:comparisonNumber
                    indices = Util.getIndices(modelNamesAndIndex,raw_target_data.textdata(iComparison,:));
                    
                    if (length(indices)<2)
                        disp(['Corrupt Target file in Data Target=' target_name]);
                        obj.addTarget2ErrorList(target_name); 
                        continue;
                    end
                    
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
                target_data.comparison_tables = comparisonTables ;
                target_data.fileNames = modelNamesAndIndex;
                obj.addTarget(target_data);
            end
        end

        function fileNamesList = getFileList(path, pattern)
            fileNamesList = dir([path filesep pattern]);
            isub = [fileNamesList(:).isdir]; %# returns logical vector
            fileNamesList = {fileNamesList(~isub).name};
            fileNamesList = sort(fileNamesList);
        end  
        

        function weights = getWeights(CMT)
            cm1 = sum(CMT);
            cm2 = 1./cm1;
            cm3 = cm2./sum(cm2);
            weights = cm3;
        end
    	
	function data=loadDecoyComparisons(path)
		data.pairs = {};
		data.values = [];
		i=1;
		fd = fopen(path);
		line = fgetl(fd);
		while ischar(line)
			cols = strsplit(line);
			if (length(cols) == 4)
				data.pairs(i,:)=cols(1:2);
				data.values(i,:)=str2double(cols(3:4));
				
				i=i+1;
			end
			line = fgetl(fd);
		end
		fclose(fd);
	end
	function targetsData=loadAllCMT(path,name,pattern,fields)
        	results_files = RawData2M.getFileList(path, pattern);
		targetsData = {};
	        targetNameExtractionAgent = regexprep(pattern, '*', '');            
        	for iTarget = 1:length(results_files)       
	        	disp(iTarget);
	        	target_name = regexprep(results_files{iTarget}, targetNameExtractionAgent, '');
			targetsData{iTarget}.name=target_name;
			rawTargetData=RawData2M.loadDecoyComparisons([path filesep results_files{iTarget}])
			targetsData{iTarget}.pairs=rawTargetData.pairs;
			targetsData{iTarget}.values=rawTargetData.values;
		end
	    
	end
	function targetsCMT=applyLoadAllCMT(path)
		targetsCMT=RawData2M.loadAllCMT(path, 'ComparisonTable', '*_similarityTable.dat', RawData2M.default_ComparisonTable_fields);		
		save('targetsCMT.mat','targetsCMT');
	end
	
	function tList=getTargetsList(targetsCMT)
		tList = {};
		for iTarget = 1:length(targetsCMT)
			tList=[tList targetsCMT{iTarget}.name];
		end
	end

	function comparison_tables = createTargetCMT(experimentData,iTarget,targetsCMT,iCMTTarget)
                [nRecords v] = size(targetsCMT{iCMTTarget}.values);
		l = experimentData.targetsData{iTarget}.nModels;
		comparison_tables = {};
		for iCompare = 1:v
			comparison_tables{iCompare}.CMTable = zeros(l,l);
		end
                %exDBModelNames = regexprep(experimentData.targetsData{iTarget}.fileNames,'.scwrl.out.[1-9]|.xml','');%the file names
                exDBModelNames = regexprep(experimentData.targetsData{iTarget}.fileNames,'.out.*.xml','.pdb');%the file names
                exDBModelNames = regexprep(exDBModelNames,'[T|R][0-9][0-9][0-9][0-9].','');%the file names

		for iRecord = 1:nRecords
			%find indeces of decoys
			try
				iDecoys = Util.getIndices(exDBModelNames,targetsCMT{iCMTTarget}.pairs(iRecord,:));
			%insert value of comparison to the table
				for iCompare = 1:v
					comparison_tables{iCompare}.CMTable(iDecoys(1),iDecoys(2)) = targetsCMT{iCMTTarget}.values(iRecord,iCompare);
				end
			catch e
				%disp(getReport(e));
			end
		end
	end

	function experimentData = addComparisonTablesToExD2(experimentDataBase, targetsCMT)
            	experimentData = ExperimentData([experimentDataBase.name '_CMTs']);
	        experimentData.flagComparisonTable = 1;
		
		tList = RawData2M.getTargetsList(targetsCMT);
		for iTarget = 1:length(experimentDataBase.targetsData)	
			try
				iCMTTarget = Util.getIndices(tList,{experimentDataBase.targetsData{iTarget}.targetName});
                     		values = experimentDataBase.targetsData{iTarget}.values;
	                        fields = experimentDataBase.targetsData{iTarget}.fields;
                     
           		        name = experimentDataBase.targetsData{iTarget}.targetName;
		                fileNames = experimentDataBase.targetsData{iTarget}.fileNames;
                     
				comparison_tables = RawData2M.createTargetCMT(experimentDataBase,iTarget,targetsCMT,iCMTTarget);
                     
                    		experimentData.targetsData{iTarget} = TargetData(name,fileNames,experimentData,fields,values,comparison_tables,RawData2M.default_ComparisonTable_fields);
				
			catch e
				disp(getReport(e))
				experimentData.targetsData{iTarget} = [];
			end
		end
		experimentData.targetsData = experimentData.targetsData(~cellfun(@isempty,experimentData.targetsData));
        end
	function newExD = addComparisonTablesToExD(experimentData, pathToTable)
		RawData2M.loadDecoyComparisons(pathToTable);
	end
	function removeEmptyCMTsTargets(experimentData)
		for iTarget = 1:length(experimentData.targetsData)
			if isempty(experimentData.targetsData{iTarget}.comparisonTables)
				experimentData.targetsData{iTarget} = [];
			end
		end
		experimentData.targetsData = experimentData.targetsData(~cellfun(@isempty, experimentData.targetsData));
	end

	function experimentData = updateClusters(experimentData)
		rate = 0.95
		for iTarget = 1:length(experimentData.targetsData)
		    table1 = experimentData.targetsData{iTarget}.comparisonTables{2}.CMTable>rate;
		    neighbors = {};
		    for iModel = 1:experimentData.targetsData{iTarget}.nModels
        
		        neighbors{iModel} = find(table1(iModel,:));
        
		    end
		    experimentData.targetsData{iTarget}.comparisonTables{2}.neighbors = neighbors;
		    experimentData.targetsData{iTarget}.comparisonTables{2}.weights = RawData2M.getWeights(experimentData.targetsData{iTarget}.comparisonTables{2}.CMTable);
		    experimentData.targetsData{iTarget}.comparisonTables{2}.rate = rate;
		end
		rate = 1
		for iTarget = 1:length(experimentData.targetsData)
		    table1 = experimentData.targetsData{iTarget}.comparisonTables{1}.CMTable<=rate;
		    neighbors = {};
		    for iModel = 1:experimentData.targetsData{iTarget}.nModels
        
		        neighbors{iModel} = find(table1(iModel,:));
        
		    end
		    experimentData.targetsData{iTarget}.comparisonTables{1}.neighbors = neighbors;
		    experimentData.targetsData{iTarget}.comparisonTables{1}.rate = rate;
		end
	end 
    end
    
end

