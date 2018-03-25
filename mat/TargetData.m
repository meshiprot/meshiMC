classdef TargetData < handle
    properties
        targetName
        statistics
        fileNames
        fields
        values
        comparisonTables
        fieldsComparisonTables
        experimentData
        nModels
    end
    
    methods
        %% Constructor
        function obj = TargetData(targetName,fileNames,experimentData,fields,values,comparisonTables,fieldsComparisonTables)
            if (nargin == 0) 
                error('This is weird');
            end
            obj.targetName = targetName;            
            if ((nargin >= 5) || (nargin == 3))
                obj.fileNames  = fileNames;
                obj.experimentData = experimentData;
                if (nargin >= 5)
                    obj.fields  = fields;
                    obj.values  = values;
                    obj.nModels = size(values,1);
                end
            end
            if (nargin>5)
                obj.comparisonTables = comparisonTables;
                obj.fieldsComparisonTables = fieldsComparisonTables;
            else
                obj.comparisonTables = {};
                obj.fieldsComparisonTables = {};
            end
        end
        
        %% Column-wize methods
        function dup = duplicate(obj,experimentData)
            if (experimentData.flagComparisonTable == 0)
                dup = TargetData(obj.targetName,obj.fileNames,experimentData,obj.fields,obj.values);
            else
                dup = TargetData(obj.targetName,obj.fileNames,experimentData,obj.fields,obj.values,obj.comparisonTables,obj.fieldsComparisonTables);
            end
        end
        
        function extracted = extract(obj,fields,experimentData,flagComparisonTable)
            extracted = TargetData(obj.targetName,obj.fileNames,experimentData); 
            allFields  = obj.fields;
            if (isempty(allFields))
                error('This is weird');
            end
            indices           = Util.getIndices(allFields,fields);
            extracted.fields  = fields;
            extracted.values  = obj.values(:,indices);
            extracted.nModels = obj.nModels;
            
            if (flagComparisonTable == 1)
                extracted.fieldsComparisonTables = obj.fieldsComparisonTables;
                extracted.comparisonTables = obj.comparisonTables; 
            end
            
        end
        
        function difference = subtract(obj,fieldsToSubtract,experimentData,newFieldName)
            difference        = TargetData(obj.targetName,obj.fileNames,experimentData); 
            difference.fields = newFieldName;
            indices           = Util.getIndices(obj.fields,fieldsToSubtract);
            difference.values = obj.values(:,indices(1))-obj.values(:,indices(2));
            difference.nModels = obj.nModels;
        end 
        
        function concatenated = concatenate(obj, otherTargetData,experimentData)
            concatenated         = TargetData(obj.targetName,obj.fileNames,experimentData);
            concatenated.fields  = {obj.fields{:} otherTargetData.fields{:}};
            concatenated.values  = [obj.values otherTargetData.values];
            concatenated.nModels = obj.nModels;
        end
        
        function diff = differentFields(obj,other)
            if (length(obj.fields) ~= length(other.fields))
                diff = true;
            else 
                diff = false;
                for i = 1:length(obj.fields)
                    if (~strcmp(obj.fields{i},other.fields{i}))
                        diff = true;
                    end
                end
            end
        end
        
        function modified = applySigmoidToGdtTs(obj,sigmoid,objectiveFunction,ComparisonFlagField)
            index = Util.getIndices(obj.fields,{objectiveFunction});
            modified = TargetData(obj.targetName,obj.fileNames,obj.experimentData,obj.fields,obj.values,obj.comparisonTables,obj.fieldsComparisonTables);
            modified.values(:,index) = sigmoid.sigmoid(modified.values(:,index));
        end
        
        function extracted = extractModelsByNameList(obj, modelsNameList)
		
		%modelIndices = cellfun(@(x) find(findstr(obj.fileNames,x)) ,modelsNameList,'UniformOutput', false);
		modelIndices = cellfun(@(x) find(~cellfun(@isempty,strfind(obj.fileNames,x))) ,modelsNameList,'UniformOutput', false);
	
		modelIndices = cell2mat(modelIndices);
		
		extracted = TargetData(obj.targetName,[],obj,obj.fields,[]);
	        extracted.values = -9999.9999*ones(size(obj.values));
        	extracted.fileNames = {};
	        extracted.fileNames = obj.fileNames(modelIndices);
        	extracted.values(modelIndices,:)  = obj.values(modelIndices,:);
	        
        	extracted.fileNames(cellfun('isempty',extracted.fileNames)) = [];
	        extracted.values(extracted.values(:,1) == -9999.9999,:) = [];
        	extracted.nModels = length(extracted.fileNames);
	end
	function extracted = extractModels(obj,pattern,experimentData)
            extracted = TargetData(obj.targetName,[],experimentData,obj.fields,[]);
            extracted.values = -9999.9999*ones(size(obj.values));
            extracted.fileNames = {};
            for iModel = length(obj.fileNames):-1:1
                if (~isempty(strfind(obj.fileNames{iModel},pattern)))
                    extracted.fileNames{iModel} = obj.fileNames{iModel};
                    extracted.values(iModel,:)  = obj.values(iModel,:);
                end
            end
            extracted.fileNames(cellfun('isempty',extracted.fileNames)) = [];
            extracted.values(extracted.values(:,1) == -9999.9999,:) = [];
            extracted.nModels = length(extracted.fileNames);
        end
        
        function filterByGdt_ts(obj, minimalGdt_ts,maximalGdt_ts)
            gdtIndex = strcmp('gdt_ts',obj.fields);
            gdt = obj.values(:,gdtIndex);
            models2remove =  (gdt < minimalGdt_ts) | (gdt > maximalGdt_ts);
            obj.values(models2remove,:) = [];
            obj.fileNames(models2remove) = [];
            obj.nModels = length(obj.values(:,1));
            for iCompare = 1:length(obj.fieldsComparisonTables)
                obj.comparisonTables{iCompare}.CMTable(models2remove,:) = [];
                obj.comparisonTables{iCompare}.CMTable(:,models2remove) = [];

            end

        end
	
	
	function chooseNModels(obj, numberOfModels)
	   if (obj.nModels < numberOfModels)
		disp(['Error - Not Enough Models' obj.targetName]);
		exit;
	   end
	   models2remove = randperm(obj.nModels,obj.nModels - numberOfModels);           
	   obj.values(models2remove,:) = [];
           obj.fileNames(models2remove) = [];
           obj.nModels = length(obj.values(:,1));
           for iCompare = 1:length(obj.fieldsComparisonTables)
               obj.comparisonTables{iCompare}.CMTable(models2remove,:) = [];
               obj.comparisonTables{iCompare}.CMTable(:,models2remove) = [];
	       if isfield(obj.comparisonTables{iCompare},'neighbors')
			obj.comparisonTables{iCompare}.neighbors(models2remove) = [];
			obj.comparisonTables{iCompare}.neighbors(~isempty(obj.comparisonTables{iCompare}.neighbors));
	       end
           end
	end

	%creates a line of a table of size numberOfModelsXnumberOfModelsX3 s.t. each z-cell represents a line od [similarity(Mi,Mj),Prediction(Mi,Native),Prediction(Mj,Native)]
	function [line,trueQuality] = createSimilarityPredictionLiner(obj, numberOfModels, predFunction, trueQualityFunction, similarityFunction)
		%Choose random variables
		iModels=randi(obj.nModels,numberOfModels,1)';
		iPred = Util.getIndices(obj.fields,{predFunction});
		iQuality = Util.getIndices(obj.fields,{trueQualityFunction});
		iSimilarity = Util.getIndices(obj.fieldsComparisonTables,{similarityFunction});
		
		subCMTable = obj.comparisonTables{iSimilarity}.CMTable(iModels,iModels);
		inputTable = zeros(numberOfModels,numberOfModels,3);
		trueQuality = zeros(numberOfModels,1);
		for iTable = 1:numberOfModels
			for jTable = 1:numberOfModels
				inputTable(iTable,jTable,1) = obj.comparisonTables{iSimilarity}.CMTable(iModels(iTable),iModels(jTable));
				inputTable(iTable,jTable,2) = obj.values(iModels(iTable),iPred);
				inputTable(iTable,jTable,3) = obj.values(iModels(jTable),iPred);
			end
				trueQuality(iTable) = obj.values(iModels(iTable),iQuality);
		end
		line = inputTable(:);

	end

	function [line,trueQuality] = createFeaturesPredictionLiner(obj, numberOfModels, trueQualityFunction)
		%Choose random variables
		iModels=randi(obj.nModels,numberOfModels,1)';
		iPred = Util.getIndices(obj.fields,CoefOptimization.DefaultFeaturesToLearn);
		iQuality = Util.getIndices(obj.fields,{trueQualityFunction});
		
		inputTable = zeros(numberOfModels,length(CoefOptimization.DefaultFeaturesToLearn));
		trueQuality = zeros(numberOfModels,1);
		for iTable = 1:numberOfModels
			line(iTable,:) = obj.values(iTable,iPred);
			trueQuality(iTable) = obj.values(iModels(iTable),iQuality);
		end

	end

        function filterRandDuplicateModels(obj)
	  fileNames = obj.fileNames;
	  ModelGroupName = regexp(fileNames,'.out');
	  ModelNames = unique(cellfun(@(x,y) x(1:y-1), fileNames, ModelGroupName,'UniformOutput',false));
	  ModelsOfName = cellfun(@(x) regexp(fileNames,x),ModelNames,'UniformOutput',false);
	  %ModelsOfName = cellfun(@(x) fileNames(~cellfun(@isempty,x)),ModelsOfName,'UniformOutput',false);
	  ModelsOfName = cellfun(@(x) find(~cellfun(@isempty,x)),ModelsOfName,'UniformOutput',false);
	  NewModelsIndex = cellfun(@(x) x(randi(length(x))),ModelsOfName);
	
	
	  %models2remove = ~NewModelsIndex;
	  models2remove = setdiff(1:length(fileNames),NewModelsIndex);
          obj.values(models2remove,:) = [];
	  obj.fileNames(models2remove) = [];
	  obj.nModels = length(obj.values(:,1));
	  for iCompare = 1:length(obj.fieldsComparisonTables)
		obj.comparisonTables{iCompare}.CMTable(models2remove,:) = [];
		obj.comparisonTables{iCompare}.CMTable(:,models2remove) = [];
    	  end
        end

        %% Target_Data - this function calculates statistics between objective_function and each feature in scoreList.
        % it saves the data in target.statisctics
        % used for collection of score data.
        % statistics = Correlation, Best predicted model
        function stats = calculateStatistics(obj, objective_function, scoreList,delta_flag)
            list = [objective_function scoreList];
            indices = cellfun(@(score) Util.getIndices(obj.fields, {score}), list);
            score_values = obj.values(:,indices);
	    if (nargin > 3)
	     if (delta_flag == 1)
	 	delta_index = Util.getIndices(obj.fields, {['delta_' objective_function]});
	 	score_values(:,1) = score_values(:,1)+obj.values(:,delta_index);
	     end
            end
            
            %the fields the statistics are gathered.
            obj.statistics.ref_fields = list;
            obj.statistics.iref_fields = indices;
            
            %Correlation
            %model_target_corr = corr(score_values,'type','spearman');
            model_target_corr = corr(score_values);
            obj.statistics.target_corr = model_target_corr(1,:);
            
            %Loss - indices of the best predicted model in the target (for
            %each score function.
            obj.statistics.iMaxPredictedModelScore = cellfun(@(v,m) find(v == m) , num2cell(score_values,1), num2cell(max(score_values)),'UniformOutput',false);
            obj.statistics.iMaxPredictedModelScore = cellfun(@(x) x(1) , obj.statistics.iMaxPredictedModelScore);
            obj.statistics.MaxPredictedModelScore = score_values(obj.statistics.iMaxPredictedModelScore,1)';
            
            %Enrichment
            
            for iScore = 1:length(list)
                enrch5(iScore)   = Util.getEnrichment(score_values(:,iScore),score_values(:,1),0.95);
                enrch10(iScore)   = Util.getEnrichment(score_values(:,iScore),score_values(:,1),0.9);
                
                [perror, correlation, enrichment] = Util.getErrorAndCorr(score_values(:,iScore), score_values(:,1));

                error(iScore) = perror;  
                median_score(iScore) = median(score_values(:,iScore));
                
            end
            obj.statistics.median = median_score;
            obj.statistics.error = error;
            obj.statistics.enrichment10 = enrch10;
            obj.statistics.enrichment5 = enrch5;
            
        end
        
        
        function calculateCMTWeights(obj)
            assert(~isempty(obj.fieldsComparisonTables))
                
            for iField = 1:length(obj.fieldsComparisonTables)
                obj.comparisonTables{iField}.weights = Util.getCMTWeights(obj.comparisonTables{iField}.CMTable);
            end
        end
        
        
    end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Static methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Static = true)
        function fields = getEndFields(file)    
            xmlStruct = xml2struct(file);
            children = xmlStruct.Children;
            conformations = TargetData.select(children,'end');
            conformation = conformations(1);
            attributes = conformation.Attributes;
            children   = conformation.Children;
            evenChildren = TargetData.select(children,'all');
            fields  =  {attributes.Name  evenChildren.Name};
        end
        
        function fields = getStartFields(file)
            fields = TargetData.getFields(file,'start');
        end
        
        function fields = getFields(file, flag)
            xmlStruct = xml2struct(file);
            children = xmlStruct.Children;
            conformations = TargetData.select(children,flag);
            conformation = conformations(1);
            attributes = conformation.Attributes;
            children   = conformation.Children;
            evenChildren = TargetData.select(children,'all');
            fields  =  {attributes.Name  evenChildren.Name};
        end
            
             
        % Generator
        function targetData = getTargetData(baseDir,target,fields, selectionFlag,weadFlag, isFilesInPath)
            %selection flag = start/end/all
	    filespath = [baseDir '/' target];
            if (nargin > 5)
               if (isFilesInPath == true)
                    filesPath=baseDir;
               end
            end
	    disp (['filesPath=' filesPath])
            fields = TargetData.removeFieldsFileNamesAndValue(fields);
            targetData = TargetData(target);
            if (isdir(filesPath))
                dirString = [filesPath '/*.xml'];
                fileList = dir(dirString);
                numberOfFiles = size(fileList,1);
                errors = cell(numberOfFiles);
                targetData.fileNames = cell(numberOfFiles,1);
                data     = cell(numberOfFiles,1);
                targetData.values   = zeros(numberOfFiles,length(fields));
                targetData.nModels  = numberOfFiles;
                for iFile =1:numberOfFiles
                    fileName = fileList(iFile).name;
                    if (isempty(strfind('.N.',fileName)))
                        disp(fileName);
                        try
                            [modelData errorMessage] = TargetData.extractConformationsDataFromXml( [filesPath '/' fileName],fields, selectionFlag );
				if (isempty(errorMessage{:}))
                                vs = str2double({modelData{:}{:}.value});
                                errorMessage = TargetData.checkValues(vs,fileName,fields);
                                if (isempty(errorMessage))
                                    targetData.values(iFile,:) = vs;
                                    targetData.fileNames(iFile) = {fileName};
                                    data(iFile) = modelData;
                                else
                                    errors{iFile} = errorMessage;
                                end
                            else
                                    errors{iFile} = errorMessage;
                            end
                        catch errorMessage
                            disp(errorMessage.message)
                        end
                    end
                end
                targetData.fields    = fields;
                bad = cellfun(@(X)isempty(X),targetData.fileNames);
                targetData.values(bad,:)        = [];
                targetData.fileNames = targetData.fileNames(~bad);
                targetData.nModels   = size(targetData.values,1);
 
                
                 
                if (~isempty(targetData))
                    if (weadFlag)
                        [targetData errors] = TargetData.wead(targetData);
                    end
                    save([target '.' selectionFlag '.targetData.mat'],'targetData');
                end
                if (~isempty(errors))
                    save([target '.' selectionFlag '.errors.mat'],'errors');
                end
            end
        end
        
        function [targetData errors] = wead(targetData)
            zScores   = zscore(targetData.values);
            badValues = (zScores>6) | (zScores < -6);
            errors = {};
            while(~isempty(find(badValues,1)))
                outliers = sum(badValues,2)>0;
                errors = {errors{:} ['(outliers ' {targetData.fileNames(outliers)} ')']};
                targetData.values    = targetData.values(~outliers,:);
                targetData.nModels   = size(targetData.values,1);
                targetData.fileNames = targetData.fileNames(~outliers,:);
                zScores   = zscore(targetData.values);
                badValues = (zScores>6) | (zScores < -6);
            end
        end
        
        
        function errorMessage = checkValues(values, fileName,fields)
            nanElelents = isnan(values);
            if (~isempty(find(nanElelents,1)))
                disp(fields(nanElelents));
                errorMessage = ['NaN in ' fileName];
            else
                if (~isempty(find(isinf(values),1)))
                    errorMessage = ['Inf in ' fileName];
                else
                    errorMessage = '';
                end
            end
            
        end
        
        function [conformationsData errorMessages] = extractConformationsDataFromXml( file,fields, selectionFlag )    
            xmlStruct = xml2struct(file);
            children = xmlStruct.Children;
            conformations = TargetData.select(children,selectionFlag);
            conformationsData = cell(1,length(conformations));
            errorMessages     = cell(1,length(conformations));
            for iConf = 1:length(conformations)
                conformation = conformations(iConf);
                [conformationData errorMessage] = TargetData.getConformationData(conformation,fields,selectionFlag);
                conformationsData(iConf) = {conformationData};
                errorMessages(iConf) = {errorMessage};
                conformationsData(iConf) = {conformationsData(iConf)};
            end
        end
        
        function [conformationData errorMessage] = getConformationData(conformation,fields,selectionFlag)
            attributes = conformation.Attributes;
            children   = conformation.Children;
            evenChildren = TargetData.select(children,'all');
            childrenAttributes = [evenChildren.Attributes];
            names  =  {attributes.Name  evenChildren.Name};
            childrenValues = {attributes.Value childrenAttributes.Value};
            newNames  = cell(1,length(fields));
            newValues = cell(1,length(fields));
            fileIndex = find(strcmp('fileName',names));
            if (isempty(fileIndex))
                errorMessage = 'Cannot find file name.';
                conformationData = [];
                return;
            end
            fileName = childrenValues{fileIndex(1)};
            
            errorMessage = TargetData.check(childrenValues,selectionFlag, fileName);
            if (~isempty(errorMessage))
                conformationData = [];
                return;
            end
            for i = 1:length(fields)
                j = find(strcmp(fields{i},names));
                if (length(j)<=0)
                    errorMessage = [ 'cannot find field ' fields{i} ' in ' fileName];
                    conformationData = [];
                    return;
                end
                if  (length(j) > 1)
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %     The three lines below were used in meshi.8.17/CASP8_9 in which by mistake there were two energy fields
                     %                if (strcmp(fields{i},'energy'))
                     %                   j = j(2);
                     %                end
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    errorMessage = [errorMessage 'field ' fields{i} ' ocures more than once in ' filename]; %#ok<AGROW>
                    conformationData = [];
                    return;
                end
                newNames{i} = fields{i};
                newValues{i} = childrenValues{j};
            end
            conformationData = struct('name',newNames,'value',newValues);
        end

        function selected = select(array,flag)
            n = length(array);% Weird array of odd length with empty odd positions.
            if (strcmp(flag,'end'))
                selected = struct(array(n-1));
            elseif (strcmp(flag,'start'))
                selected = struct(array(2));
            elseif (strcmp(flag,'all'))
                selected = array(2:2:n);
            else
                error(['unknown flag ' flag]);
            end
        end
        
        function errorMessage = check(values,selectionFlag, fileName)
            if(~strcmp(selectionFlag,'end'))
                errorMessage = '';
            else
                if (isempty(find(strcmp('MCM_END',values),1)))
                    errorMessage = ['Cannot find MCM_END in ' fileName];
                else
                    errorMessage = '';
                end
            end
        end
        
        function fields = removeFieldsFileNamesAndValue(fields)
            fields(strcmp('fileName',fields) | strcmp('value',fields))= [];
        end
    end

    
end

