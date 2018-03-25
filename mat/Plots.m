
classdef Plots < handle


    methods(Access=public, Static=true)
        
        % plot predicted values with added weights
        % iTarget - the target number to plot
        % iModel - the model number to plot
        function values = plot1data(validationData, iTarget, iModel)
            
            if (~Plots.check(validationData,iTarget,iModel))
                return;
            end
            %%collected weight of predictors.
            values(:,1) = validationData.targetsData{iTarget}.statistics.predictions(iModel,:);
            values(:,2) = validationData.targetsData{iTarget}.statistics.weights(iModel,:);

            values = sortrows(values,1);
            for i = 2:length(values(:,2))
                values(i,2) = values(i,2) + values(i-1,2);
            end
            values(:,2) = values(:,2)/values(end,2);
        end
        function addGroupScore(experimentData)
		for iTarget = 1:length(experimentData.targetsData)
			valIndex = length(experimentData.targetsData{iTarget}.fields) +1;
			experimentData.targetsData{iTarget}.fields = [experimentData.targetsData{iTarget}.fields 'GroupScore']
			for iDecoy = 1:length(experimentData.targetsData{iTarget}.fileNames)
				iStr = findstr(experimentData.targetsData{iTarget}.fileNames{iDecoy},'TS');
				experimentData.targetsData{iTarget}.values(iDecoy,valIndex) = 1/str2num(experimentData.targetsData{iTarget}.fileNames{iDecoy}(iStr+2));
			end
		end
	end 
	% receives validationData and checks for each GROUP checks how many times their decoys were the 1st Decoy of the TARGET. Also, check how many times MESHI-score chosen the group to be its first choise.
	function [histo] = firstDecoyByGroup(experimentData)
		%First - create the group names participating in this DB.
		experimentData.filterByGdt_ts(0,0.99999,10);
		experimentData.calculateStatistics('gdt_ts',{'predictionWeightedMedian'},1);
		histo = [];
		histo.groups = [];
		group_names = {'BAKER', 'MULTICOM', 'RaptorX', 'Pcons','HHpred','FALCON','Jiang','FFAS','RBO','chuo','Distill','SAM-','Bhageerath','Phyre','panther','3D-JIGSAW','raptor','intfold','coma','confuzz','mufold','loopp'};
		group_names = lower(group_names)
		for iTarget = 1:length(experimentData.targetsData)
			
			histo.groups = [histo.groups; lower(experimentData.targetsData{iTarget}.fileNames)];
			histo.groups = regexprep(histo.groups,'_ts.*.xml','');
			%histo.groups = regexprep(histo.groups,'.*.\.N\..*xml','Native');	
			histo.groups = regexprep(histo.groups,'_al.*.xml','');	
			for iGroup = 1:length(group_names)
				 histo.groups = regexprep( histo.groups,[group_names{iGroup} '.*'],group_names{iGroup});
			end
			histo.groups = regexprep(histo.groups,'-server.*.xml','');	
			histo.groups = unique(histo.groups);
		end
		% Count how many times each group was successful	
		histo.group_count = zeros(length(histo.groups),1)';
		for iTarget = 1:length(experimentData.targetsData)
			firstDecoyName = lower(experimentData.targetsData{iTarget}.fileNames(experimentData.statistics.iMaxPredictedModelScore(iTarget,1)));
			firstDecoyName = regexprep(firstDecoyName,'_ts.*.xml','');
			firstDecoyName = regexprep(firstDecoyName,'_al.*.xml','');
			for iGroup = 1:length(group_names)
                                 firstDecoyName = regexprep(firstDecoyName,[group_names{iGroup} '.*'],group_names{iGroup});
                        end
			firstDecoyName = regexprep(firstDecoyName,'-server.*.xml','');
			groupIndex = Util.getIndices(histo.groups, {firstDecoyName});
			histo.group_count(groupIndex) = histo.group_count(groupIndex) +1;
		end	
		histo.meshi_count = zeros(length(histo.groups),1)';
		for iTarget = 1:length(experimentData.targetsData)
			firstDecoyName = lower(experimentData.targetsData{iTarget}.fileNames(experimentData.statistics.iMaxPredictedModelScore(iTarget,2)));
			firstDecoyName = regexprep(firstDecoyName,'_ts.*.xml','');
			firstDecoyName = regexprep(firstDecoyName,'_al.*.xml','');
			for iGroup = 1:length(group_names)
                                 firstDecoyName = regexprep(firstDecoyName,[group_names{iGroup} '.*'],group_names{iGroup});
                        end
			firstDecoyName = regexprep(firstDecoyName,'-server.*.xml','');
			groupIndex = Util.getIndices(histo.groups, {firstDecoyName});
			histo.meshi_count(groupIndex) = histo.meshi_count(groupIndex) +1;
		end	
	end

 
        
        % This function return histo. which contains the data for a plot
        % with a grouped weights for each range of predicted score between
        % [1/k-1, 1/k]
        % iTarget - the target number to plot
        % iModel - the model number to plot
        function [histo points]= plot2data(validationData,iTarget,iModel,k)
            %%collected weight of predictors in range.
            if (~Plots.check(validationData,iTarget,iModel))
                return;
            end
            values(:,1) = validationData.targetsData{iTarget}.statistics.predictions(iModel,:);
            values(:,2) = validationData.targetsData{iTarget}.statistics.weights(iModel,:);
            values = values(values(:,2)>1,:);
            
            values = sortrows(values,1);
            j = 1;
            for i = 1:k
                histo(i,1) = i/k;
                histo(i,2) = 0;
                while ((j < length(values(:,2))) && (values(j,1) < i/k))
                    histo(i,2) = histo(i,2) + values(j,2);
                    j=j+1;
                end

            end
            histo(:,2) = histo(:,2)/sum(values(:,2));
            
            
            points.predScore.x = validationData.targetsData{iTarget}.values(iModel,Util.getIndices(validationData.targetsData{iTarget}.fields, {'predictionWeightedMedian'}));
            points.maxPredScore.x = validationData.targetsData{iTarget}.values(iModel,Util.getIndices(validationData.targetsData{iTarget}.fields, {'predictionMaxWeightedMedian'}));
            points.score.x = validationData.targetsData{iTarget}.values(iModel,1);
            
            points.maxPredScore.y = histo(max(find(histo(:,1)<points.maxPredScore.x)),2);
            points.predScore.y = histo(max(find(histo(:,1)<points.predScore.x)),2);
            points.score.y = histo(max(find(histo(:,1)<points.score.x)),2);
            
        end
        
        
        
        function f = plot2(values,points,add)
            if (add)
                hold on;
            end
            plot(values(:,1),values(:,2),'b',1,max(values(:,2))+0.1);
            hold on;         
            plot(points.score.x ,points.score.y,'*r');
            plot(points.predScore.x ,points.predScore.y,'*b');
            f = plot(points.maxPredScore.x ,points.maxPredScore.y,'*g');
            hold off;
        end
        
        
        function plotW_H(COrun,iTarget,iModel)
            [values points]= Plots.plot2data(COrun.VD,iTarget,iModel,100);
            
            [f sc ratio] = Plots.Plot_wsVwd(COrun,iTarget,iModel,0);
            
            f2 = figure
            
            ax1 = subplot(2,1,1);
            Plots.plot2(values,points,0);
            
            ax2 = subplot(2,1,2);
            
            
            set(sc,'parent',ax2) 
            close(f)
            
            set(f2,'position',[500 200 500 700])
            
        end
        
        

        
        %% all - is validation data with predicted score values for each of the models.

        function ComparisonTable(all,headline)
            
            %colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD><B>',text,'</B></TD></TR> </table></html>'];
            %colergen = @(color,text) ['<html><table border=0 width=400 ><TR><TD><B>',text,'</B></TD></TR> </table></html>'];
            
            M_Loss = median(all.statistics.MaxPredictedModelScore - meshgrid(all.statistics.MaxPredictedModelScore(:,1),1:length(all.statistics.ref_fields))');
            M_Loss = Plots.bold_max(M_Loss);
            Loss_IQR = Plots.IqrBold(all.statistics.MaxPredictedModelScore - meshgrid(all.statistics.MaxPredictedModelScore(:,1),1:length(all.statistics.ref_fields))','max')


            M_Corr = median(all.statistics.target_corr);
            Corr_IQR = Plots.IqrBold(all.statistics.target_corr,'max')
            
            
            M_ENR5 = median(all.statistics.enrichment5);
            M_ENR5 = Plots.bold_max(M_ENR5);
            ENR5_IQR = Plots.IqrBold(all.statistics.enrichment5,'max')
            
            M_ENR10 = median(all.statistics.enrichment10);
            %M_ENR10 = Plots.bold_max(M_ENR10);
            ENR10_IQR = Plots.IqrBold(all.statistics.enrichment10,'max')
            
            M_ERROR = median(all.statistics.error);    
            Error_IQR = Plots.IqrBold(all.statistics.error,'min');
            Error_IQR{3} = 'Not Relevant'
            %M_ERROR = Plots.bold_min(M_ERROR);
            %M_ERROR{3} = 'Not Relevant';         
            %f = figure('Position',[500 500 500 500]);

            % Column names and column format
            columnname = all.statistics.ref_fields(2:end);
            columnname = {'MESHI','GOAP','DFire','QMean'};
            %columnname = {'<html><center>MESHI</center></html>', ...
            %           '<html><center>GOAP</center></html>', ...
            %           '<html><center>DFire</center></html>', ...
            %           '<html><center>QMean</center></html>'};
            columnname = {'<html><center>MESHI+GOAP</center></html>', ...
                        '<html><center>MESHI</center></html>', ...
                       '<html><center>GOAP</center></html>', ...
                       '<html><center>QMean</center></html>'};

            %rowname = {'IQR Corr';'Median Corr';'Median Loss';'Median ENR5';'Median ENR10'};
            rowname = {'Median Error';'Median (IQR) Corr';'Median Loss';'Median ENR5'};

            
            
            
            
            %str = colergen('#FF0000',M_Loss{1});
            %M_Loss = {colergen('#FF0000',M_Loss{2}) M_Loss{3:end}};
            %disp(M_Loss);
            
            figurePos = [650 400 600 300];
            hFigure = figure('Position',figurePos ,'Color','w');
            set(gcf,'numbertitle','off','name','Protein Prediction Score Statistics') % See the help for GCF

            
            boxPos = [100 10 434 20];
            color=get(hFigure,'Color');
            %MyBox = uicontrol('style','text','BackgroundColor',color);
            MyBox = uicontrol('style','text','BackgroundColor','w');
            set(MyBox,'String','Table1 - QA methods staticstics comparisons CASP10 database');
            set(MyBox,'Position',boxPos);
            
            boxPos = [70 135 434 30];
            color=get(hFigure,'Color');
            %MyBox = uicontrol('style','text','BackgroundColor',color);
            MyBox = uicontrol('style','text','BackgroundColor','w');
            %set(MyBox,'String','MESHI meta parameters: CA size = 100 (randomly chosen from a 10 cv run), MaxSteps = 1000, temp = 0.0001, Coefs = 15');
            set(MyBox,'String',headline);
            set(MyBox,'Position',boxPos);
            
            
            
            tablePos = [40 35 549 94];
            hTable = uitable(hFigure,'Data',[Error_IQR; Corr_IQR; Loss_IQR; ENR5_IQR],...
                                     'RowName',rowname,...
                                     'ColumnName',columnname,...
                                     'ColumnWidth',{90},...
                                     'Position',tablePos);
            
            
            
        end
        
        function mIQRData = IqrBold(data,order)
            %colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD><B>',text,'</B></TD></TR> </table></html>'];
            colergen = @(color,text) ['<html><table border=0 width=400 ><TR><TD><B>',text,'</B></TD></TR> </table></html>'];
           
            M_data = median(data);
            
            
            idata=0;
            %find min index
            if (strcmp(order,'min')==1)
                idata = find(min(M_data(2:end))==M_data) - 1;
            end
            if (strcmp(order,'max')==1)
                idata = find(max(M_data(2:end))==M_data) - 1;
            end
            
            %combain IQR ans data
            data_IQR = iqr(data);
            for index = 1:length(M_data)-1
                
                M_data_iqr{index} = [num2str(M_data(index+1)) ' (' num2str(data_IQR(index+1)) ')'];
            end
            
            
            
            M_data_iqr{idata} = colergen('#FF0000',M_data_iqr{idata});
            
            mIQRData = M_data_iqr;

        end
        
        
        % Plot to discribe the weight of true score of a model and a
        % predicted score of a model accross all configurations in COrun.
        function [f fsc ratio] = Plot_wsVwd(COrun,iTarget,iModel,visible)
            wd = COrun.VD.targetsData{iTarget}.statistics.weights(iModel,:);


            for iCon = 1:length(COrun.CA)
                wScore = COrun.CA{iCon}.sigmoid.backSigmoid(COrun.VD.targetsData{iTarget}.values(iModel,1));
                wScore = COrun.CA{iCon}.sigmoid.derivative(wScore);%the weight (derivative of sigmoid) of the objective function 
                ws(iCon) = wScore;
            end

            %ws = ws/sum(ws);
            [ws' wd'];
            ratio = wd.*ws;
            c = abs(COrun.VD.targetsData{iTarget}.statistics.predictions(iModel,:) - COrun.VD.targetsData{iTarget}.values(iModel,1));
            for iColor = 1:length(c)
                c(iColor) = floor(c(iColor)*20)/20;
            end
            
            
            f = figure;
            
            if (nargin>3)
                if (visible==0)
                    set(f,'visible','off');
                end
            end
            %new = abs(COrun.VD.targetsData{iTarget}.statistics.predictions(iModel,find(ratio==max(ratio))) - COrun.VD.targetsData{iTarget}.values(iModel,1))
            fsc = scatter(wd,ws,25,c,'filled');
            
            %legend(['Model Gdt = ' num2str(COrun.VD.targetsData{iTarget}.values(iModel,1))],['Mean(W) = ' num2str(mean(wd))])
            boxPos = [100 350 200 40];
           
            %MyBox = uicontrol('style','text','BackgroundColor',color);
            sortedWeights = sort(wd);
            topSorted = sortedWeights(floor(length(sortedWeights)*0.9));
            mmm = mean(wd(find(wd>topSorted)));
            MyBox = uicontrol('style','text','BackgroundColor','w');
            iWMedian = Util.getIndices(COrun.VD.targetsData{iTarget}.fields,{'predictionWeightedMedian'});
            iWMaxMedian = Util.getIndices(COrun.VD.targetsData{iTarget}.fields,{'predictionMaxWeightedMedian'});
            set(MyBox,'String',['Model Gdt = ' num2str(COrun.VD.targetsData{iTarget}.values(iModel,1)) ' Mean(W) = ' num2str(mean(wd)) ' Mean(W(m)) = ' num2str(mmm) ' WMedianScore = ' num2str(COrun.VD.targetsData{iTarget}.values(iModel,iWMedian))  ' WMaxMedianScore = ' num2str(COrun.VD.targetsData{iTarget}.values(iModel,iWMaxMedian))   ]);
            set(MyBox,'Position',boxPos);
            
            
            m = find(COrun.VD.targetsData{iTarget}.statistics.weights(iModel,:) == max(COrun.VD.targetsData{iTarget}.statistics.weights(iModel,:)));
            m1 = find(COrun.VD.targetsData{iTarget}.statistics.weights(iModel,:) > 0.98 );
            m2 = find(COrun.VD.targetsData{iTarget}.statistics.weights(iModel,:) < 1.02 );
            m3 = setdiff(m1,m2);
            m4 = setdiff(m1,m3);

            avg = mean(COrun.VD.targetsData{iTarget}.statistics.predictions(iModel,m4));
            avg = abs(avg - COrun.VD.targetsData{iTarget}.values(iModel,1));

            new = [abs(COrun.VD.targetsData{iTarget}.statistics.predictions(iModel,m) - COrun.VD.targetsData{iTarget}.values(iModel,1)) COrun.VD.targetsData{iTarget}.statistics.weights(iModel,m)];
            old = abs(COrun.VD.targetsData{iTarget}.values(iModel,86) - COrun.VD.targetsData{iTarget}.values(iModel,1));

            old_error = abs(COrun.VD.targetsData{iTarget}.values(iModel,86) - COrun.VD.targetsData{iTarget}.values(iModel,1));
            new_error = abs(COrun.VD.targetsData{iTarget}.values(iModel,90) - COrun.VD.targetsData{iTarget}.values(iModel,1));
        end
        
        %plot the ROC curve of iTarget in the validation data COrun.VD
        % errorRate - the error number which above it is failure and
        % beneath it is success.
        % predictionRef - is the name of the prediction that the error is
        % calculated from.
        function [iW iFail iSuccess rate] = plotROC_weights(COrun, iTarget, errorRate, predictionRef)
            iPred = Util.getIndices(COrun.VD.targetsData{1}.fields,{predictionRef})
            
            
            for iModel = 1:COrun.VD.targetsData{iTarget}.nModels
                wd = COrun.VD.targetsData{iTarget}.statistics.weights(iModel,:);
                sortedWeights = sort(wd);
                topSorted = sortedWeights(floor(length(sortedWeights)*0.9));
                avg10Weights(iModel) = mean(wd(find(wd>topSorted)));
            end
            %avg10Weights = median(COrun.VD.targetsData{iTarget}.statistics.weights')'
            i = 1;
            for iW = 3:-0.1:0.5
                iModels = find(avg10Weights>iW);
                error = abs(COrun.VD.targetsData{iTarget}.values(iModels,1) - COrun.VD.targetsData{iTarget}.values(iModels,iPred));
                currentModels = COrun.VD.targetsData{iTarget}.values(iModels,:);
                iFail(i) = sum(error>errorRate);
                iSuccess(i) = sum(error<errorRate);
                
                
                if (iModels>0)
                    
                    index = find(error<errorRate);
                    if (length(index)>0)
                        maxSuccessModels(i) = max(currentModels(index,1));
                    else
                        maxSuccessModels(i) = 0;
                    end
                    if (iW >2)
                    index = index 
                    iModels = iModels'
                    end
                    bestModelScore(i) = max(COrun.VD.targetsData{iTarget}.values(iModels,1));
                    [m index] = max(currentModels(:,iPred));
                    [x index] = max(m == currentModels(:,iPred));
                    
                    bestPredictedModelScore(i) = currentModels(index,1);
                    
                    
                    indexSet = find(error<errorRate);
                    if (length(indexSet) >0 && length(index)>0 )
                        
                        if (ismember(indexSet,index))
                            bestPredictedModelScoreIndex(i) = 1;
                        else
                            bestPredictedModelScoreIndex(i) = 0;
                        end
                    else
                        bestPredictedModelScoreIndex(i) = 0;
                    end
                        
                else
                    bestModelScore(i) = 0;
                    bestPredictedModelScore(i) = 0;
                     bestPredictedModelScoreIndex(i) = 0;
                end
                i = i+1;
            end
            iW = (3:-0.1:0.5)';
            iFail = iFail';
            iSuccess = iSuccess';
            rate = (iSuccess'./iFail')';
            
            [iW iSuccess iFail rate bestModelScore' bestPredictedModelScore' bestPredictedModelScoreIndex' maxSuccessModels']  
            plot(iFail'/COrun.VD.targetsData{iTarget}.nModels,iSuccess'/COrun.VD.targetsData{iTarget}.nModels);
        end
        
         %plot the ROC curve of all targets in the validation data COrun.VD
        % errorRate - the error number which above it is failure and
        % beneath it is success.
        % predictionRef - is the name of the prediction that the error is
        % calculated from.
        function [iW success fail rate] = plotROC_Aweights(COrun, errorRate, predictionRef)
            iPred = Util.getIndices(COrun.VD.targetsData{1}.fields,{predictionRef})
            
            totalModels = 0;
            for iTarget = 1:length(COrun.VD.targetsData)
                [iW iF iS] = Plots.plotROC_weights(COrun,iTarget, errorRate, predictionRef);
                iFail(iTarget,:) = iF;
                iSucess(iTarget,:) = iS;
                totalModels = totalModels + COrun.VD.targetsData{iTarget}.nModels;
            end
            fail = sum(iFail)'/totalModels;
            success = sum(iSucess)'/totalModels;
            
            
            rate = success./fail;
            [iW success fail rate]
            plot(fail,success);
        end
        
        function Plot_MetaParamExploration()
            f = figure('Position',[450 200 800 700],'Color','w');
            
            %Error plots
            Plots.plot_metaP('Median Error',1,'ConfigArrayPTest_Coefs15_MStep1000RES.mat',...
                        2,'end','No. of Predictors', {'10' '100' '500' '1.0e+3'});
            Plots.plot_metaP('Median Error',2,'MaxStepsPTest_ConfArr10RES.mat',...
                        2,5,'No. of Optimization Steps', { '10' '100' '1000' '2000' });
            Plots.plot_metaP('Median Error',3,'CoefNum_PTest10_1000RES.mat',...
                        2,6,'No. of Featurs', { '5' '10' '15' '20' '25'});

            %Correlation plots
            Plots.plot_metaP('Median Corr',4,'ConfigArrayPTest_Coefs15_MStep1000RES.mat',...
                        2,'end','No. of Predictors', {'10' '100' '500' '1.0e+3'});
            Plots.plot_metaP('Median Corr',5,'MaxStepsPTest_ConfArr10RES.mat',...
                        2,5,'No. of Optimization Steps', { '10' '100' '1000' '2000' });
            Plots.plot_metaP('Median Corr',6,'CoefNum_PTest10_1000RES.mat',...
                        2,6,'No. of Featurs', { '5' '10' '15' '20' '25'});

            boxPos = [225 235 400 20];
            color=get(f,'Color');
            %MyBox = uicontrol('style','text','BackgroundColor',color);
            MyBox = uicontrol('style','text','BackgroundColor','w');
            %set(MyBox,'String','MESHI meta parameters: CA size = 100 (randomly chosen from a 10 cv run), MaxSteps = 1000, temp = 0.0001, Coefs = 15');
            set(MyBox,'String','Figure 1 - Meta parameters exploration, error and correlation box plot statistics.');
            set(MyBox,'Position',boxPos);
            
        end
        function Plot_MetaParamExploration_OnlyOneParam()
            f = figure('Position',[100 150 600 200],'Color','w');
            
            %Error plots
            Plots.plot_metaP_onlyOneParameter('error','Median Error',1,'ConfigArrayPTest_Coefs15_MStep1000RES.mat',...
                        [1,2,3,5],'No. of Predictors', {'1' '10' '100' '1000'});
            
            %Correlation plots
            Plots.plot_metaP_onlyOneParameter('corr','Median Corr',2,'ConfigArrayPTest_Coefs15_MStep1000RES.mat',...
                        [1,2,3,5],'No. of Predictors', {'1' '10' '100' '1000'});
            

            boxPos = [80 200 20 100];
            color=get(f,'Color');
            %MyBox = uicontrol('style','text','BackgroundColor',color);
            MyBox = uicontrol('style','text','BackgroundColor','w');
            %set(MyBox,'String','MESHI meta parameters: CA size = 100 (randomly chosen from a 10 cv run), MaxSteps = 1000, temp = 0.0001, Coefs = 15');
            set(MyBox,'String','Figure 1 - Meta parameters exploration, error and correlation box plot statistics.');
            set(MyBox,'Position',boxPos);
            
        end
        
        function Plot_StartEnd_DataComparison()
            f = figure('Position',[100 150 600 200],'Color','w');
            
            %Error plots
            Plots.plot_metaP_onlyOneParameter('error','Median Error',1,'boxplotMatrix.mat',...
                        [1:2],'Energy Minimizatrion', {'without' 'with'});
            
            %Correlation plots
            Plots.plot_metaP_onlyOneParameter('corr','Median Corr',2,'boxplotMatrix.mat',...
                        [1:2],'Energy Minimizatrion', {'without' 'with'});
            

            boxPos = [80 200 20 100];
            color=get(f,'Color');
            %MyBox = uicontrol('style','text','BackgroundColor',color);
            MyBox = uicontrol('style','text','BackgroundColor','w');
            %set(MyBox,'String','MESHI meta parameters: CA size = 100 (randomly chosen from a 10 cv run), MaxSteps = 1000, temp = 0.0001, Coefs = 15');
            set(MyBox,'String','Figure 1 - Meta parameters exploration, error and correlation box plot statistics.');
            set(MyBox,'Position',boxPos);
            
        end
        
        function Plot_ErrorCorr_DataComparison(boxplotMatrix)
            f = figure('Position',[100 150 600 200],'Color','w');
            
            %Error plots
            Plots.plot_metaP_onlyOneParameter('error','Median Error',1,boxplotMatrix,...
                        [1:4],'', {'' '' '' ''});
            
            %Correlation plots
            Plots.plot_metaP_onlyOneParameter('corr','Median Corr',2,boxplotMatrix,...
                        [1:4],'', {'' '' '' ''});
            

            boxPos = [80 200 20 100];
            color=get(f,'Color');
            %MyBox = uicontrol('style','text','BackgroundColor',color);
            MyBox = uicontrol('style','text','BackgroundColor','w');
            %set(MyBox,'String','MESHI meta parameters: CA size = 100 (randomly chosen from a 10 cv run), MaxSteps = 1000, temp = 0.0001, Coefs = 15');
            set(MyBox,'String','Figure 1 - Meta parameters exploration, error and correlation box plot statistics.');
            set(MyBox,'Position',boxPos);
            
        end
        
        function e = PlotTargetModels_ScorePred(validationData,iTarget)
            if ((length(validationData.targetsData) < iTarget) || (iTarget < 0))
                e = -1;
                return;
            end
            f = figure('Color','w');
            index = Util.getIndices(validationData.targetsData{iTarget}.fields,{'gdt_ts' 'predictionWeightedMedian'});
            
            %enrichment box
            x = 0.9;
            objective = (validationData.targetsData{iTarget}.values(:,index(1)));
            sortedObjective     = sort(objective);
            topObjective        = sortedObjective(round(length(objective)*x))
            
            prediction = validationData.targetsData{iTarget}.values(:,index(2));
            sortedPredictions   = sort(prediction);
            topPrediction       = sortedPredictions(round(length(objective)*x));
            topIndices          = find(prediction>=topPrediction);

            
            topObjectiveIndices = find(objective<topObjective);

            
            enr_points = setdiff(topIndices,topObjectiveIndices);

            
            hold on;
            scatter(validationData.targetsData{iTarget}.values(:,index(1)),validationData.targetsData{iTarget}.values(:,index(2)),'b','filled');
            scatter(validationData.targetsData{iTarget}.values(enr_points,index(1)),validationData.targetsData{iTarget}.values(enr_points,index(2)),'g','filled');

            ylabel('Score');
            xlabel('Gdt ts');
            ylim([0 1]);
            xlim([0 1]);
            rectangle()
            
            
            %rectangle('Position',[topObjective topPrediction (1-topObjective) (1-topPrediction)])
            e = 0;
        end
        



    end
    methods(Access=private, Static=true)
        
        % 1 - transfering list to a str cell array
        % 2 - bolding the max value in the array
        function newList = bold_max(list)
            %colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD><B>',text,'</B></TD></TR> </table></html>'];
            colergen = @(color,text) ['<html><table border=0 width=400 ><TR><TD><B>',text,'</B></TD></TR> </table></html>'];
            
            
            list = list(2:end);
        
            iMax_list = find(max(list)==list);
            list = arrayfun(@num2str, list, 'unif', 0);
            if (length(iMax_list) == 1)
                list{iMax_list} = colergen('#FF0000',list{iMax_list});
            else
                for iMax = iMax_list 
                    list{iMax} = colergen('#FF0000',list{iMax});
                end
            end
            newList = list;
        end
        function newList = bold_min(list)
            %colergen = @(color,text) ['<html><table border=0 width=400 bgcolor=',color,'><TR><TD><B>',text,'</B></TD></TR> </table></html>'];
            colergen = @(color,text) ['<html><table border=0 width=400 ><TR><TD><B>',text,'</B></TD></TR> </table></html>'];
            
            
            list = list(2:end);
        
            iMin_list = find(min(list)==list);
            list = arrayfun(@num2str, list, 'unif', 0);
            if (length(iMin_list) == 1)
                list{iMin_list} = colergen('#FF0000',list{iMin_list});
            else
                for iMin = iMin_list 
                    list{iMin} = colergen('#FF0000',list{iMin});
                end
            end
            newList = list;
        end
        
        function plot_metaP_onlyOneParameter(type,yLabel,subNumber,data,range,meta_param_name,xticks)
           load(data);

           data_matrix = boxplotMatrix.errorWeightedMedian(:,range);
           if (strcmp(type ,'error')) 
               data_matrix = boxplotMatrix.errorWeightedMedian(:,range);
               
           else if (strcmp(type ,'corr')) 
                   data_matrix = boxplotMatrix.corrWeightedMedian(:,range);
               end
           end


           subplot(1,2,subNumber);

           boxplot(data_matrix,xticks ,'notch','on');
           ylabel(yLabel);
           if (strcmp(type ,'error')) 
               
               %BreakPlot([1:4],data_matrix,0.18,0.3,'RPatch');
               %ylim([0.1 0.3])
           end
           xlabel(meta_param_name);
           text_h = findobj(gca, 'Type', 'text');


            for cnt = 1:length(text_h)
                set(text_h(cnt),    'FontSize', 8);
            end
        end
        
        
        function plot_metaP(type,subNumber,data,rangeS,rangeE,meta_param_name,xticks)
           load(data);
           if (strcmp(rangeE,'end'))
               rangeE = length(boxplotMatrix.errorWeightedMedian(1,:));
           end
           data_matrix = boxplotMatrix.errorWeightedMedian(:,rangeS:rangeE);
           if (strcmp(type ,'error')) 
               data_matrix = boxplotMatrix.errorWeightedMedian(:,rangeS:rangeE);
           else if (strcmp(type ,'corr')) 
                   data_matrix = boxplotMatrix.corrWeightedMedian(:,rangeS:rangeE);
               end
           end


           subplot(3,3,subNumber);

           boxplot(data_matrix,xticks ,'notch','on');
           ylabel(type);
           %ylim([0.1 0.2])
           xlabel(meta_param_name);
           text_h = findobj(gca, 'Type', 'text');


            for cnt = 1:length(text_h)
                set(text_h(cnt),    'FontSize', 8);
            end
        end
        
        
        function maxWeights_predictions = getMaxWeightValues(validationData, iTarget)
            
            for iModel = 1:validationData.targetsData{iTarget}.nModels
                values = Plots.plot2data(validationData, iTarget,iModel,100);
                maxWeights_predictions(iModel) = values(find(values(:,2)== max(values(:,2)),1));
            end
            
        end
        
        function model_list = getkScoreModelList(validationData, iTarget,maxScore)
            model_list = find(validationData.targetsData{iTarget}.values(:,1)<maxScore);
        end

        
        
        
        function vrf = check(validationData, iTarget,iModel)
            if (iTarget <= length(validationData.targetsData) && (iModel <= validationData.targetsData{iTarget}.nModels))
                vrf = 1;
                return;
            end
            vrf = 0;
        end
        
        
    end
end
