classdef Mexico2013
    
    properties
    end
    
    methods(Static=true)
        function A
            Mexico2013.A2('atomicPairwisePMFSumma');
        end
        function A2(feature)
            correlations = [];
            cd d:\Work\selectionProject\testSet\start
            load validationData.mat;
            Mexico2013.origGdt(validationData);
            correlations = Mexico2013.A1(validationData,feature,1,'.',{'Server models'},correlations,[-4000 11000],[0.7 0.60 0.2 0.35],'A1','A2','A3');
            
            cd d:\Work\selectionProject\testSet\all
            load validationData.mat
            Mexico2013.origGdt(validationData);
            correlations = Mexico2013.A1(validationData,feature,2,'.', {'Server models' 'Minimized models'},correlations, ...
                                                                        [-6000 1000],[0.2 0.25 0.4 0.35],[],'A4','A5');
            
            indices = Util.getIndices(validationData.targetsData{1}.fields,{feature,'contacts8'});
            for i = 1:length(validationData.targetsData)
                validationData.targetsData{i}.values(:,indices(1)) = validationData.targetsData{i}.values(:,indices(1)).*...
                                                                     validationData.targetsData{i}.values(:,indices(2));
            end
             correlations = Mexico2013.A1(validationData,feature,3,'.', {'Server models' 'Minimized models', 'Normalized'},correlations, [-20000 2000],...
                                          [0.2 0.25 0.4 0.35],[],[],'A6');
           
        end
        
        function origGdt(validationData)
                        indices = Util.getIndices(validationData.targetsData{1}.fields,{'gdt_ts','delta_gdt_ts'});
            for i = 1:length(validationData.targetsData)
                validationData.targetsData{i}.values(:,indices(1)) = validationData.targetsData{i}.values(:,indices(1))-...
                                                                     validationData.targetsData{i}.values(:,indices(2));
                if (sum(validationData.targetsData{i}.values(:,indices(1))>1)~=0)
                    display('*******************************************')
                    display(validationData.targetsData{i});
                    j = validationData.targetsData{i}.values(:,indices(1))>1;
                    validationData.targetsData{i}.fileNames{j};
                    display('*******************************************')
                end
            end                
        end
        
        function correlations = A1(validationData,feature,n,format,insetLabel,correlations,ylim, insetPosition,fileName1,fileName2,fileName3) 
            [gdt,energy] = validationData.plotFields({'gdt_ts' feature},format,100,1,1);
            Mexico2013.styleFigure('Test set (CASP10 models)','GDT\_TS (%)', 'Energy (Summa & Levitt 2007)',[0 0 0],[0 100], ylim);
            Mexico2013.printFigure(gcf,fileName1);
                
            mainAxes = gca;
            c = corr([gdt energy]);
            c = c(1,2);
            disp(['Overall correlation = ' num2str(c)]);
            
             energy_gdt = validationData.extract({'gdt_ts' feature},'temp');
            for i = length(energy_gdt.targetsData):-1:1
                c = corr([energy_gdt.targetsData{i}.values(:,1) energy_gdt.targetsData{i}.values(:,2)]);
                correlations(i,n) = c(1,2);
            end
            
            input('press')
            [sc ,ic] = sort(correlations(:,n));
            disp([sc(1) sc(end)]);
            disp({validationData.targetsData{ic(1)}.targetName validationData.targetsData{ic(end)}.targetName}); 
             axes(mainAxes);
             ph = plot(energy_gdt.targetsData{ic(1)}.values(:,1)*100, energy_gdt.targetsData{ic(1)}.values(:,2),'go');
             set(ph,'MarkerSize',7,'MarkerFaceColor',[0 1 0])
             ph = plot(energy_gdt.targetsData{ic(end)}.values(:,1)*100, energy_gdt.targetsData{ic(end)}.values(:,2),'ro');
             set(ph,'MarkerSize',7,'MarkerFaceColor',[1 0 0])
            Mexico2013.printFigure(gcf,fileName2);
            hold
            input('press');
            inset = axes('position',insetPosition);            
            boxplot(inset, correlations,'notch','on','labels',insetLabel);
            %set(inset,'Color','none');
            Mexico2013.styleFigure('','', 'Correlation',[1 1 1], [0 (n+1)],[-1 1]);
            Mexico2013.printFigure(gcf,fileName3);
            set(inset,'Visible','off');
            input('press');
            axes(mainAxes);
        end
        
        function styleFigure(titelMain, titleX,titleY,colorX,limX,limY,ax)
            fig = gcf;
            if (nargin == 6)
                ax  = gca;
            end
            set(fig,'color','w');
            set(ax,'box','off');
            hx = xlabel(ax,titleX);
            hy = ylabel(ax,titleY);
            ht = title(ax,titelMain);
            set(findobj(ax,'Type','text'),'FontSize',15,'FontName','Times New Roman', 'Rotation',-15,'VerticalAlignment','Cap','HorizontalAlignment','left');
            set(ax,'xcolor',colorX);
            set([hx hy],'FontSize',15); 
            set(ht,'FontSize',20)
            set([hx hy ht],'FontName','Times New Roman');
            if (~isempty(limX))
                xlim(limX);
            end
            if ((nargin == 6) && (~isempty(ylim)))
                ylim(limY)
            end
        end
        
        function B
            features = {'contacts8' 'contacts15' 'atomicPairwisePMFSumma'  'ramachandranSidechain' 'compositePropensity'...
                        'hydrogenBondsPairs' 'solvateEnergy' ...
                        'cooperativeSummaNonPolar'  'solvationSTD'  'ramachSTD' 'bSS'   };
            featureNames = {'contacts8'  'contacts15' 'SummaLevitt' 'Ramachandran' 'Propensity'...
                        'HB-Pairs' 'solvation' ...
                        'coopSummaNP'  'solvationSTD'  'ramachSTD' 'bSS'   };
            cd d:\Work\selectionProject\testSet\all
            load validationData.mat;
            Mexico2013.origGdt(validationData);
            gdtIndex = Util.getIndices(validationData.targetsData{1}.fields,{'gdt_ts'});
            contactsIndex = Util.getIndices(validationData.targetsData{1}.fields,{'contacts8'});
            indices = Util.getIndices(validationData.targetsData{1}.fields,features);
            for i = length(features):-1:1
                for j = length(validationData.targetsData):-1:1
                    c = corr([validationData.targetsData{j}.values(:,gdtIndex) validationData.targetsData{j}.values(:,indices(i))]);
                    correlations(1,j,i) = c(1,2);
                    c = corr([validationData.targetsData{j}.values(:,gdtIndex) validationData.targetsData{j}.values(:,indices(i))./...
                                                                               validationData.targetsData{j}.values(:,contactsIndex)]);
                    correlations(2,j,i) = c(1,2);
                    c = corr([validationData.targetsData{j}.values(:,gdtIndex) validationData.targetsData{j}.values(:,indices(i)).*...
                                                                               validationData.targetsData{j}.values(:,contactsIndex)]);
                    correlations(3,j,i) = c(1,2);
                end
            end
                medians(1,:) = median(correlations(1,:,:));
                medians(2,:) = median(correlations(2,:,:));
                medians(3,:) = median(correlations(3,:,:));
                [~, iMax] = max(abs(medians));
                for i = length(features):-1:1
                   for j = length(validationData.targetsData):-1:1
                      bestCorrelations(j,i) = correlations(iMax(i),j,i);
                   end
                end
            boxplot(bestCorrelations,'notch','on','labels',featureNames);
            Mexico2013.styleFigure('Sample features','','Correlation',[1 1 1],[0 13],[-1 1]);
            set(gca,'Position',[0.11 0.16 0.9 0.75]);
            Mexico2013.printFigure(gcf,'B');
        end
        
        
        function C(target,fig1,fig2,fig3,fig4,fig5,fig6) 
            cd d:\Work\selectionProject\testSet\all
            load validationData.mat;
            figure;
            a = axes;
            Mexico2013.styleFigure('Weighted Histogram of 3000 independant scores', ...
                         'Score (%)','Frequncy (%)',[0 0 0],[0 1]);
            hold
            input('press 1');
            Mexico2013.C1(a,target,1,fig1,fig2,fig3,[0 0 1],{Mexico2013.getModelName(validationData.targetsData{target}.fileNames{1}) 'GDT\_TS' '<Score>'});
            input('press 2');
            Mexico2013.C1(a,target,4,[],[],fig4,[1 0 0],{Mexico2013.getModelName(validationData.targetsData{target}.fileNames{1}) 'GDT\_TS' '<Score>' ...
                                                     Mexico2013.getModelName(validationData.targetsData{target}.fileNames{4}) 'GDT\_TS' '<Score>'});
            input('press 3');
            Mexico2013.C1(a,target,3,[],[],fig5,[0 1 0],{Mexico2013.getModelName(validationData.targetsData{target}.fileNames{1}) 'GDT\_TS' '<Score>' ...
                                                     Mexico2013.getModelName(validationData.targetsData{target}.fileNames{4}) 'GDT\_TS' '<Score>' ... 
                                                     Mexico2013.getModelName(validationData.targetsData{target}.fileNames{3}) 'GDT\_TS' '<Score>'});
            inset = axes('Position',[0.55 0.5 0.3 0.4]);
            indices = Util.getIndices(validationData.targetsData{1}.fields,{'gdt_ts' 'predictionWeightedMedian'});
            plot(validationData.targetsData{target}.values(:,indices(1))*100, validationData.targetsData{target}.values(:,indices(2))*100,'c.')
            c = corr([validationData.targetsData{target}.values(:,indices(1))*100, validationData.targetsData{target}.values(:,indices(2))*100]);
            disp(['Correlation = ' num2str(c(1,2))]);
            e = Util.getEnrichment(validationData.targetsData{target}.values(:,indices(2)),validationData.targetsData{target}.values(:,indices(1)));
            disp(['Enrichment = ' num2str(e)]);
            hold
            plot([75 100],[75 100],'k:')
            plot(validationData.targetsData{target}.values(1,indices(1))*100, validationData.targetsData{target}.values(1,indices(2))*100,'*')
            plot(validationData.targetsData{target}.values(4,indices(1))*100, validationData.targetsData{target}.values(4,indices(2))*100,'r*')
            plot(validationData.targetsData{target}.values(3,indices(1))*100, validationData.targetsData{target}.values(3,indices(2))*100,'g*')
            Mexico2013.styleFigure('', ...
                         'GDT\_TS (%)','<Score> (%)',[0 0 0],[]);
            Mexico2013.printFigure(gcf,fig6);
        end
        
        function C1(ha,target,model,fileName1,fileName2, fileName3,color,names)
            cd d:\Work\selectionProject\testSet\all
            load validationData.mat;
            validationData.name
            Mexico2013.origGdt(validationData);

            %disp(validationData.targetsData{target}.fileNames{model});
            predictions = validationData.targetsData{target}.statistics.predictions(model,:);
            weights     = validationData.targetsData{target}.statistics.predictions(model,:);
            gdtIndex = Util.getIndices(validationData.targetsData{1}.fields,{'gdt_ts'});
            gdt = validationData.targetsData{target}.values(model,gdtIndex);
            wm = Configuration.weightedMedian(predictions,weights);
            disp(['GDT = ' num2str(gdt) ' weighted median = ' num2str(wm)]);
            %[sp ip] = sort(predictions);
            totalWeights = sum(weights);
            for i = 50:-1:1
                bin(i) = sum(weights((predictions>(i-1)/50) & (predictions<=(i)/50)))*100/totalWeights;
                binCenters(i) = (i-0.5)/50;
            end
            %binSize = totalWeights/40000;
            %edges = 0:binSize:1;
            %N= histc(predictions,edges);
            %plot(cumsum(weights(ip)),sp);
            plot(ha,binCenters,bin,'Color',color);
            Mexico2013.printFigure(gcf,fileName1);
            input('press')
            plot(ha,[gdt gdt],[bin(round(50*gdt)) bin(round(50*gdt))],'Marker','O','Color',color,'MarkerFaceColor',color,'MarkerSize',10);
            Mexico2013.printFigure(gcf,fileName2);
            plot(ha,[wm wm],[bin(round(50*wm)) bin(round(50*wm))],'Marker','P','Color',color,'MarkerFaceColor',color,'MarkerSize',10);
            lh = legend(names,'Location', 'NorthWest','Box','off','EdgeColor',[1 1 1]);
            Mexico2013.printFigure(gcf,fileName3);
            
        end
        
        function D
            dirs = {'d:\Work\selectionProject\testSet\all' 'd:\Work\selectionProject\testSet\fm_tbm' ...
                    'd:\Work\selectionProject\Casp10Results_tatiana\qmean&observed' 'D:\Work\selectionProject\Casp10Results_tatiana\dfire&observed(energy+)Casp10'};
            for i = 4:-1:1
                cd(dirs{i});
                load validationData.mat;
                Mexico2013.origGdt(validationData);
                indices = Util.getIndices(validationData.targetsData{1}.fields,{'gdt_ts' 'predictionWeightedMedian'});
                for j = length(validationData.targetsData):-1:1
                    gdt   = validationData.targetsData{j}.values(:,indices(1));
                    score = validationData.targetsData{j}.values(:,indices(2));
                    c = corr([gdt score]);
                    correlations(j,i) = c(1,2);
                    diff              = gdt-score;
                    errors(j,i)       = sqrt(mean(diff.*diff));
                    enrichments(j,i)   = Util.getEnrichment(score,gdt);
                end
            end
            erh = subplot(1,3,1);
            coh = subplot(1,3,2);
            enh = subplot(1,3,3);
       
            boxplot(erh,errors(:,1),'notch','on','labels',{'All targets'});
            Mexico2013.styleFigure('Error', '','(%)',[0 0 0],[],[],erh)
            boxplot(coh,correlations(:,1),'notch','on','labels',{'All targets'});
            Mexico2013.styleFigure('Correlation', '','',[0 0 0],[],[],coh)
            boxplot(enh,enrichments(:,1),'notch','on','labels',{'All targets'});
            Mexico2013.styleFigure('5% enrichment', '','',[0 0 0],[],[],enh)
            Mexico2013.printFigure(gcf,'D1');
            input('press');
            boxplot(erh,errors(:,1:2),'notch','on','labels',{'All targets' 'FM + TBM'});
            Mexico2013.styleFigure('Error', '','(%)',[0 0 0],[],[],erh)
            boxplot(coh,correlations(:,1:2),'notch','on','labels',{'All targets' 'FM + TBM'});
            Mexico2013.styleFigure('Correlation', '','',[0 0 0],[],[],coh)
            boxplot(enh,enrichments(:,1:2),'notch','on','labels',{'All targets' 'FM + TBM'});
            Mexico2013.styleFigure('5% enrichment', '','',[0 0 0],[],[],enh)
            Mexico2013.printFigure(gcf,'D2');
            input('press');
            boxplot(erh,errors(:,1:3),'notch','on','labels',{'All targets' 'FM + TBM' 'qMean'});
            Mexico2013.styleFigure('Error', '','(%)',[0 0 0],[],[],erh)
            boxplot(coh,correlations(:,1:3),'notch','on','labels',{'All targets' 'FM + TBM' 'qMean'});
            Mexico2013.styleFigure('Correlation', '','',[0 0 0],[],[],coh)
            boxplot(enh,enrichments(:,1:3),'notch','on','labels',{'All targets' 'FM + TBM' 'qMean'});
            Mexico2013.styleFigure('5% enrichment', '','',[0 0 0],[],[],enh)
            Mexico2013.printFigure(gcf,'D3');
            input('press');
            boxplot(erh,errors(:,1:3),'notch','on','labels',{'All targets' 'FM + TBM' 'qMean'});
            Mexico2013.styleFigure('Error', '','(%)',[0 0 0],[],[],erh)
            boxplot(coh,correlations,'notch','on','labels',{'All targets' 'FM + TBM' 'qMean' 'dFire'});
            Mexico2013.styleFigure('Correlation', '','',[0 0 0],[],[],coh)
            boxplot(enh,enrichments,'notch','on','labels',{'All targets' 'FM + TBM' 'qMean' 'dFire'});
            Mexico2013.styleFigure('5% enrichment', '','',[0 0 0],[],[],enh)
            Mexico2013.printFigure(gcf,'D4');
        end
        
        function E(servers,fileName)
            cd d:\Work\selectionProject\testSet\fm_tbm
            load validationData.mat;
            Mexico2013.origGdt(validationData);
            for i = 2:-1:1
                for j = 3:-1:1
                    ahs(i,j) = subplot(2,3,(i-1)*3+j);
                    serverData = validationData.extractModels(servers{i,j});
                    indices = Util.getIndices(validationData.targetsData{1}.fields,{'gdt_ts' 'predictionWeightedMedian'});
                    for l = length(serverData.targetsData):-1:1
                        table = [];
                        target = serverData.targetsData{l};
                        target.fileNames{:};
                        if (target.nModels> 0)
                            table(:,1) = target.values(:,indices(1)); % real gdt
                            table(:,2) = (length(target.fileNames):-1:1)'; % server rank
                            r = (length(target.fileNames):-1:1)';
                            [s, si] = sort(target.values(:,indices(2)),'descend');
                            table(si,3) = r; %my rank
                            table(:,4) = target.values(:,indices(2)); % my score
                            c = corr(table,'Type','Spearman');
                            correlation(l,:) = c(1,2:3);
                        end
                    end
                    boxplot(ahs(i,j),correlation,'notch','on','labels',{servers{i,j} 'Our score'});
                    Mexico2013.styleFigure('', '','Spearman correlation',[0 0 0],[],[],ahs(i,j))
                end
            end
            set(findobj(gcf,'Type','text'),'FontSize',10);
            input('press')
            Mexico2013.printFigure(gcf,'E');
        end            
            
        
        function Z
            cd d:\Work\selectionProject\testSet\fm_tbm
            load validationData.mat;
            Mexico2013.origGdt(validationData);
            errors = [];
            allInterdeciles = [];
            allEntropies = [];
            hold
            for i = 1:length(validationData.targetsData)
                if(mod(i,10) == 0)
                    disp(i)
                end
                [medians, interdeciles, entropies] = Configuration.weightedMedian(validationData.targetsData{i}.statistics.predictions,...
                                                                     validationData.targetsData{i}.statistics.predictions);
                errors = [errors;  abs(medians - validationData.targetsData{i}.values(:,8))];  
                allInterdeciles = [allInterdeciles; interdeciles]; 
                allEntropies = [allEntropies; entropies];
%                plot(errors,allEntropies,'.');
%                input('press');
            end
            h1 = subplot(1,2,1);
            h2 = subplot(1,2,2);
            plot(h1,errors,allInterdeciles,'.');
            corr([errors,allInterdeciles])
            plot(h2,errors,allEntropies,'.');
            corr([errors,allEntropies])
            
        end
        
        
        function printFigure(fig,fileName)
            if (~isempty(fileName))
                figure(fig);
                cd d:\Work\selectionProject\Mexico2013
                set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 25 15])
                print('-dtiff',fileName, '-r600');
            end
        end
        
        function partialCorrelation(feature)
            features = {'contacts8' 'contactsHr' 'contacts15' 'atomicPairwisePMFSumma' 'hydrogenBonds' 'ramachandranSidechain' 'compositePropensity'...
                        'hydrogenBondsPairs' 'solvateEnergy' ...
                        'cooperativeSummaNonPolar'  'solvationSTD'  'ramachSTD' 'bSS'  'hSS' };
            cd d:\Work\selectionProject\testSet\all
            load validationData.mat;
            Mexico2013.origGdt(validationData);
            gdtIndex = Util.getIndices(validationData.targetsData{1}.fields,{'gdt_ts'});
            featureIndex = Util.getIndices(validationData.targetsData{1}.fields,{feature});
            indices = Util.getIndices(validationData.targetsData{1}.fields,features);
            for i = length(features):-1:1
                for j = length(validationData.targetsData):-1:1
                    c = partialcorr([ validationData.targetsData{j}.values(:,featureIndex) validationData.targetsData{j}.values(:,indices(i))],...
                                  validationData.targetsData{j}.values(:,gdtIndex));
                    correlations(j,i) = c(1,2);
                end
            end
             boxplot(correlations,'notch','on','labels',features);
            Mexico2013.styleFigure('Sample features','','Correlation',[1 1 1],[1 1 1],[0 14],[-1 1]);
        end
        
        function modelName = getModelName(fileName)
            dotIndices = strfind(fileName,'.');
            temp = fileName((dotIndices(1)+1):(dotIndices(3)-1));
            underScoreIndex = strfind(temp,'_');
            modelName = [temp(1:underScoreIndex-1) '\' temp(underScoreIndex:end)];
        end
        
        
        function qmean(validationData)
            cd 'D:\Work\selectionProject\Casp10Results_tatiana\qmean&observed';
            cd 'D:\Work\selectionProject\Casp10Results_tatiana\dfire&observed(energy+)Casp10'
            for i = 1:length(validationData.targetsData)
                target = validationData.targetsData{i};
                l = length(target.values(1,:));
                targetName = target.targetName(1:8);
                fid = fopen([targetName '.txt'],'r');
                line = 'xxx';
                while(ischar(line))
                    line = fgetl(fid);
                    if (ischar(line)) 
                        pdbIndex = strfind(line,'pdb');
                        fileName = line(1:pdbIndex+2);
                        index = -1;
                        for j = 1:length(target.fileNames)
                            %target.fileNames{j}(1:end-14)
                            if (strcmp(fileName(1:end-4),target.fileNames{j}(1:end-14)))
                                gdtAndScore = sscanf(line(pdbIndex+3:end),'%f %f');
                                index = j;
                            end
                        end
                        if (index ~= -1)
                            %if (abs(gdtAndScore(1)-target.values(index,8))>0.01)
                                %target
                                %j
                                %[gdtAndScore(1) target.values(index,8)]'
                            %   error('this is weird');
                            %end
                            target.values(index,l) = gdtAndScore(2);
                        end
                    end
                end
                %plot(gdt,score,'.');
                %input('press');
            end
        end

    end
  
    
end

