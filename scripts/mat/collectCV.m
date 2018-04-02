%********************************************************************************************
% This file specifies the parameters for the CoefOptimization.repeatedCoefOptimization program.
% It will be added to the CopyrightAndReadme field of the resulted experimentData object.
% ********************************************************************************************
function cv = collect_configurations(nFolds,nParams)
 tic
	pv
	metrics = {'error','target_corr','enrichment5','enrichment10'};
	cv = {};
	for iMetric = 1:length(metrics)
		cv.(metrics{iMetric}) = NaN(nFolds,nParams);
	end
	cv.loss = NaN(nFolds,nParams);
	for iFold = 1:nFolds
		for iParam = 1:nParams
			load(['fold' num2str(iFold) '.' num2str(iParam) filesep 'vd.mat']);
			validationData.calculateStatistics(objectiveFunction,{'predictionWeightedMedian'});
			validationData.filterByGdt_ts(0,0.9999,10);                                                                  
			validationData.calculateStatistics('gdt_ts',{'predictionWeightedMedian'},1);                                 

			for iMetric = 1:length(metrics)
		                cv.(metrics{iMetric})(iFold,iParam) = median(validationData.statistics.(metrics{iMetric})(:,2));
		        end
			cv.loss(iFold,iParam) = median(validationData.statistics.MaxPredictedModelScore(:,1) - validationData.statistics.MaxPredictedModelScore(:,2));
		end
	end
	save('cv.mat','cv');
 toc
 end
