


exD = ExperimentData('assembly');
tList = dir('ExTargets_VD');
isub = [tList(:).isdir]; %# reiturns logical vector                                 ExTargets_VD
isub = [~isub];
tList = {tList(isub).name}; 
load experimentData.mat
originEX = experimentData;
 
for iTarget = 1:length(tList)
	disp(tList{iTarget});
	load (['ExTargets_VD/' tList{iTarget}]);
	exD.targetsData{iTarget} = experimentData.targetsData{1};
	iOriginTarget = Util.getIndices(originEX.getTargetNames(),{exD.targetsData{iTarget}.targetName});
	exD.targetsData{iTarget}.comparisonTables = originEX.targetsData{iOriginTarget}.comparisonTables;
	exD.targetsData{iTarget}.fieldsComparisonTables = originEX.targetsData{iOriginTarget}.fieldsComparisonTables;
end


experimentData = exD.duplicate('Validation');
save('validationData.mat','experimentData');

experimentData = exD;
save('validationData_allData.mat','experimentData');

