

function vd = GModelsVD(validationData,minimalGdt_ts,maximalGdt_ts)
    vd = ExperimentData('GModelsVD');
    for iTarget = 1:length(validationData.targetsData)
        
        
        target = validationData.targetsData{iTarget}.duplicate(validationData);
        target.filterByGdt_ts(minimalGdt_ts,maximalGdt_ts);
        [w l] = size(target.values);
        if (w>1) 
            vd.targetsData{length(vd.targetsData)+1} = target;
        end
    end

end