function indices = getIndices(list,keys)
        listSize = length(keys);
        indices = ones(1,listSize);
        for iKey=1:listSize
            key = keys{iKey};
            index = find(strcmp(list,key));
            if (isempty(index)) 
                error(['key not found ' key]);
            end
            if (length(index) > 1) 
                error(['Key found more than once ',key]);
            end
            indices(iKey)= index;
        end
end
