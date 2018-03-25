


function collectTestData(name)
	cd(name);
	test_folder = dir();
	isub = [test_folder(:).isdir]; %# returns logical vector
	param_group_number_folders = {test_folder(isub).name};
	param_group_number_folders(ismember(param_group_number_folders,{'.','..','logs'})) = [];

 	ca = {};
	
		

		
	flist = dir(['configuration.*']);	
	flist = {flist(:).name};
		
		%disp(['the current directory: ' pwd]);
		%disp(flist);
	i=1;
	for file_number = 1:length(flist)
		load([flist{file_number}]);
		for ci = 1:length(configuration)
			ca{i} =  configuration{ci};	
			i = i+1;
		end
	end

	
	cd .. ;
	save(['ca.mat'],'ca');
		
		

end
