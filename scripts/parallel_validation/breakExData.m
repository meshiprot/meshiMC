


addpath /fastspace/users/siditom/MESHI_cloud/MY_MATLAB13_tm;


load experimentData.mat;
exBase = experimentData;
tList = experimentData.getTargetNames;
fid = fopen('ExTargets','wt')
for i=1:length(tList)
  fprintf(fid, [tList{i} '\n']);
end
fclose(fid);
