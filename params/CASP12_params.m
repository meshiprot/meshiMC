% exPath='../experimentData.mat';
% CASP12 params
exPath = '/home/cluster/users/siditom/data/meshi_db_v2.0/experimentData.meshiDBv2.0.mat';
numberOfCoefsRange = 10;
normalizers = {'length' 'contacts8' 'contacts15'};
nonZeroCoefsNum = 15;

cvJobs = 50;
numOfConfigsInArray = 20;

MCMCSteps = 3000;
% MCMCSteps = 1000;
initialTemp = 0.001;
% seed = 1;

learningFraction = 0.5;
objectiveFunction = 'gdt_ts';
deltaObjective = 'delta_gdt_ts';
configEnergy = ConfigEnergyErrorsAndMedianCorrelation; 
featuresToLearn={'energy','bondEnergy','angleEnergy','planeEnergy','ramachandranSidechain','ramachSTD','cooperativeZRamachandranSidechain','cooperativeZstdRamachandranSidechain','ramachandranCore','compositePropensity','cooperativeZPropensity','cooperativeZstdPropensity','atomicPairwisePMFSumma','summaStd','excludedVolume','cooperativeZSumma','cooperativeSummaPolar','cooperativeSummaNonPolar','cooperativeSummaNeutral','cooperativeSummaPolarNN_OO','cooperativeSummaPolarBb','cooperativeZstdSumma','cooperativeStdSummaPolar','cooperativeStdSummaNonPolar','cooperativeStdSummaNeutral','cooperativeStdSummaPolarNN_OO','cooperativeStdSummaPolarBb','solvationEnergy','solvationSCpolar','solvationSCcarbon','solvationBBpolar','solvationBBcarbon','solvationHB','solvationBuriedHB','solvationSTD','solvationEntropy','atomEnvironmentEnergy','atomEnvironmentPropensity','atomEnvironmentEnergySS','atomEnvironmentPropensitySS','solvateEnergy','solvateSCpolar','solvateSCcarbon','solvateBBpolar','solvateSTD','hydrogenBonds','hydrogenBondsPairs','hydrogenBondsAnglesHOC','hydrogenBondsAnglesOHN','contactsAndSASA','contacts12','contacts14','contacts14Core','SASAratio','rg','N_RGhSS','hSS','bSS','cSS','hSShCoil','hSSbSS','hSSbCoil','hSScSS','hSScCoil','goap','dfire','goap_ag','contacts8','contacts11','contacts15','contactsHr','flatRamachEnergy','samudralaEnergy','secondaryStructureFraction','ssCompatibility','sasaCompatibility','coverage','scwrl','one','length','nAtoms'}; % CASP12 scorefunctiontraining, with no iterativefunctions
