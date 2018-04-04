% exPath='../experimentData.mat';
% CASP13 params

exPath = '/home/cluster/users/siditom/data/meshi_db_v2.0/experimentData.meshiDBv2.0.mat';
numberOfCoefsRange = 10;
normalizers = {'length' 'contacts8' 'contacts15' 'coverage'};
nonZeroCoefsNum = 15;

cvJobs = 500;
numOfConfigsInArray = 2;

MCMCSteps = 3000;
% MCMCSteps = 1000;
initialTemp = 0.001;
% seed = 1;

learningFraction = 0.5;
objectiveFunction = 'gdt_ts8';
deltaObjective = 'delta_gdt_ts';
configEnergy = ConfigEnergyErrorCorrelationEnrichment;
featuresToLearn = {'energy','distanceConstraints','bondEnergy','angleEnergy','planeEnergy','outOfPlaneEnergy','ramachandranSidechain','ramachSTD','cooperativeZRamachandranSidechain','cooperativeZstdRamachandranSidechain','ramachandranCore','compositePropensity','cooperativeZPropensity','cooperativeZstdPropensity','atomicPairwisePMFSumma','summaStd','excludedVolume','cooperativeZSumma','cooperativeSummaPolar','cooperativeSummaNonPolar','cooperativeSummaNeutral','cooperativeSummaPolarNN_OO','cooperativeSummaPolarBb','cooperativeZstdSumma','cooperativeStdSummaPolar','cooperativeStdSummaNonPolar','cooperativeStdSummaNeutral','cooperativeStdSummaPolarNN_OO','cooperativeStdSummaPolarBb','solvationEnergy','solvationSCpolar','solvationSCcarbon','solvationBBpolar','solvationBBcarbon','solvationHB','solvationBuriedHB','solvationSTD','solvationEntropy','bbPolarN','bbCarbonN','scPolarN','scCarbonN','atomEnvironmentEnergy','atomEnvironmentEnergyBackbone','atomEnvironmentEnergyPolar','atomEnvironmentEnergyPositive','atomEnvironmentEnergyNegative','atomEnvironmentEnergyAromatic','atomEnvironmentEnergyAliphatic','atomEnvironmentPropensity','atomEnvironmentPropensityBackbone','atomEnvironmentPropensityPolar','atomEnvironmentPropensityPositive','atomEnvironmentPropensityNegative','atomEnvironmentPropensityAromatic','atomEnvironmentPropensityAliphatic','atomEnvironmentEnergySS','atomEnvironmentEnergyCoil','atomEnvironmentEnergyHelix','atomEnvironmentEnergySheet','atomEnvironmentPropensitySS','atomEnvironmentPropensityCoil','atomEnvironmentPropensityHelix','atomEnvironmentPropensitySheet','solvateEnergy','solvateSCpolar','solvateSCcarbon','solvateBBpolar','solvateSTD','hydrogenBonds','hydrogenBondsPairs','hydrogenBondsAnglesHOC','hydrogenBondsAnglesOHN','contactsAndSASA','contacts12','contacts14','contacts14Core','SASAratio','rg','N_RGhSS','EhSShCoil','EhSSbSS','EhSSbCoil','EhSScSS','EhSScCoil','hSS','bSS','cSS','hSShCoil','hSSbSS','hSSbCoil','hSScSS','hSScCoil','goap','dfire','goap_ag','conservationContacts8','contacts8','conservationRgRatio','conservation_H_RgRatio','conservationContacts11','contacts11','conservationContacts15','contacts15','conservationContactsHr','contactsHr','conservationRgRatioHr','tetherEnergy','flatRamachEnergy','tetherAll','samudralaEnergy','ssCompatibility','SaturatedSsCompatibility','helixCompatibility','SheetCompatibility','sasaCompatibility','saturatedSasaCompatibility','secondaryStructureFraction','saturatedSsFraction','HelixFraction','SheetFraction','coverage','sScoverage','nAtoms','scwrl','length','one','deepCNF3Compatibility','deepCNF8Compatibility'}; % CASP13
