% SetDataPathandFilenames

% Debug mode.
ExperimentInfo.debug = 0;
% Slice visualization step.
ExperimentInfo.sliceStep = 4;
% Path to datasets.
ExperimentInfo.dataPath = '/home/makrogianniss/Data/MRIConvertOutput/';
% Algorithm parameters for:
% Leg selection.
ExperimentInfo.LegSelectionThreshold = 300 / 1500;  %400/1500
% Shading correction.
ExperimentInfo.tophatTransform = 0;
ExperimentInfo.strelSizeFactorTopHat = 16; %8;
% Number of scales in deformable segmentation.
ExperimentInfo.nScales = 3;
% Structural element size for definition of initial mask.
ExperimentInfo.snakemaskseSize = 10;
% Remove the bone.
ExperimentInfo.removeBone = 1;
ExperimentInfo.boneremovalThreshold = 100;
ExperimentInfo.strelSizeFactorBoneRemoval = 128;
% Clustering.
ExperimentInfo.clusterNonZeroValuesOnly = 1;
% Post-clustering bone removal.
ExperimentInfo.areaThreshold = 20;
% Distance metric learning.
DistanceLearningMethods = {'', 'PCA', 'ISOMAP', 'Laplacian', 'LLE', 'LTSA', 'MDS'};
ExperimentInfo.DistanceLearningMethod = DistanceLearningMethods{1};
ExperimentInfo.OutputDimensionality = 2;
% Clustering algorithms.
ClusteringMethods = {'K-MEANS', 'FCM', 'MeanShift'};
ExperimentInfo.ClusteringMethod = ClusteringMethods{1};
ExperimentInfo.MeanShiftBandwidth = 200;  % 200

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ExperimentInfo.SubjectID = 'BLSA_1703_07';  % hip prosthesis (?)
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_1703_07\311250577_301_T1_3D_80_SLICE_20091113\BLSA_1703_07_311250577_301_T1_3D_80_SLICE_20091113.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_1703_07\311250577_401_WS_T1_3D_80_SLICE_20091113\BLSA_1703_07_311250577_401_WS_T1_3D_80_SLICE_20091113.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_1703_07\311250577_501_FS_T1_3D_80_SLICE_20091113\BLSA_1703_07_311250577_501_FS_T1_3D_80_SLICE_20091113.hdr';
% ExperimentInfo.ProcessedLeg = 'right';

% ExperimentInfo.SubjectID = 'BLSA_1927_05';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_1927_05\311852000_201_T1_3D_80_SLICE_20091120\BLSA_1927_05_311852000_201_T1_3D_80_SLICE_20091120.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_1927_05\311852000_301_WS_T1_3D_80_SLICE_20091120\BLSA_1927_05_311852000_301_WS_T1_3D_80_SLICE_20091120.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_1927_05\311852000_401_FS_T1_3D_80_SLICE_20091120\BLSA_1927_05_311852000_401_FS_T1_3D_80_SLICE_20091120.hdr';
% ExperimentInfo.ProcessedLeg = 'left';

% % BLSA_5012_18
% ExperimentInfo.SubjectID = 'BLSA_5012_18';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_5012_18/310571886_701_T1_3D_80_SLICE_20091105/BLSA_5012_18_310571886_701_T1_3D_80_SLICE_20091105.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_5012_18/310571886_501_WS_T1_3D_80_SLICE_20091105/BLSA_5012_18_310571886_501_WS_T1_3D_80_SLICE_20091105.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_5012_18/310571886_601_FS_T1_3D_80_SLICE_20091105/BLSA_5012_18_310571886_601_FS_T1_3D_80_SLICE_20091105.hdr';
% ExperimentInfo.ProcessedLeg = 'left';
% ExperimentInfo.boneremovalThreshold = 130;


% ExperimentInfo.SubjectID = 'BLSA_5189_14';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_5189_14\311674832_301_T1_3D_80_SLICE_20091118\BLSA_5189_14_311674832_301_T1_3D_80_SLICE_20091118.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_5189_14\311674832_401_WS_T1_3D_80_SLICE_20091118\BLSA_5189_14_311674832_401_WS_T1_3D_80_SLICE_20091118.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_5189_14\311674832_501_FS_T1_3D_80_SLICE_20091118\BLSA_5189_14_311674832_501_FS_T1_3D_80_SLICE_20091118.hdr';
% ExperimentInfo.ProcessedLeg = 'left';

% ExperimentInfo.SubjectID = 'BLSA_530_14';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_5302\310487230_201_T1_3D_80_SLICE_20091104\BLSA_5302_310487230_201__T1_3D_80_SLICE_20091104.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_5302\310487230_301_WS_T1_3D_80_SLICE_20091104\BLSA_5302_310487230_301__WS_T1_3D_80_SLICE_20091104.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_5302\310487230_401_FS_T1_3D_80_SLICE_20091104\BLSA_5302_310487230_401__FS_T1_3D_80_SLICE_20091104.hdr';
% ExperimentInfo.ProcessedLeg = 'left';

% ExperimentInfo.SubjectID = 'BLSA_5748_04';  % Unequal matrix sizes. 
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_5748_04/310640752_301_T1_3D_80_SLICE_20091106/BLSA_5748_04_310640752_301_T1_3D_80_SLICE_20091106.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_5748_04/310640752_401_WS_T1_3D_80_SLICE_20091106/BLSA_5748_04_310640752_401_WS_T1_3D_80_SLICE_20091106.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_5748_04/310640752_501_FS_T1_3D_80_SLICE_20091106/BLSA_5748_04_310640752_501_FS_T1_3D_80_SLICE_20091106.hdr';
% ExperimentInfo.ProcessedLeg = 'left';

% ExperimentInfo.SubjectID = 'BLSA_6063_03';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%    'BLSA_6063_03/311594904_401_T1_3D_80_SLICE_20091117/BLSA_6063_03_311594904_401_T1_3D_80_SLICE_20091117.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%    'BLSA_6063_03/311594904_501_WS_T1_3D_80_SLICE_20091117/BLSA_6063_03_311594904_501_WS_T1_3D_80_SLICE_20091117.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%    'BLSA_6063_03/311594904_601_FS_T1_3D_80_SLICE_20091117/BLSA_6063_03_311594904_601_FS_T1_3D_80_SLICE_20091117.hdr';
% ExperimentInfo.ProcessedLeg = 'left';
% ExperimentInfo.LegSelectionThreshold = 400 / 1500;  %400/1500

% BLSA_7879_01
% ExperimentInfo.SubjectID = 'BLSA_7879_01';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_7879_01/309971509_201_T1_3D_150_SLICE_20091029/BLSA_7879_01_309971509_201_T1_3D_150_SLICE_20091029.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_7879_01/309971509_301_WS_T1_3D_150_SLICE_20091029/BLSA_7879_01_309971509_301_WS_T1_3D_150_SLICE_20091029.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_7879_01/309971509_401_FS_T1_3D_150_SLICE_20091029/BLSA_7879_01_309971509_401_FS_T1_3D_150_SLICE_20091029.hdr';
% ExperimentInfo.ProcessedLeg = 'left';
% % Top hat parameters.
% ExperimentInfo.tophatTransform = 0;
% ExperimentInfo.strelSizeFactorTopHat = 16;
% ExperimentInfo.LegSelectionThreshold = 350 / 1500;  %400/1500
% 
% BLSA_5898_06
ExperimentInfo.SubjectID = 'BLSA_5898_06';
ExperimentInfo.FirstUseableSlice = 30; 
ExperimentInfo.LastUseableSlice = 50;
ExperimentInfo.LRShiftWS = -0.961; %-0.961;
ExperimentInfo.LRShiftFS = 2.272; %2.272;
ExperimentInfo.NonSuppressedFilename = ...
    'BLSA_5898_06/312283400_301_T1_3D_80_SLICE_20091125/BLSA_5898_06_312283400_301_T1_3D_80_SLICE_20091125.hdr';
ExperimentInfo.WaterSuppressedFilename = ...
    'BLSA_5898_06/312283400_401_WS_T1_3D_80_SLICE_20091125/BLSA_5898_06_312283400_401_WS_T1_3D_80_SLICE_20091125.hdr';
ExperimentInfo.FatSuppressedFilename = ...
    'BLSA_5898_06/312283400_501_FS_T1_3D_80_SLICE_20091125/BLSA_5898_06_312283400_501_FS_T1_3D_80_SLICE_20091125.hdr';
ExperimentInfo.ProcessedLeg = 'left';


% BLSA_6108_03
% ExperimentInfo.SubjectID = 'BLSA_6108_03';
% ExperimentInfo.FirstUseableSlice = 47; 
% ExperimentInfo.LastUseableSlice = 48;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
%ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%    'BLSA_6108_03/309258806_201_TI_3D_95_SLICE_20091021/BLSA_6108_03_309258806_201_TI_3D_95_SLICE_20091021.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%    'BLSA_6108_03/309258806_301_WS_T1_3D_95_SLICE_20091021/BLSA_6108_03_309258806_301_WS_T1_3D_95_SLICE_20091021.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%    'BLSA_6108_03/309258806_401_FS_T1_3D_95_SLICE_20091021/BLSA_6108_03_309258806_401_FS_T1_3D_95_SLICE_20091021.hdr';
% ExperimentInfo.ProcessedLeg = 'right';

% % BLSA_1698_06
% ExperimentInfo.SubjectID = 'BLSA_1698_06';
% ExperimentInfo.FirstUseableSlice = 47; 
% ExperimentInfo.LastUseableSlice = 48;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_1698_06/309357122_201_TI_3D_95_SLICE_20091022/BLSA_1698_06_309357122_201_TI_3D_95_SLICE_20091022.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_1698_06/309357122_301_WS_T1_3D_95_SLICE_20091022/BLSA_1698_06_309357122_301_WS_T1_3D_95_SLICE_20091022.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_1698_06/309357122_401_FS_T1_3D_95_SLICE_20091022/BLSA_1698_06_309357122_401_FS_T1_3D_95_SLICE_20091022.hdr';
% ExperimentInfo.ProcessedLeg = 'left';

% BLSA_6030_04
% ExperimentInfo.SubjectID = 'BLSA_6030_04';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_6030_04/313427769_301_T1_3D_80_SLICE_20091208/BLSA_6030_04_313427769_301_T1_3D_80_SLICE_20091208.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_6030_04/313427769_401_WS_T1_3D_80_SLICE_20091208/BLSA_6030_04_313427769_401_WS_T1_3D_80_SLICE_20091208.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_6030_04/313427769_501_FS_T1_3D_80_SLICE_20091208/BLSA_6030_04_313427769_501_FS_T1_3D_80_SLICE_20091208.hdr';
% ExperimentInfo.ProcessedLeg = 'right';
% ExperimentInfo.LegSelectionThreshold = 350 / 1500;  %400/1500
% 
% % BLSA_1632_09
% ExperimentInfo.SubjectID = 'BLSA_1632_09';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_1632_09/313494902_301_T1_3D_80_SLICE_20091209/BLSA_1632_09_313494902_301_T1_3D_80_SLICE_20091209.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_1632_09/313494902_401_WS_T1_3D_80_SLICE_20091209/BLSA_1632_09_313494902_401_WS_T1_3D_80_SLICE_20091209.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_1632_09/313494902_501_FS_T1_3D_80_SLICE_20091209/BLSA_1632_09_313494902_501_FS_T1_3D_80_SLICE_20091209.hdr';
% ExperimentInfo.ProcessedLeg = 'left';
% ExperimentInfo.nScales = 2;
% ExperimentInfo.snakemaskseSize = 1;
% 
% % BLSA_5519_07
% ExperimentInfo.SubjectID = 'BLSA_5519_07';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_5519_07/313513983_301_T1_3D_80_SLICE_20091209/BLSA_5519_07_313513983_301_T1_3D_80_SLICE_20091209.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_5519_07/313513983_401_WS_T1_3D_80_SLICE_20091209/BLSA_5519_07_313513983_401_WS_T1_3D_80_SLICE_20091209.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_5519_07/313513983_501_FS_T1_3D_80_SLICE_20091209/BLSA_5519_07_313513983_501_FS_T1_3D_80_SLICE_20091209.hdr';
% ExperimentInfo.ProcessedLeg = 'right';

% % BLSA_1531_06
% ExperimentInfo.SubjectID = 'BLSA_1531_06';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_1531_06/313578979_401_T1_3D_80_SLICE_20091210/BLSA_1531_06_313578979_401_T1_3D_80_SLICE_20091210.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_1531_06/313578979_501_WS_T1_3D_80_SLICE_20091210/BLSA_1531_06_313578979_501_WS_T1_3D_80_SLICE_20091210.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_1531_06/313578979_601_FS_T1_3D_80_SLICE_20091210/BLSA_1531_06_313578979_601_FS_T1_3D_80_SLICE_20091210.hdr';
% ExperimentInfo.ProcessedLeg = 'right';
% ExperimentInfo.nScales = 1;
% ExperimentInfo.snakemaskseSize = 2;
% ExperimentInfo.removeBone = 1;
% ExperimentInfo.boneremovalThreshold = 10;

% % BLSA_1945_04
% ExperimentInfo.SubjectID = 'BLSA_1945_04';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_1945_04/313598045_301_T1_3D_80_SLICE_20091210/BLSA_1945_04_313598045_301_T1_3D_80_SLICE_20091210.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_1945_04/313598045_401_WS_T1_3D_80_SLICE_20091210/BLSA_1945_04_313598045_401_WS_T1_3D_80_SLICE_20091210.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_1945_04/313598045_501_FS_T1_3D_80_SLICE_20091210/BLSA_1945_04_313598045_501_FS_T1_3D_80_SLICE_20091210.hdr';
% ExperimentInfo.ProcessedLeg = 'right';
% ExperimentInfo.nScales = 2;
% ExperimentInfo.snakemaskseSize = 4;

% % BLSA_6053_06
% ExperimentInfo.SubjectID = 'BLSA_6053_06';
% ExperimentInfo.FirstUseableSlice = 39; 
% ExperimentInfo.LastUseableSlice = 40;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_6053_06/313665975_301_T1_3D_80_SLICE_20091211/BLSA_6053_06_313665975_301_T1_3D_80_SLICE_20091211.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_6053_06/313665975_401_WS_T1_3D_80_SLICE_20091211/BLSA_6053_06_313665975_401_WS_T1_3D_80_SLICE_20091211.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_6053_06/313665975_501_FS_T1_3D_80_SLICE_20091211/BLSA_6053_06_313665975_501_FS_T1_3D_80_SLICE_20091211.hdr';
% ExperimentInfo.ProcessedLeg = 'left';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # To Do (in Ken's dB)
%BLSA_4903_01, pronounced inhomogeneity artifacts.
% ExperimentInfo.SubjectID = 'BLSA_4903_01';
% ExperimentInfo.FirstUseableSlice = 47; 
% ExperimentInfo.LastUseableSlice = 48;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_4903_01/309345284_201_TI_3D_95_SLICE_20091022/BLSA_4903_01_309345284_201_TI_3D_95_SLICE_20091022.hdr';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_4903_01/309345284_301_WS_T1_3D_95_SLICE_20091022/BLSA_4903_01_309345284_301_WS_T1_3D_95_SLICE_20091022.hdr';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_4903_01/309345284_401_FS_T1_3D_95_SLICE_20091022/BLSA_4903_01_309345284_401_FS_T1_3D_95_SLICE_20091022.hdr';
% ExperimentInfo.ProcessedLeg = 'left';
% ExperimentInfo.LegSelectionThreshold = 200 / 1500;  %400/1500
% ExperimentInfo.nScales = 1;
% ExperimentInfo.snakemaskseSize = 2;
% 
% BLSA_6096_03, calf only.
% ExperimentInfo.SubjectID = 'BLSA_6096_03';
% ExperimentInfo.FirstUseableSlice = 47; 
% ExperimentInfo.LastUseableSlice = 48;
% ExperimentInfo.LRShiftWS = -0.961; %-0.961;
% ExperimentInfo.LRShiftFS = 2.272; %2.272;
% ExperimentInfo.NonSuppressedFilename = ...
%     'BLSA_6096_03\';
% ExperimentInfo.WaterSuppressedFilename = ...
%     'BLSA_6096_03\';
% ExperimentInfo.FatSuppressedFilename = ...
%     'BLSA_6096_03\';
% ExperimentInfo.ProcessedLeg = 'left';


ExperimentInfo.filenames = {ExperimentInfo.NonSuppressedFilename, ...
    ExperimentInfo.FatSuppressedFilename, ...
    ExperimentInfo.WaterSuppressedFilename};

