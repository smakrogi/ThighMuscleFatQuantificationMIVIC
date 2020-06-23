function ExperimentInfo = SetDataPathandFilenamesThighNew(InfoFiles, Path)
% SetDataPathandFilenamesThighNew

% Debug mode.
ExperimentInfo.debug = 0;
% Command line mode.
ExperimentInfo.commandlineMode = 0;
% Slice visualization step.
ExperimentInfo.sliceStep = 2;
% Algorithm parameters for:
% Rescale factor before analysis in feature space.
ExperimentInfo.rescaleFactor = 1.0;
% Number of clusters.
ExperimentInfo.nClusters = 2;
% Shading correction.
ExperimentInfo.tophatTransform = 0;
ExperimentInfo.strelSizeFactorTopHat = 16; %8;
% Structural element size for definition of initial mask.
ExperimentInfo.snakemaskseSize = 4; % 1,4,10
% Gaussian smoothing sigma before snake segmentation.
ExperimentInfo.snakeGaussianSigma = 0.75;
% Remove the bone.
ExperimentInfo.removeBone = 1;
% ExperimentInfo.boneremovalThreshold = 100;
ExperimentInfo.strelSizeFactorBoneRemoval = 256; % 128;
ExperimentInfo.boneAreaThreshold = 0.028;
ExperimentInfo.standarddevCorticalBoneRemoval = 1.0; % 1.5
ExperimentInfo.rgIterationsCorticalBoneRemoval = 2;
ExperimentInfo.seedingIterationsCorticalBoneRemoval = 12;
% Clustering.
ExperimentInfo.clusterNonZeroValuesOnly = 1;
% Warp FS and WS to NS before clustering.
ExperimentInfo.warpWSandFStoNS = 0;
% Post-clustering bone removal.
% ExperimentInfo.areaThreshold = 20;
% Distance metric learning.
ExperimentInfo.DistanceLearningMethods = ...
    {'', 'PCA', 'KPCA', 'PPCA', 'ISOMAP', 'Laplacian', 'LLE', 'LTSA', 'MDS'};
ExperimentInfo.DistanceLearningMethod = ...
    ExperimentInfo.DistanceLearningMethods{1};
ExperimentInfo.OutputDimensionality = 2;
% Clustering algorithms.
ExperimentInfo.ClusteringMethods = ...
    {'K-MEANS', 'FCM', 'emgmmbaycls', 'ncutclus', 'MeanShift', 'emgmmquad'};
ExperimentInfo.ClusteringMethod = ExperimentInfo.ClusteringMethods{3};
ExperimentInfo.MeanShiftBandwidth = 200;  % 200
ExperimentInfo.GMMnDistributions = 2;
ExperimentInfo.CovarianceMatrixTypes = ...
{'Full', 'Diag', 'Spherical'};
ExperimentInfo.CovarianceMatrixType = ...
    ExperimentInfo.CovarianceMatrixTypes{1};
% SAT extraction.
ExperimentInfo.SATExtractionAlgorithms = {'Snake', 'GAC'};
ExperimentInfo.SATExtractionAlgorithm = ExperimentInfo.SATExtractionAlgorithms{1};
% Snakes parameters.
ExperimentInfo.snakeExternalForces = { 'gradient', 'gvf', 'dtransform'};
ExperimentInfo.snakeExternalForce =  ExperimentInfo.snakeExternalForces{1};
ExperimentInfo.gvfInputs = {'Canny', 'Parzen', 'Gradient'};
ExperimentInfo.gvfInput = ExperimentInfo.gvfInputs{1};
ExperimentInfo.SupplementaryDataFilenameSuffix = '_SupplementaryData.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    % Path to datasets.
    ExperimentInfo.dataPath = 'Z:\workspace\Preprocessing\Experiments_20_slices.20100219\';
elseif isunix
    % Path to datasets.
    ExperimentInfo.dataPath = '/home/makrogianniss/workspace/Preprocessing/Experiments_20_slices.20100219/';
end

ExperimentInfo.FirstUseableSlice = 8; %2; %9; %20100219_Experiments
ExperimentInfo.LastUseableSlice = 12; %15; %11;

% ExperimentInfo.FirstUseableSlice = 2;   %20100224_Experiments
% ExperimentInfo.LastUseableSlice = 4;

ExperimentInfo.LRShiftWS = 0.0; % -0.961;
ExperimentInfo.LRShiftFS = 0.0; % 2.272;

% Parse subject and image lists.
ExperimentInfo.ProcessedLeg = 'left';
SubjectID_List_Prefix = 'SubjectIDs_List_Nii_Preprocessed';
NSImage_List_Prefix = 'T13D_List_Nii_Preprocessed';
WSImage_List_Prefix = 'WST13D_List_Nii_Preprocessed';
FSImage_List_Prefix = 'FST13D_List_Nii_Preprocessed';

foundFiles = 0;
for i=1:length(InfoFiles)
       % Compare filenames with text patterns.
       % SubjectIDs_List_Nii_Preprocessed_
       if ~isempty(strfind(InfoFiles{i}, SubjectID_List_Prefix))
            % Read-in the subject list.
            [ExperimentInfo.SubjectID, ExperimentInfo.dataPath] = ...
                ReadThighInfoFile([Path, InfoFiles{i}]);   
            % C = textscan(fid, '%s %s %s %s %s', 'delimiter', '/');
            foundFiles = foundFiles + 1;
       % WST13D_List_Nii_Preprocessed_
       elseif ~isempty(strfind(InfoFiles{i}, WSImage_List_Prefix))
            ExperimentInfo.WaterSuppressedFilename = ...
                ReadThighInfoFile([Path, InfoFiles{i}]);   
            foundFiles = foundFiles + 1;
       % FST13D_List_Nii_Preprocessed_
       elseif ~isempty(strfind(InfoFiles{i}, FSImage_List_Prefix))
            ExperimentInfo.FatSuppressedFilename = ...
                ReadThighInfoFile([Path, InfoFiles{i}]);      
            foundFiles = foundFiles + 1;
       % T13D_List_Nii_Preprocessed_
       elseif ~isempty(strfind(InfoFiles{i}, NSImage_List_Prefix))
            ExperimentInfo.NonSuppressedFilename = ...
                ReadThighInfoFile([Path, InfoFiles{i}]);       
            foundFiles = foundFiles + 1;
       end
end

if foundFiles ~= 4
    error('Runtime error: Could not read all info files.')
end

% Initialize these two parameters.
ExperimentInfo.dataPath = [ExperimentInfo.dataPath filesep];
ExperimentInfo.nSubjects = length ( ExperimentInfo.SubjectID) ;
ExperimentInfo.nScales = cell(ExperimentInfo.nSubjects, 1);
ExperimentInfo.LegSelectionThreshold = cell(ExperimentInfo.nSubjects, 1);
for(i=1:ExperimentInfo.nSubjects)
    % Number of scales in deformable segmentation.
    ExperimentInfo.nScales{i} = 3;
    % Leg selection.
    ExperimentInfo.LegSelectionThreshold{i} = 0.3;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [StringCell, firstLine] = ReadThighInfoFile(filename)
    fid = fopen( filename );
    line = fgetl(fid);
    i = 0;
    while ischar(line)
        if i==0, firstLine = line;
        else StringCell{i} = line;
        end
        line = fgetl(fid);
        i=i+1;
    end
    fclose(fid);
end