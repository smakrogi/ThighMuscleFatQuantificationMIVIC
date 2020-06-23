function [stats, TissueLabelVolume2, infoString] = CalculateThighMuscleandAdiposeVolumes...
 (subjectNumber, ExperimentInfo, tissueLabels, DistributionModels, varargin)
% syntax: [stats, TissueLabelVolume2, infoString] = CalculateThighMuscleandAdiposeVolumes...
%   (subjectNumber, ExperimentInfo, tissueLabels, DistributionModels,
%   varargin);
% 3T NIA/NIH, by S.K. Makrogiannis. 

% threshold = 50;
% strelSizeFactor = 128;

% Read volumes with intermediate results.
TissueLabelVolume = load_untouch_nii([varargin{1}]);
TissueLabelVolume = uint16(TissueLabelVolume.img);
OneLegVolume = load_untouch_nii([varargin{2}]);
OneLegVolume = uint16(OneLegVolume.img);
SupplementaryData = LoadSupplementaryData(subjectNumber, ExperimentInfo);
BoundingBox = SupplementaryData.BoundingBox;    
InternalSATSurfaceMask = load_untouch_nii([varargin{3}]);
InternalSATSurfaceMask = uint16(InternalSATSurfaceMask.img);
BoneLabelVolume = load_untouch_nii([varargin{4}]);
BoneLabelVolume = uint16(BoneLabelVolume.img);
% Create NoBoneMask.
NoBoneMask = uint16(OneLegVolume>0) - uint16(BoneLabelVolume > 0);


InternalSATSurfaceMask = uint16(InternalSATSurfaceMask);

% Read voxel size.
[nonsuppressedHeader, filetype] = ...
    load_nii_hdr([ExperimentInfo.dataPath ExperimentInfo.NonSuppressedFilename{subjectNumber}]);
nonsuppressedVoxelVolume = prod(nonsuppressedHeader.dime.pixdim(2:4));
nSlices = ExperimentInfo.LastUseableSlice - ...
    ExperimentInfo.FirstUseableSlice + ...
    1; 
voxelHeight = nonsuppressedHeader.dime.pixdim(4);
voxelSize = nonsuppressedHeader.dime.pixdim(2:4);
description = nonsuppressedHeader.hist.descrip;
origin = [nonsuppressedHeader.hist.qoffset_x nonsuppressedHeader.hist.qoffset_y nonsuppressedHeader.hist.qoffset_z];

% Subtract muscle and bone from one leg to show SAT.
InternalSATSurfaceMask = InternalSATSurfaceMask == 1;
SATVolume = OneLegVolume .*uint16(1-InternalSATSurfaceMask);

% Produce tissue masks inside subcute region.
MuscleMask = uint16(TissueLabelVolume == tissueLabels.muscle);
InterMFMask = uint16(TissueLabelVolume == tissueLabels.fat);
% List all object labels.
labels = unique(TissueLabelVolume(:));
% Remove zero label.
[temp1, temp2, labels] = find( labels );
clear temp1; clear temp2;
[temp1, temp2, found] = find( labels ~= tissueLabels.muscle & labels ~= tissueLabels.fat);
if(~isempty(found))
    otherLabels = zeros(size(found));
    for i=1:length(temp1)
        otherLabels(i) = labels(temp1(i), temp2(i));
    end
    IntraMFMask = uint16(TissueLabelVolume == otherLabels);
else
    IntraMFMask = uint16(zeros(size(TissueLabelVolume)));
end
clear temp1; clear temp2;
BoneMask = uint16(InternalSATSurfaceMask) .* uint16(1 - NoBoneMask);


% Compute basic tissue variables.
stats.SATVoxels = sum( SATVolume(:)>0 );
stats.SATVolume = stats.SATVoxels * nonsuppressedVoxelVolume;
stats.AverageSATArea = stats.SATVolume / (voxelHeight * nSlices);
stats.MuscleVoxels = sum( MuscleMask(:) );
stats.MuscleVolume = stats.MuscleVoxels * nonsuppressedVoxelVolume;
stats.AverageMuscleArea = stats.MuscleVolume / (voxelHeight * nSlices);
stats.InterMFVoxels = sum( InterMFMask(:) );
stats.InterMFVolume = stats.InterMFVoxels * nonsuppressedVoxelVolume;
stats.AverageInterMFArea = stats.InterMFVolume / (voxelHeight * nSlices);
stats.IntraMFVoxels = sum( IntraMFMask(:) );
stats.IntraMFVolume = stats.IntraMFVoxels * nonsuppressedVoxelVolume;
stats.AverageIntraMFArea = stats.IntraMFVolume / (voxelHeight * nSlices);
stats.BoneVoxels = sum( BoneMask(:) );
stats.BoneVolume = stats.BoneVoxels * nonsuppressedVoxelVolume;
stats.AverageBoneArea = stats.BoneVolume / (voxelHeight * nSlices);


% Parse variables to string and display to console.
infoString = [];
infoString = [infoString sprintf('\nSubject ID, %s\n', ExperimentInfo.SubjectID{subjectNumber})];
infoString = [infoString sprintf('Voxel size, [%.3f %.3f %.3f] (mm), voxel volume, %.3f (mm^3)\n', ...
    nonsuppressedHeader.dime.pixdim(2:4), nonsuppressedVoxelVolume)];
infoString = [infoString sprintf('First used slice, %d, Last used slice, %d \n', ...
    ExperimentInfo.FirstUseableSlice, ExperimentInfo.LastUseableSlice)];
infoString = [infoString sprintf('Using %s leg... \n', ExperimentInfo.ProcessedLeg)];
% infoString = [infoString sprintf('Total number of Muscle voxels, %d\n', stats.MuscleVoxels)];
% infoString = [infoString sprintf('Muscle volume, %.3f (mm^3)\n', stats.MuscleVolume)];
infoString = [infoString sprintf('Muscle area ave., %.3f (mm^2)\n', stats.AverageMuscleArea)];
% infoString = [infoString sprintf('Total number of SAT voxels, %d\n', stats.SATVoxels)];
% infoString = [infoString sprintf('SAT volume, %.3f (mm^3)\n', stats.SATVolume)];
infoString = [infoString sprintf('SAT area ave., %.3f (mm^2)\n', stats.AverageSATArea)];
% infoString = [infoString sprintf('Total number of InterMF voxels, %d\n', stats.InterMFVoxels)];
% infoString = [infoString sprintf('InterMF volume, %.3f\n', stats.InterMFVolume)];
infoString = [infoString sprintf('InterMF area ave., %.3f (mm^2)\n', stats.AverageInterMFArea)];
% infoString = [infoString sprintf('Total number of IntraMF voxels, %d\n', stats.IntraMFVoxels)];
% infoString = [infoString sprintf('IntraMF volume, %.3f\n', stats.IntraMFVolume)];
infoString = [infoString sprintf('IntraMF area ave., %.3f (mm^2)\n', stats.AverageIntraMFArea)];
% infoString = [infoString sprintf('Total number of Bone voxels, %d\n', stats.BoneVoxels)];
% infoString = [infoString sprintf('Bone volume, %.3f\n', stats.BoneVolume)];
infoString = [infoString sprintf('Bone area ave., %.3f (mm^2)\n', stats.AverageBoneArea)];
fprintf( infoString );


% Store voxel counts in csv format.
imagefileprefix = [ExperimentInfo.SubjectID{subjectNumber}, '_', ExperimentInfo.DistanceLearningMethod, ...
        '_', ExperimentInfo.ClusteringMethod];

% fid=fopen([imagefileprefix '_ThighQuantification.csv'], 'wb', 'l');
% fprintf(fid, infoString2);
% fclose(fid);


% Create slice montages.
if ( ExperimentInfo.commandlineMode )
    AxialSliceMontage(OneLegVolume, ExperimentInfo.sliceStep);
    saveas(gcf, [ExperimentInfo.SubjectID{subjectNumber} '_oneleg_montage'], 'png')
    
    AxialSliceMontage(SATVolume, ExperimentInfo.sliceStep);
    saveas(gcf, [imagefileprefix '_sat_segmentation_montage'], 'png')
    
    AxialSliceMontage(OneLegVolume.*uint16(MuscleMask), ExperimentInfo.sliceStep);
    saveas(gcf, [imagefileprefix, '_muscle_segmentation_montage'], 'png')
    
    AxialSliceMontage(OneLegVolume.*uint16(InterMFMask), ExperimentInfo.sliceStep);
    saveas(gcf, [imagefileprefix, '_InterMF_segmentation_montage'], 'png')

    AxialSliceMontage(OneLegVolume.*uint16(IntraMFMask), ExperimentInfo.sliceStep);
    saveas(gcf, [imagefileprefix, '_IntraMF_segmentation_montage'], 'png')
end


% Create second label volume and save in nifti format.
InterMF_Label = 2;
IntraMF_Label = 3;
Muscle_Label  = 4;
TissueLabelVolume2 = uint16(SATVolume > 0) + ...
    InterMF_Label * InterMFMask + ...
    IntraMF_Label * IntraMFMask + ...
    Muscle_Label * MuscleMask;
% TissueLabelVolume2 = uint16(TissueLabelVolume2) + uint16(TissueLabelVolume);
TissueLabelVolume2 = uint16(TissueLabelVolume2) ;

if ( ExperimentInfo.commandlineMode )
    AxialSliceMontage(uint16(TissueLabelVolume2), ExperimentInfo.sliceStep);
    saveas(gcf, [imagefileprefix, '_label_segmentation_montage'], 'png')
end
% TissueLabelVolume2 = TransposeVolume(TissueLabelVolume2);
TissueLabelVolume2Nii = make_nii(TissueLabelVolume2, voxelSize, origin, 4, ...
    description);
save_nii(TissueLabelVolume2Nii, [imagefileprefix '_TissueLabelVolume2.hdr']);
clear TissueLabelVolume2Nii;


% Calculate stats and display histograms of each tissue type.
I_SAT = find(SATVolume(:)>0);
I_Muscle = find(MuscleMask(:));
I_InterMF = find(InterMFMask(:));
I_IntraMF = find(IntraMFMask(:));
I_Bone =  find(BoneMask(:));
clear *Volume;
stats = NonSuppressedIntensityStats(subjectNumber, ExperimentInfo, stats, ...
    I_SAT, I_Muscle, I_InterMF, I_IntraMF, I_Bone, BoundingBox, imagefileprefix);
stats = FatSuppressedIntensityStats(subjectNumber, nonsuppressedHeader, ...
    ExperimentInfo, stats, I_SAT, I_Muscle, I_InterMF, I_IntraMF, I_Bone, BoundingBox, imagefileprefix);
stats = WaterSuppressedIntensityStats(subjectNumber, nonsuppressedHeader, ...
    ExperimentInfo, stats, I_SAT, I_Muscle, I_InterMF, I_IntraMF, I_Bone, BoundingBox, imagefileprefix);


% Extract Gaussian model parameters.
if ~isempty(DistributionModels)
    stats = ExtractCovarianceMatrixElements(DistributionModels, stats);
end


% Add intensity statistics to infostring.
infoString = [infoString sprintf('Muscle NS intensity ave, %.3f\n', stats.nonsuppressed_muscle_mean)];
infoString = [infoString sprintf('SAT NS intensity ave, %.3f\n', stats.nonsuppressed_sat_mean)];
infoString = [infoString sprintf('InterMF NS intensity ave, %.3f\n', stats.nonsuppressed_intermf_mean)];
infoString = [infoString sprintf('IntraMF NS intensity ave, %.3f\n', stats.nonsuppressed_intramf_mean)];
infoString = [infoString sprintf('Bone NS intensity ave, %.3f\n', stats.nonsuppressed_bone_mean)];
infoString = [infoString sprintf('Muscle F intensity ave, %.3f\n', stats.fatsuppressed_muscle_mean)];
infoString = [infoString sprintf('SAT FS intensity ave, %.3f\n', stats.fatsuppressed_sat_mean)];
infoString = [infoString sprintf('InterMF FS intensity ave, %.3f\n', stats.fatsuppressed_intermf_mean)];
infoString = [infoString sprintf('IntraMF FS intensity ave, %.3f\n', stats.fatsuppressed_intramf_mean)];
infoString = [infoString sprintf('Bone FS intensity ave, %.3f\n', stats.fatsuppressed_bone_mean)];
infoString = [infoString sprintf('Muscle WS intensity ave, %.3f\n', stats.watersuppressed_muscle_mean)];
infoString = [infoString sprintf('SAT WS intensity ave, %.3f\n', stats.watersuppressed_sat_mean)];
infoString = [infoString sprintf('InterMF WS intensity ave, %.3f\n', stats.watersuppressed_intermf_mean)];
infoString = [infoString sprintf('IntraMF WS intensity ave, %.3f\n', stats.watersuppressed_intramf_mean)];
infoString = [infoString sprintf('Bone WS intensity ave, %.3f\n', stats.watersuppressed_bone_mean)];

% Write all statistics to csv text file.
WriteStatisticsToFile(ExperimentInfo, subjectNumber, stats);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newstats = NonSuppressedIntensityStats(subjectNumber, ExperimentInfo, stats, ...
    I_SAT, I_Muscle, I_InterMF, I_IntraMF, I_Bone, BoundingBox, imagefileprefix)

% Read the non saturated volume.
% % infoString = sprintf('Non suppressed volume: %s \n',[ExperimentInfo.dataPath ExperimentInfo.filename]);
% % fprintf( infoString );
% NonSuppressed = analyze75read([ExperimentInfo.dataPath ExperimentInfo.filename]);
NonSuppressed = load_untouch_nii([ExperimentInfo.dataPath ExperimentInfo.NonSuppressedFilename{subjectNumber}]);
NonSuppressed = NonSuppressed.img(:, :, ExperimentInfo.FirstUseableSlice:ExperimentInfo.LastUseableSlice);
% NonSuppressed = TransposeVolume(NonSuppressed);
NonSuppressed = NonSuppressed(BoundingBox.xmin:BoundingBox.xmax, ...
    BoundingBox.ymin:BoundingBox.ymax, :);

NonSuppressed = NonSuppressed(:);
MuscleNonSuppressed = NonSuppressed(I_Muscle);
stats.nonsuppressed_muscle_mean = mean(MuscleNonSuppressed);
SAT_NonSuppressed = NonSuppressed(I_SAT);
stats.nonsuppressed_sat_mean = mean(SAT_NonSuppressed);
InterMF_NonSuppressed = NonSuppressed(I_InterMF);
stats.nonsuppressed_intermf_mean = mean(InterMF_NonSuppressed);
IntraMF_NonSuppressed = NonSuppressed(I_IntraMF);
stats.nonsuppressed_intramf_mean = mean(IntraMF_NonSuppressed);
Bone_NonSuppressed = NonSuppressed(I_Bone);
stats.nonsuppressed_bone_mean = mean(Bone_NonSuppressed);


if ( ExperimentInfo.commandlineMode )
    figure, subplot(141), hist(single(MuscleNonSuppressed),100), title('NS-Muscle'), grid on, axis square;
    subplot(142), hist(single(SAT_NonSuppressed),100), title('NS-SAT'), grid on, axis square;
    subplot(143), hist(single(InterMF_NonSuppressed),100), title('NS-InterMF'), grid on, axis square;
    subplot(144), hist(single(IntraMF_NonSuppressed),100), title('NS-IntraMF'), grid on, axis square;    
    saveas(gcf, [imagefileprefix, '_ns_tissue_histograms.png']);
end

infoString = sprintf('Non-suppressed intensity averages: muscle, %.3f, SAT, %.3f, InterMF, %.3f, IntraMF, %.3f, Bone, %.3f\n', ...
    stats.nonsuppressed_muscle_mean, stats.nonsuppressed_sat_mean, stats.nonsuppressed_intermf_mean, ...
    stats.nonsuppressed_intramf_mean, stats.nonsuppressed_bone_mean);
fprintf( infoString );

newstats = stats;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newstats = FatSuppressedIntensityStats(subjectNumber,nonsuppressedHeader, ExperimentInfo, stats, ...
    I_SAT, I_Muscle, I_InterMF, I_IntraMF, I_Bone, BoundingBox, imagefileprefix)

% Read the fat saturated volume.
FatSuppressed = load_untouch_nii([ExperimentInfo.dataPath ExperimentInfo.FatSuppressedFilename{subjectNumber}]);
FatSuppressed = FatSuppressed.img(:, :, ExperimentInfo.FirstUseableSlice:ExperimentInfo.LastUseableSlice);
% FatSuppressed = TransposeVolume(FatSuppressed);
referenceMatrixSize =  nonsuppressedHeader.dime.dim(2:3);
resizeFactor = referenceMatrixSize(1) / size(FatSuppressed,1);
if( resizeFactor ~= 1 )
    fprintf('Need to resize fat volume to match non-suppressed.\n');
    
    FatSuppressed = ResizeVolumeMatrixSize(FatSuppressed, referenceMatrixSize);
end
FatSuppressed = ShiftVolumeLR(FatSuppressed, ExperimentInfo.LRShiftFS * resizeFactor);
FatSuppressed = FatSuppressed(BoundingBox.xmin:BoundingBox.xmax, ...
    BoundingBox.ymin:BoundingBox.ymax, :);

FatSuppressed = FatSuppressed(:);
MuscleFatSuppressed = FatSuppressed(I_Muscle);
stats.fatsuppressed_muscle_mean = mean(MuscleFatSuppressed);
SAT_FatSuppressed = FatSuppressed(I_SAT);
stats.fatsuppressed_sat_mean = mean(SAT_FatSuppressed);
InterMF_FatSuppressed = FatSuppressed(I_InterMF);
stats.fatsuppressed_intermf_mean = mean(InterMF_FatSuppressed);
IntraMF_FatSuppressed = FatSuppressed(I_IntraMF);
stats.fatsuppressed_intramf_mean = mean(IntraMF_FatSuppressed);
Bone_FatSuppressed = FatSuppressed(I_Bone);
stats.fatsuppressed_bone_mean = mean(Bone_FatSuppressed);

if ( ExperimentInfo.commandlineMode )
    figure, subplot(141), hist(single(MuscleFatSuppressed),100), title('FS-Muscle'), grid on, axis square;
    subplot(142), hist(single(SAT_FatSuppressed),100), title('FS-SAT'), grid on, axis square;
    subplot(143), hist(single(InterMF_FatSuppressed),100), title('FS-InterMF'), grid on, axis square;
    subplot(144), hist(single(IntraMF_FatSuppressed),100), title('FS-IntraMF'), grid on, axis square;
    saveas(gcf, [imagefileprefix, '_fs_tissue_histograms.png']);
end

infoString = sprintf('Fat-suppressed intensity averages: muscle, %.3f, SAT, %.3f, InterMF, %.3f, IntraMF, %.3f, Bone, %.3f\n', ...
    stats.fatsuppressed_muscle_mean, stats.fatsuppressed_sat_mean, stats.fatsuppressed_intermf_mean, ...
    stats.fatsuppressed_intramf_mean, stats.fatsuppressed_bone_mean);
fprintf( infoString );

newstats = stats;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newstats = WaterSuppressedIntensityStats(subjectNumber, nonsuppressedHeader, ExperimentInfo, stats, ...
    I_SAT, I_Muscle, I_InterMF, I_IntraMF, I_Bone, BoundingBox, imagefileprefix)

% Read the water saturated volume.

WaterSuppressed = load_untouch_nii([ExperimentInfo.dataPath ExperimentInfo.WaterSuppressedFilename{subjectNumber}]);
WaterSuppressed = WaterSuppressed.img(:, :, ExperimentInfo.FirstUseableSlice:ExperimentInfo.LastUseableSlice);
% WaterSuppressed = TransposeVolume(WaterSuppressed);
referenceMatrixSize =  nonsuppressedHeader.dime.dim(2:3);
resizeFactor = referenceMatrixSize(1) / size(WaterSuppressed,1);
if( resizeFactor ~= 1 )
    fprintf('Need to resize water volume to match non-suppressed.\n');
    WaterSuppressed = ResizeVolumeMatrixSize(WaterSuppressed, referenceMatrixSize);
end
WaterSuppressed = ShiftVolumeLR(WaterSuppressed, ExperimentInfo.LRShiftWS * resizeFactor);
WaterSuppressed = WaterSuppressed(BoundingBox.xmin:BoundingBox.xmax, ...
    BoundingBox.ymin:BoundingBox.ymax, :);

WaterSuppressed = WaterSuppressed(:);
MuscleWaterSuppressed = WaterSuppressed(I_Muscle);
stats.watersuppressed_muscle_mean = mean(MuscleWaterSuppressed);
SAT_WaterSuppressed = WaterSuppressed(I_SAT);
stats.watersuppressed_sat_mean = mean(SAT_WaterSuppressed);
InterMF_WaterSuppressed = WaterSuppressed(I_InterMF);
stats.watersuppressed_intermf_mean = mean(InterMF_WaterSuppressed);
IntraMF_WaterSuppressed = WaterSuppressed(I_IntraMF);
stats.watersuppressed_intramf_mean = mean(IntraMF_WaterSuppressed);
Bone_WaterSuppressed = WaterSuppressed(I_Bone);
stats.watersuppressed_bone_mean = mean(Bone_WaterSuppressed);

if ( ExperimentInfo.commandlineMode )
    figure, subplot(141), hist(single(MuscleWaterSuppressed),100), title('WS-Muscle'), grid on, axis square;
    subplot(142), hist(single(SAT_WaterSuppressed),100), title('WS-SAT'), grid on, axis square;
    subplot(143), hist(single(InterMF_WaterSuppressed),100), title('WS-InterMF'), grid on, axis square;
    subplot(144), hist(single(IntraMF_WaterSuppressed),100), title('WS-IntraMF'), grid on, axis square;    
    saveas(gcf, [imagefileprefix, '_ws_tissue_histograms.png']);
end

infoString = sprintf('Water-suppressed intensity averages: muscle, %.3f, SAT, %.3f, InterMF, %.3f, IntraMF, %.3f, Bone, %.3f\n', ...
    stats.watersuppressed_muscle_mean, stats.watersuppressed_sat_mean, stats.watersuppressed_intermf_mean, ...
    stats.watersuppressed_intramf_mean, stats.watersuppressed_bone_mean);
fprintf( infoString );

newstats = stats;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newstats = ExtractCovarianceMatrixElements(DistributionModels, stats)

% Get number of groups or tissues.
nComponents = length(DistributionModels.Pclass);

% Extract matrix elements.
if nComponents == 2 % has no intramf model.
    % Muscle.
    stats.muscle.meanFS = DistributionModels.Pclass{1}.Mean(1);
    stats.muscle.meanWS = DistributionModels.Pclass{1}.Mean(2);
    stats.muscle.covarFSFS = DistributionModels.Pclass{1}.Cov(1,1);
    stats.muscle.covarFSWS = DistributionModels.Pclass{1}.Cov(1,2);
    stats.muscle.covarWSWS = DistributionModels.Pclass{1}.Cov(2,2);
    % InterMF.
    stats.intermf.meanFS = DistributionModels.Pclass{2}.Mean(1);
    stats.intermf.meanWS = DistributionModels.Pclass{2}.Mean(2);
    stats.intermf.covarFSFS = DistributionModels.Pclass{2}.Cov(1,1);
    stats.intermf.covarFSWS = DistributionModels.Pclass{2}.Cov(1,2);
    stats.intermf.covarWSWS = DistributionModels.Pclass{2}.Cov(2,2);
    % IntraMF.
    stats.intramf.meanFS = [];
    stats.intramf.meanWS = [];
    stats.intramf.covarFSFS = [];
    stats.intramf.covarFSWS = [];
    stats.intramf.covarWSWS = [];
elseif nComponents == 3 % has intramf model.
     % Muscle.
    stats.muscle.meanFS = DistributionModels.Pclass{1}.Mean(1);
    stats.muscle.meanWS = DistributionModels.Pclass{1}.Mean(2);
    stats.muscle.covarFSFS = DistributionModels.Pclass{1}.Cov(1,1);
    stats.muscle.covarFSWS = DistributionModels.Pclass{1}.Cov(1,2);
    stats.muscle.covarWSWS = DistributionModels.Pclass{1}.Cov(2,2);
    % InterMF.
    stats.intermf.meanFS = DistributionModels.Pclass{2}.Mean(1);
    stats.intermf.meanWS = DistributionModels.Pclass{2}.Mean(2);
    stats.intermf.covarFSFS = DistributionModels.Pclass{2}.Cov(1,1);
    stats.intermf.covarFSWS = DistributionModels.Pclass{2}.Cov(1,2);
    stats.intermf.covarWSWS = DistributionModels.Pclass{2}.Cov(2,2);
    % IntraMF.
    stats.intramf.meanFS = DistributionModels.Pclass{3}.Mean(1);
    stats.intramf.meanWS = DistributionModels.Pclass{3}.Mean(2);
    stats.intramf.covarFSFS = DistributionModels.Pclass{3}.Cov(1,1);
    stats.intramf.covarFSWS = DistributionModels.Pclass{3}.Cov(1,2);
    stats.intramf.covarWSWS = DistributionModels.Pclass{3}.Cov(2,2);
else
     % Muscle.
    stats.muscle.meanFS = [];
    stats.muscle.meanWS = [];
    stats.muscle.covarFSFS = [];
    stats.muscle.covarFSWS = [];
    stats.muscle.covarWSWS = [];
    % InterMF.
    stats.intermf.meanFS = [];
    stats.intermf.meanWS = [];
    stats.intermf.covarFSFS = [];
    stats.intermf.covarFSWS = [];
    stats.intermf.covarWSWS = [];
    % IntraMF.
    stats.intramf.meanFS = [];
    stats.intramf.meanWS = [];
    stats.intramf.covarFSFS = [];
    stats.intramf.covarFSWS = [];
    stats.intramf.covarWSWS = [];
end

newstats = stats;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutputVolume = ResizeVolumeMatrixSize(InputVolume, matrixSize)
% Resize matrix size before creating the feature vectors.
% matrixSize = [size(ReferenceVolume,1) size(ReferenceVolume,2)];
% OutputVolume = ReferenceVolume;
for i=1:size(InputVolume, 3)
    OutputVolume(:,:,i) = imresize(InputVolume(:,:,i), matrixSize);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutputVolume = ShiftVolumeLR(InputVolume, shiftVector)

% Generate deformation fields with meshgrid and add the shiftVector.
volumeSize = size(InputVolume);
[Y, X] = meshgrid(1:volumeSize(1), 1:volumeSize(2));
Y = Y + shiftVector(1);

% Apply interpolation to each slice.
OutputVolume = zeros(size(InputVolume));
for i=1:volumeSize(3)
    OutputVolume(:,:,i) = interp2(double(InputVolume(:,:,i)), Y, X);
end
OutputVolume = uint16(OutputVolume);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteStatisticsToFile(ExperimentInfo, subjectNumber, stats)

% Store stats in csv format.
infoString2 = [];
infoString2 = [infoString2 sprintf('Subject ID, First used slice, Last used slice, Leg, ')];
infoString2 = [infoString2 sprintf('Muscle area ave., SAT area ave., InterMF area ave., IntraMF area ave., Bone area ave., ')];
infoString2 = [infoString2 sprintf('Muscle NS intensity ave., SAT NS intensity ave., InterMF NS intensity ave., IntraMF NS intensity ave.,Bone NS intensity ave., ')];
infoString2 = [infoString2 sprintf('Muscle FS intensity ave., SAT FS intensity ave., InterMF FS intensity ave., IntraMF FS intensity ave., Bone FS intensity ave., ')];
infoString2 = [infoString2 sprintf('Muscle WS intensity ave., SAT WS intensity ave., InterMF WS intensity ave., IntraMF WS intensity ave., Bone WS intensity ave., ')];
infoString2 = [infoString2 sprintf('Muscle FS EM  ave., Muscle WS EM ave., Muscle FSFS EM covar, Muscle FSWS EM covar, Muscle WSWS EM covar, ')];
infoString2 = [infoString2 sprintf('InterMF FS EM ave., InterMF WS EM ave., InterMF FSFS EM covar, InterMF FSWS EM covar, InterMF WSWS EM covar, ')];
infoString2 = [infoString2 sprintf('IntraMF FS EM ave., IntraMF WS EM ave., IntraMF FSFS EM covar, IntraMF FSWS EM covar, IntraMF WSWS EM covar\n')];
infoString2 = [infoString2 sprintf( '%s, %d, %d, %s, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n', ...
    ExperimentInfo.SubjectID{subjectNumber}, ExperimentInfo.FirstUseableSlice, ExperimentInfo.LastUseableSlice, ExperimentInfo.ProcessedLeg,  ...
    stats.AverageMuscleArea, stats.AverageSATArea, stats.AverageInterMFArea, stats.AverageIntraMFArea, stats.AverageBoneArea, ...
    stats.nonsuppressed_muscle_mean, stats.nonsuppressed_sat_mean, stats.nonsuppressed_intermf_mean, stats.nonsuppressed_intramf_mean, stats.nonsuppressed_bone_mean, ...
    stats.fatsuppressed_muscle_mean, stats.fatsuppressed_sat_mean, stats.fatsuppressed_intermf_mean, stats.fatsuppressed_intramf_mean, stats.fatsuppressed_bone_mean, ...
    stats.watersuppressed_muscle_mean, stats.watersuppressed_sat_mean, stats.watersuppressed_intermf_mean, stats.watersuppressed_intramf_mean, stats.watersuppressed_bone_mean, ...
    stats.muscle.meanFS, stats.muscle.meanWS, stats.muscle.covarFSFS, stats.muscle.covarFSWS, stats.muscle.covarWSWS, ...
    stats.intermf.meanFS, stats.intermf.meanWS, stats.intermf.covarFSFS, stats.intermf.covarFSWS, stats.intermf.covarWSWS, ...
    stats.intramf.meanFS, stats.intramf.meanWS, stats.intramf.covarFSFS, stats.intramf.covarFSWS, stats.intramf.covarWSWS)];

% Store in csv format.
imagefileprefix = [ExperimentInfo.SubjectID{subjectNumber}, '_', ExperimentInfo.DistanceLearningMethod, ...
        '_', ExperimentInfo.ClusteringMethod];

fid=fopen([imagefileprefix '_ThighQuantification.csv'], 'wb', 'l');
fprintf(fid, infoString2);
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SupplementaryData = LoadSupplementaryData(subjectNumber, ExperimentInfo)

    imagefileprefix = [ExperimentInfo.SubjectID{subjectNumber}, '_', ...
    ExperimentInfo.DistanceLearningMethod, '_', ...
    ExperimentInfo.ClusteringMethod];
    SupplementaryData = [];
    
    % Check if a supplementary data file already exists.
    fid = fopen([imagefileprefix,  ExperimentInfo.SupplementaryDataFilenameSuffix]);
    if fid~=-1
        fclose( fid );
        load([imagefileprefix,  ExperimentInfo.SupplementaryDataFilenameSuffix], 'SupplementaryData');
    end

end
