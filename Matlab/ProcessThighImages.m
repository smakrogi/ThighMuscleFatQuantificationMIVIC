function [TissueLabelVolumeOutput, tissueLabels, DistributionModels] = ...
    ProcessThighImages(subjectNumber, ExperimentInfo, varargin)
% syntax: [TissueLabelVolumeOutput, fatLabel] = ProcessThighImages...
% (subjectNumber, ExperimentInfo, varargin);
% 3T NIA/NIH, by S.K. Makrogiannis. 

% Read-in the volumes.
% info = analyze75info([ExperimentInfo.dataPath ExperimentInfo.filenames{1}]);
info = load_nii_hdr([ExperimentInfo.dataPath ExperimentInfo.NonSuppressedFilename{subjectNumber}]);
voxelSize = info.dime.pixdim(2:4);
origin = [info.hist.qoffset_x info.hist.qoffset_y info.hist.qoffset_z];
description = info.hist.descrip;
% infoString = sprintf('Volume descriptor: %s\n', info.Filename);
% fprintf( infoString );
% infoString = sprintf('Subject ID: %s\n', info.PatientID);
% fprintf( infoString );
% infoString = sprintf('Volume descriptor: %s\n', info.Descriptor);
% fprintf( infoString );
% infoString = sprintf('Exposure date, time: %s, %s\n', info.ExposureDate, info.ExposureTime);
% fprintf( infoString );
% infoString = sprintf('First useable slice: %d, ...
% Last useable slice: %d\n', ExperimentInfo.FirstUseableSlice,...
%     ExperimentInfo.LastUseableSlice);
% fprintf( infoString );

infoString = sprintf('Reading Non suppressed volume... \n');
fprintf( infoString );
% NonSuppressed = analyze75read([ExperimentInfo.dataPath ExperimentInfo.NonSuppressedFilename{subjectNumber}]);
NonSuppressed = load_untouch_nii([ExperimentInfo.dataPath ExperimentInfo.NonSuppressedFilename{subjectNumber}]);
% NonSuppressed = load_untouch_nii([ExperimentInfo.dataPath ExperimentInfo.NonSuppressedFilename{subjectNumber}], ...
% [ExperimentInfo.FirstUseableSlice:ExperimentInfo.LastUseableSlice], 1);
NonSuppressed = NonSuppressed.img(:,:,...
    ExperimentInfo.FirstUseableSlice:ExperimentInfo.LastUseableSlice);
% NonSuppressed = TransposeVolume(NonSuppressed);

infoString = sprintf('Reading Fat suppressed volume... \n');
fprintf( infoString );
% FatSuppressed = analyze75read([ExperimentInfo.dataPath ExperimentInfo.FatSuppressedFilename{subjectNumber}]);
FatSuppressed = load_untouch_nii([ExperimentInfo.dataPath ExperimentInfo.FatSuppressedFilename{subjectNumber}]);
FatSuppressed = FatSuppressed.img(:,:,...
    ExperimentInfo.FirstUseableSlice:ExperimentInfo.LastUseableSlice);
% FatSuppressed = TransposeVolume(FatSuppressed);
resizeFactor = size(NonSuppressed,1) / size(FatSuppressed,1);
if( numel(FatSuppressed) ~= numel(NonSuppressed) )
    fprintf('Need to resize fat volume to match non-suppressed.\n');
    FatSuppressed = ResizeVolumeMatrixSize(FatSuppressed, NonSuppressed);
end
FatSuppressed = ShiftVolumeLR(FatSuppressed, ExperimentInfo.LRShiftFS * resizeFactor);

infoString = sprintf('Reading Water suppressed volume... \n');
fprintf( infoString );
%WaterSuppressed = analyze75read([ExperimentInfo.dataPath ExperimentInfo.WaterSuppressedFilename{subjectNumber}]);
WaterSuppressed = load_untouch_nii([ExperimentInfo.dataPath ExperimentInfo.WaterSuppressedFilename{subjectNumber}]);
WaterSuppressed = WaterSuppressed.img(:,:,...
    ExperimentInfo.FirstUseableSlice:ExperimentInfo.LastUseableSlice);
% WaterSuppressed = TransposeVolume(WaterSuppressed);
resizeFactor = size(NonSuppressed,1) / size(WaterSuppressed,1);
if( numel(WaterSuppressed) ~= numel(NonSuppressed) )
    fprintf('Need to resize water volume to match non-suppressed.\n');
    WaterSuppressed = ResizeVolumeMatrixSize(WaterSuppressed, NonSuppressed);
end
WaterSuppressed = ShiftVolumeLR(WaterSuppressed, ExperimentInfo.LRShiftWS * resizeFactor);

% Apply non-linear diffusion filtering.
anisodiffusionParameters = [2 0.0375 3];  % [5 0.0375 3];
NonSuppressed = matitk('FGAD', anisodiffusionParameters, single(NonSuppressed));
NonSuppressed = uint16( NonSuppressed );
FatSuppressed = matitk('FGAD', anisodiffusionParameters, single(FatSuppressed));
FatSuppressed = uint16( FatSuppressed );
WaterSuppressed = matitk('FGAD', anisodiffusionParameters, single(WaterSuppressed));
WaterSuppressed = uint16( WaterSuppressed );

% Top hat transform to alleviate inhomogeneity artifacts.
if ExperimentInfo.tophatTransform
    [NonSuppressed, ...
        FatSuppressed, ...
        WaterSuppressed] = ...
        TopHatTransformMultiSlice(NonSuppressed, ...
        FatSuppressed, ...
        WaterSuppressed, ...
        ExperimentInfo, ...
        ExperimentInfo.strelSizeFactorTopHat, ...
        1);
end

% If SAT has been cut off then use the extra info.
if(~isempty(varargin))
    % Read-in InternalSATSurfaceMask and pass bounding box info.
    InternalSATSurfaceVolume = load_untouch_nii([varargin{1}]);
    InternalSATSurfaceVolume = uint16(InternalSATSurfaceVolume.img);
    OneLegVolume = load_untouch_nii([varargin{2}]);
    OneLegVolume = OneLegVolume.img;
    % Read supplementary data file.
    SupplementaryData = LoadSupplementaryData(subjectNumber, ExperimentInfo);
    BoundingBox = SupplementaryData.BoundingBox;    
    
    % Crop volumes to include one leg only.
    NonSuppressed = NonSuppressed(BoundingBox.xmin:BoundingBox.xmax, ...
        BoundingBox.ymin:BoundingBox.ymax, :);
    FatSuppressed = FatSuppressed(BoundingBox.xmin:BoundingBox.xmax, ...
        BoundingBox.ymin:BoundingBox.ymax, :);
    WaterSuppressed = WaterSuppressed(BoundingBox.xmin:BoundingBox.xmax, ...
        BoundingBox.ymin:BoundingBox.ymax, :);

    % Generate SAT mask and compute WS and FS modes/averages for scaling.
    stats = ComputeAverageSATIntensities(OneLegVolume, ...
        InternalSATSurfaceVolume, ...
        NonSuppressed, ...
        FatSuppressed, ...
        WaterSuppressed);
%     stats = ComputeSATModes(OneLegVolume, ...
%         InternalSATSurfaceVolume, ...
%         NonSuppressed, ...
%         FatSuppressed, ...
%         WaterSuppressed);
    
    % Scale intensities w.r.t. the SAT intensity mean.
    NonSuppressed = single(NonSuppressed) / stats.nonsuppressed_sat_mean;    
    FatSuppressed = single(FatSuppressed) / stats.fatsuppressed_sat_mean;
    WaterSuppressed = single(WaterSuppressed) / stats.watersuppressed_sat_mean;
    
    % Mask the input volumes with the no SAT mask.
    NonSuppressed = (NonSuppressed) .* logical(InternalSATSurfaceVolume);
    FatSuppressed = (FatSuppressed) .* logical(InternalSATSurfaceVolume);
    WaterSuppressed = (WaterSuppressed) .* logical(InternalSATSurfaceVolume);

    % Remove cortical bone.
    if(ExperimentInfo.removeBone)
        % Read bone label volume.
        BoneLabelVolume = load_untouch_nii([varargin{3}]);
        BoneLabelVolume = BoneLabelVolume.img;
        % Create NoBoneMask.
        NoBoneMask = uint16(InternalSATSurfaceVolume>0) - uint16(BoneLabelVolume > 0);
        % Apply masking to MRI volumes.
        InternalSATSurfaceVolume = InternalSATSurfaceVolume .* NoBoneMask;
        NonSuppressed = NonSuppressed .* logical(NoBoneMask);
        FatSuppressed = FatSuppressed .* logical(NoBoneMask);
        WaterSuppressed = WaterSuppressed .* logical(NoBoneMask);
        clear NoBoneMask;
    end
    
    clear OneLegVolume;
end


% Write fat- and water-suppressed nifti volumes without SAT and bone.
if(ExperimentInfo.warpWSandFStoNS || ExperimentInfo.commandlineMode)
    % Write volumes to hard drive in nifti format.
    % NonSuppressed = TransposeVolume(NonSuppressed);
    NonSuppressedNii = make_nii(NonSuppressed, voxelSize, origin, 16, ...
        description);
    nonsuppressednosatFilename = [ExperimentInfo.SubjectID{subjectNumber} '_NonSuppressedNoSAT.hdr'];
    save_nii(NonSuppressedNii, nonsuppressednosatFilename);

    % FatSuppressed = TransposeVolume(FatSuppressed);
    FatSuppressedNii = make_nii(FatSuppressed, voxelSize, origin, 16, ...
        description);
    fatsuppressednosatFilename = [ExperimentInfo.SubjectID{subjectNumber} '_FatSuppressedNoSAT.hdr'];
    save_nii(FatSuppressedNii, fatsuppressednosatFilename);

    % WaterSuppressed = TransposeVolume(WaterSuppressed);
    WaterSuppressedNii = make_nii(WaterSuppressed, voxelSize, origin, 16, ...
        description);
    watersuppressednosatFilename = [ExperimentInfo.SubjectID{subjectNumber} '_WaterSuppressedNoSAT.hdr'];
    save_nii(WaterSuppressedNii, watersuppressednosatFilename);

    clear NonSuppressedNii FatSuppressedNii WaterSuppressedNii;
end

% Warp fat- and water-suppressed volumes to non-suppressed volume.
if(ExperimentInfo.warpWSandFStoNS)
    FatSuppressed = WarpToNonSuppressedVolume(fatsuppressednosatFilename, nonsuppressednosatFilename);
    WaterSuppressed = WarpToNonSuppressedVolume(watersuppressednosatFilename, nonsuppressednosatFilename);
end


% Downsample.
%rescaleFactor = 1/8;
if( ExperimentInfo.rescaleFactor ~= 1)
    for i=1:size(NonSuppressed, 3)
        SubNonSuppressed(:,:,i) = ...
            imresize(NonSuppressed(:,:,i), ...
            ExperimentInfo.rescaleFactor, 'Method', 'nearest');
        SubFatSuppressed(:,:,i) = ...
            imresize(FatSuppressed(:,:,i), ...
            ExperimentInfo.rescaleFactor, 'Method', 'nearest');
        SubWaterSuppressed(:,:,i) = ...
            imresize(WaterSuppressed(:,:,i), ...
            ExperimentInfo.rescaleFactor, 'Method', 'nearest');
        SubInternalSATSurfaceVolume(:,:,i) = ...
            imresize(InternalSATSurfaceVolume(:,:,i), ...
            ExperimentInfo.rescaleFactor, 'Method', 'nearest');
    end
else
    SubNonSuppressed = NonSuppressed;
    SubFatSuppressed = FatSuppressed;
    SubWaterSuppressed = WaterSuppressed;
    SubInternalSATSurfaceVolume = InternalSATSurfaceVolume;
end

% Keep original size and clear some more memory.
originalVolumeSize = size(NonSuppressed);
clear NonSuppressed FatSuppressed WaterSuppressed;



% Form feature vectors.
NonSuppressedVector = SubNonSuppressed(:);
FatSuppressedVector = SubFatSuppressed(:);
WaterSuppressedVector = SubWaterSuppressed(:);
InternalSATSurfaceVector = SubInternalSATSurfaceVolume(:);


if( ExperimentInfo.clusterNonZeroValuesOnly )
    [I,J] = find(InternalSATSurfaceVector);
    NoZerosNonSuppressedVector = NonSuppressedVector(I);
    NoZerosFatSuppressedVector = FatSuppressedVector(I);
    NoZerosWaterSuppressedVector = WaterSuppressedVector(I);
end

if ( ExperimentInfo.commandlineMode )
    figure, subplot(131), hist(single(NoZerosNonSuppressedVector),100), ...
        title('Non suppressed intensities'), grid on, axis square;
    subplot(132), hist(single(NoZerosFatSuppressedVector),100), ...
        title('Fat suppressed intensities'), grid on, axis square;
    subplot(133), hist(single(NoZerosWaterSuppressedVector),100), ...
        title('Water suppressed intensities'), grid on, axis square;
    saveas(gcf, [ExperimentInfo.SubjectID{subjectNumber} '_t1w_fat_water_histograms.png']);
end

% Form Data Matrix.
if ( ExperimentInfo.clusterNonZeroValuesOnly )
    % dataMatrix = [NoZerosNonSuppressedVector, NoZerosFatSuppressedVector, NoZerosWaterSuppressedVector];
     dataMatrix = [NoZerosFatSuppressedVector, NoZerosWaterSuppressedVector];
%    dataMatrix = [NoZerosNonSuppressedVector];
else
    % dataMatrix = [NonSuppressedVector, FatSuppressedVector, WaterSuppressedVector];
    dataMatrix = [FatSuppressedVector, WaterSuppressedVector];
end
%dataMatrix = single( dataMatrix );
dataMatrix = double( dataMatrix );

% Analyze data.
if ( ExperimentInfo.clusterNonZeroValuesOnly )
    [NoZerosTissueIdx, NoZerosOutputDataMatrix, DistributionModels] = ...
        AnalyzeThighDataInFeatureSpace( subjectNumber, dataMatrix, ExperimentInfo );
	TissueIdx = zeros(size(NonSuppressedVector));
    OutputDataMatrix = zeros(size(NonSuppressedVector));
    
    for i=1:length(NoZerosTissueIdx)
        TissueIdx(I(i),J(i)) = NoZerosTissueIdx(i);
        OutputDataMatrix(I(i),J(i)) = NoZerosOutputDataMatrix(i, 1);
    end

    clear NoZeros*;
else
    TissueIdx = AnalyzeThighDataInFeatureSpace( subjectNumber, dataMatrix, ExperimentInfo );
end

% Form segmentation mask.
TissueLabelVolume = reshape( TissueIdx, size(SubNonSuppressed) );
TissueLabelVolume = uint16(TissueLabelVolume);


% Over-sampling.
if( ExperimentInfo.rescaleFactor ~= 1)
    TissueLabelVolumeOutput = zeros(originalVolumeSize);
    for i=1:size(TissueLabelVolume, 3)
        TissueLabelVolumeOutput(:,:,i) = imresize(TissueLabelVolume(:,:,i), ...
            [originalVolumeSize(1) originalVolumeSize(2)], 'Method', 'nearest');
    end
    TissueLabelVolumeOutput = round(TissueLabelVolumeOutput);
else 
    TissueLabelVolumeOutput = TissueLabelVolume;
end
TissueLabelVolumeOutput = uint16(TissueLabelVolumeOutput);

% clear  TissueLabelVolume;


% Try to remove bone marrow after tissue segmentation.
tissueLabels = DetermineFatAndMuscleLabels(TissueLabelVolume, ...
    SubWaterSuppressed, ExperimentInfo.nClusters);


% Write volumes to hard drive in nifti format.
imagefileprefix = [ExperimentInfo.SubjectID{subjectNumber}, '_', ...
    ExperimentInfo.DistanceLearningMethod, '_', ...
    ExperimentInfo.ClusteringMethod];
% TissueLabelVolumeOutput = TransposeVolume(TissueLabelVolumeOutput);
TissueLabelVolumeNii = make_nii(TissueLabelVolumeOutput, voxelSize, origin, 4, ...
    description);
save_nii(TissueLabelVolumeNii, [imagefileprefix, '_TissueLabelVolume.hdr']);


% Also write data samples and tissue labels to drive.
NewSupplementaryData.DataMatrix = DistributionModels.Samples;
NewSupplementaryData.tissueLabels = tissueLabels;
SaveSupplementaryData(subjectNumber, ExperimentInfo, NewSupplementaryData)

% Display input/output images.
if ( ExperimentInfo.commandlineMode )
    % Generate volume masks.
    % sliceStep = 4;
    SegmentationMask = zeros([size(SubNonSuppressed), ExperimentInfo.nClusters]);
    SegmentationMask = single(SegmentationMask);
    
    for count=1:ExperimentInfo.nClusters
        SegmentationMask(:,:,:,count) = SubNonSuppressed .* ...
            logical(TissueLabelVolume==count);
    end
    
    % Display orthoviews of input and output.
    % DisplayOrthoSlices(SubNonSuppressed, uint16(SegmentationMask));
    % saveas(gcf, [ imagefileprefix, '_segmentation_orthoviews.png'], 'png')
    % Create slice montage.
    AxialSliceMontage(SubNonSuppressed, uint16(SegmentationMask), ExperimentInfo.sliceStep);
    saveas(gcf, [imagefileprefix, '_segmentation_montage'], 'png')
    
    for count=1:ExperimentInfo.nClusters
        SegmentationMask(:,:,:,count) = SubWaterSuppressed .* ...
            logical(TissueLabelVolume==count);
    end
    % DisplayOrthoSlices(SubWaterSuppressed, uint16(SegmentationMask));
    % saveas(gcf, [ imagefileprefix, '_ws_segmentation_orthoviews.png'], 'png')
    AxialSliceMontage(SubWaterSuppressed, uint16(SegmentationMask), ExperimentInfo.sliceStep);
    saveas(gcf, [imagefileprefix, '_ws_segmentation_montage'], 'png')
    
    for count=1:ExperimentInfo.nClusters
        SegmentationMask(:,:,:,count) = SubFatSuppressed .* ...
            logical(TissueLabelVolume==count);
    end
    % DisplayOrthoSlices(SubFatSuppressed, uint16(SegmentationMask));
    % saveas(gcf, [ imagefileprefix, '_fs_segmentation_orthoviews.png'], 'png')
    AxialSliceMontage(SubFatSuppressed, uint16(SegmentationMask), ExperimentInfo.sliceStep);
    saveas(gcf, [imagefileprefix, '_fs_segmentation_montage'], 'png')
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutputVolume = ResizeVolumeMatrixSize(InputVolume, ReferenceVolume)
% Resize matrix size before creating the feature vectors.
matrixSize = [size(ReferenceVolume,1) size(ReferenceVolume,2)];
OutputVolume = ReferenceVolume;

for i=1:size(ReferenceVolume,3)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tissueLabels = DetermineFatAndMuscleLabels(TissueLabelVolume, WaterSuppressed, nClusters)

% Allocate mean vector.
meanIntensity = zeros(length(nClusters));

% Find which label corresponds to fat+bone by regionprops.
for i=1:nClusters
    TissueMask = uint16(TissueLabelVolume==i);
    TissueMap =  WaterSuppressed .* logical(TissueMask);
    [temp1,temp2,V] = find(TissueMap(:));
    clear temp1; clear temp2;
    meanIntensity(i) = mean(V);
end
[temp,tissueLabels.fat] = max(meanIntensity);
[temp,tissueLabels.muscle] = min(meanIntensity);
clear temp;
clear TissueMap;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stats = ComputeSATModes(OneLegVolume, ...
    InternalSATSurfaceVolume, ...
    NonSuppressed, ...
    FatSuppressed, ...
    WaterSuppressed)

    % Use OneLegVolume and InternalSATSurfaceVolume to create a SAT mask.
    InternalSATSurfaceVolume = InternalSATSurfaceVolume == 1;
    SATVolume = uint16(OneLegVolume) .*uint16(1-InternalSATSurfaceVolume);
    
    % Compute averages over SAT region.
    I_SAT = find(SATVolume(:)>0);
    SAT_NonSuppressed = single(NonSuppressed(I_SAT));
    stats.nonsuppressed_sat_mode_value = ComputeMode(SAT_NonSuppressed);
    SAT_FatSuppressed = single(FatSuppressed(I_SAT));
    stats.fatsuppressed_sat_mode_value = ComputeMode(SAT_FatSuppressed);
    SAT_WaterSuppressed = single(WaterSuppressed(I_SAT));
    stats.watersuppressed_sat_mode_value = ComputeMode(SAT_WaterSuppressed);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mode_value = ComputeMode(Variable)
    % Computes the value of variable mode.
    nbins = 128;
    step = ( max(Variable) - min(Variable) ) / nbins;
    edges = min(Variable):step:max(Variable);
    [count, bin] = histc(Variable, edges);
    first_mode = mode(bin);
    mode_value = ( edges(first_mode)+edges(first_mode+1) ) / 2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stats = ComputeAverageSATIntensities(OneLegVolume, ...
    InternalSATSurfaceVolume, ...
    NonSuppressed, ...
    FatSuppressed, ...
    WaterSuppressed)

    % Use OneLegVolume and InternalSATSurfaceVolume to create a SAT mask.
    InternalSATSurfaceVolume = InternalSATSurfaceVolume == 1;
    SATVolume = uint16(OneLegVolume) .*uint16(1-InternalSATSurfaceVolume);
    
    % Compute averages over SAT region.
    I_SAT = find(SATVolume(:)>0);
    SAT_NonSuppressed = NonSuppressed(I_SAT);
    stats.nonsuppressed_sat_mean = mean(SAT_NonSuppressed);
    SAT_FatSuppressed = FatSuppressed(I_SAT);
    stats.fatsuppressed_sat_mean = mean(SAT_FatSuppressed);
    SAT_WaterSuppressed = WaterSuppressed(I_SAT);
    stats.watersuppressed_sat_mean = mean(SAT_WaterSuppressed);
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
        load([imagefileprefix,  ExperimentInfo.SupplementaryDataFilenameSuffix], ...
            'SupplementaryData');
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveSupplementaryData(subjectNumber, ExperimentInfo, NewSupplementaryData)

imagefileprefix = [ExperimentInfo.SubjectID{subjectNumber}, '_', ...
    ExperimentInfo.DistanceLearningMethod, '_', ...
    ExperimentInfo.ClusteringMethod];

% Check if a supplementary data file already exists.
fid = fopen([imagefileprefix,  ExperimentInfo.SupplementaryDataFilenameSuffix]);
if fid~=-1
    fclose( fid );
    load([imagefileprefix,  ExperimentInfo.SupplementaryDataFilenameSuffix], ...
        'SupplementaryData');
end

SupplementaryData.DataMatrix = NewSupplementaryData.DataMatrix;
SupplementaryData.tissueLabels = NewSupplementaryData.tissueLabels;

save([imagefileprefix,  ExperimentInfo.SupplementaryDataFilenameSuffix], ...
    'SupplementaryData');

end
