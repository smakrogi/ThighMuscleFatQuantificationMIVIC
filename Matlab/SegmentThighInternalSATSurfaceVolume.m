function [InternalSATSurfaceMask, OneLegVolume, BoundingBox, NoBoneMask] = ...
    SegmentThighInternalSATSurfaceVolume (subjectNumber, ExperimentInfo)
% SegmentInternalSubcutaneousSurfaceVolume
% 3T NIA/NIH, by S.K. Makrogiannis. 

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
% infoString = sprintf('First useable slice: %d, Last useable slice: %d\n', ExperimentInfo.FirstUseableSlice,...
%     ExperimentInfo.LastUseableSlice);
% fprintf( infoString );

infoString = sprintf('Reading input volume... \n');
fprintf( infoString );
% NonSuppressed = analyze75read([ExperimentInfo.dataPath ExperimentInfo.NonSuppressedFilename{subjectNumber}]);
InputVolume = load_untouch_nii([ExperimentInfo.dataPath ExperimentInfo.NonSuppressedFilename{subjectNumber}]);
% NonSuppressed = load_untouch_nii([ExperimentInfo.dataPath ExperimentInfo.NonSuppressedFilename{subjectNumber}], ...
%     [ExperimentInfo.FirstUseableSlice:ExperimentInfo.LastUseableSlice], 1);
InputVolume = InputVolume.img(:,:,ExperimentInfo.FirstUseableSlice:ExperimentInfo.LastUseableSlice);
% InputVolume = TransposeVolume(InputVolume);

% Gradient Anisotropic diffusion filtering.
%NonFatSuppressed = uint16( NonFatSuppressed(:,:,1:4:20) );
anisodiffusionParameters = [2 0.0375 3];
InputVolume = matitk('FGAD', anisodiffusionParameters, single(InputVolume));
InputVolume = int16( InputVolume );
fprintf('\n');

% Select one leg, corp axial slices.
OneLegVolume = zeros(size(InputVolume));
ObjectMap = zeros(size(InputVolume));
SnakeOptions = ...
    SetPrinceSnakeOptionsThigh(...
    ExperimentInfo, ...
    subjectNumber);
for i=1:size(InputVolume, 3)
    OneLegVolume(:,:,i) = SelectOneLeg(InputVolume(:,:,i), SnakeOptions);
end

% Compute bounding box and crop volume.
[BoundingBox] = FindSquareAxialBoundingBox(OneLegVolume);
OneLegVolume = OneLegVolume(BoundingBox.xmin:BoundingBox.xmax, ...
    BoundingBox.ymin:BoundingBox.ymax, :);
OneLegVolume = int16( OneLegVolume );


% Top hat transform to reduce inhomogeneity artifacts.
if(SnakeOptions.useTopHat)
    [OneLegVolume] = ...
        TopHatTransformMultiSlice(OneLegVolume, ...
        OneLegVolume, ...
        OneLegVolume, ...
        ExperimentInfo, ...
        SnakeOptions.strelSizeFactor, ...
        0);
end

% Segment out the subcutaneous fat using snakes in multiple scales.
switch(lower(ExperimentInfo.SATExtractionAlgorithm))
    case 'snake'
        strelSizeFactor = 128;
        strelSize = round(sqrt(prod(size(OneLegVolume(:,:,1)))) / strelSizeFactor);
        se = strel('disk', strelSize);
        InternalSATSurfaceMask = zeros(size(OneLegVolume));
        InternalSATSurfaceMask = uint8(InternalSATSurfaceMask);
        for i=1:size(OneLegVolume, 3)
            infoString = sprintf('\nSlice #%d \n', i);
            fprintf('%c',char(8*ones(1,length(infoString))));
            fprintf(infoString);
            InternalSATSurfaceMask(:,:,i) = ...
                SegmentInternalSubcutaneousSurface2D(single(OneLegVolume(:,:,i)), SnakeOptions);
            InternalSATSurfaceMask(:,:,i) = imerode(InternalSATSurfaceMask(:,:,i), se); %imopen
        end
    case 'gac'
        % Segment out the subcutaneous fat using level sets.
        InternalSATSurfaceMask = SegmentInternalSubcutaneousSurfaceLevelSets(OneLegVolume, ...
            ExperimentInfo);
end

% Bone removal algorithm.
[NoBoneMask, BoneLabelMask] = ...
    RemoveBoneFromThighVolumeBeforeSATExtraction...
    (OneLegVolume .* int16(InternalSATSurfaceMask), ...
    ExperimentInfo.strelSizeFactorBoneRemoval, ...
    ExperimentInfo.boneAreaThreshold, ...
    ExperimentInfo.standarddevCorticalBoneRemoval, ...
    ExperimentInfo.rgIterationsCorticalBoneRemoval, ...
    ExperimentInfo.seedingIterationsCorticalBoneRemoval);

if ( ExperimentInfo.commandlineMode )
    AxialSliceMontage(BoneLabelMask, ExperimentInfo.sliceStep);
    saveas(gcf, [ExperimentInfo.SubjectID{subjectNumber} '_bone_label_map'], 'png')
    
    %InternalSATSurfaceMask = InternalSATSurfaceMask .* uint8(NoBoneMask);
    
    % Display orthoviews of input and output.
    % DisplayOrthoSlices(OneLegVolume, OneLegVolume.*int16(InternalSATSurfaceMask));
    % saveas(gcf, [ExperimentInfo.SubjectID{subjectNumber} '_sat_segmentation_ortho'], 'png')
    % Create slice montage.
    AxialSliceMontage(OneLegVolume, OneLegVolume.*int16(InternalSATSurfaceMask), ExperimentInfo.sliceStep);
    saveas(gcf, [ExperimentInfo.SubjectID{subjectNumber} '_sat_segmentation_montage'], 'png')
    
    % % Display orthoviews of bone segmentation.
    % DisplayOrthoSlices(OneLegVolume, OneLegVolume.*int16(BoneLabelMask));
    % saveas(gcf, [ExperimentInfo.SubjectID{subjectNumber} '_bone_segmentation_ortho'], 'png')
    % % Create slice montage.
    % AxialSliceMontage(OneLegVolume, OneLegVolume.*int16(BoneLabelMask), ExperimentInfo.sliceStep);
    % saveas(gcf, [ExperimentInfo.SubjectID{subjectNumber} '_bone_segmentation_montage'], 'png')
end

% Save ouput in Nifti format.

% % Save output in Analyze 7.5 format.
% status = MICE_SaveDICOMvolume(InternalSATSurfaceMask, 'anlz', DICOMinfoStruct, 'InternalSATSurfaceMask');
% status = 3TNIA_SaveAnalyzevolume(InternalSATSurfaceMask, format, AnalyzeinfoStruct, filename);

% InternalSATSurfaceMask = TransposeVolume(InternalSATSurfaceMask);
InternalSATSurfaceMaskNii = make_nii(InternalSATSurfaceMask, voxelSize, origin, 4, ...
    description);
save_nii(InternalSATSurfaceMaskNii, [ExperimentInfo.SubjectID{subjectNumber} '_InternalSATSurfaceMask.hdr']);

% BoneLabelMask = TransposeVolume(BoneLabelMask);
BoneLabelMaskNii = make_nii(BoneLabelMask, voxelSize, origin, 4, ...
    description);
save_nii(BoneLabelMaskNii, [ExperimentInfo.SubjectID{subjectNumber} '_BoneLabelVolume.hdr']);

% OneLegVolume = TransposeVolume(OneLegVolume);
OneLegVolumeNii = make_nii(OneLegVolume, voxelSize, origin, 4, ...
    description);
save_nii(OneLegVolumeNii, [ExperimentInfo.SubjectID{subjectNumber} '_OneLegVolume.hdr']);

% Save some more data in .mat file (may be used to save free parameters in the
% future.
SaveSupplementaryData(subjectNumber, ExperimentInfo, BoundingBox);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveSupplementaryData(subjectNumber, ExperimentInfo, BoundingBox)

imagefileprefix = [ExperimentInfo.SubjectID{subjectNumber}, '_', ...
    ExperimentInfo.DistanceLearningMethod, '_', ...
    ExperimentInfo.ClusteringMethod];

% Check if a supplementary data file already exists.
fid = fopen([imagefileprefix,  ExperimentInfo.SupplementaryDataFilenameSuffix]);
if fid~=-1
    fclose( fid );
    load([imagefileprefix,  ExperimentInfo.SupplementaryDataFilenameSuffix], 'SupplementaryData');
end

SupplementaryData.BoundingBox = BoundingBox;

save([imagefileprefix,  ExperimentInfo.SupplementaryDataFilenameSuffix], ...
    'SupplementaryData');

end
