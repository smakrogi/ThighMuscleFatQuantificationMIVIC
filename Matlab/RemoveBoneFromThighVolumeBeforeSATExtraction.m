function [NoBoneMask, BoneLabelVolume] = ...
    RemoveBoneFromThighVolumeBeforeSATExtraction...
    (OneLegVolume, strelSizeFactor, ...
    boneareaThreshold, ...
    standarddevCorticalBoneRemoval, ...
    rgIterationsCorticalBoneRemoval, ...
    seedingIterationsCorticalBoneRemoval)

infoString = sprintf(['Bone removal using bone marrow selection and ' ...
'cortical bone region growing.\n']);
fprintf( infoString );

midSlice = round(size(OneLegVolume, 3)/2);

% Gaussian smoothing.
OriginalOneLegVolume = OneLegVolume;
OneLegVolume = matitk('FGA',[0.75 9], single(OneLegVolume));

% Apply a small threshold to separate the low intensity voxels
% of air and cortical bone from the other tissues.
% NoCorticalBoneMask = OneLegVolume>threshold;

% % Use Otsu thresholding to separate muscle from fat.
% BoneandFatMask = matitk('SOT', [512] , single(OneLegVolume));
% BoneandFatMask = BoneandFatMask > 0;
% Multiple Otsu thresholding.
[CorticalBoneMask, MuscleMask] = matitk('FOMT', [2, 512] , single(OneLegVolume));
CorticalBoneMask = CorticalBoneMask > 0;
CorticalBoneMask= uint16(CorticalBoneMask);
MuscleMask = MuscleMask > 0;
MuscleMask= uint16(MuscleMask);
NoCorticalBoneMask = 1 - CorticalBoneMask;
BoneandFatMask = 1 - (CorticalBoneMask + MuscleMask);

% areaThreshold = 20;

% Apply connected component analysis.
% 1. find connected components
% Apply labeling in 3D.
LabelVolume = bwlabeln(BoneandFatMask, 6);

% Compute the region properties of the middle slice.
s = regionprops(LabelVolume(:, :, midSlice), 'Area', 'Centroid', 'Eccentricity',...
    'Perimeter', 'EulerNumber');
centroid =  cat(1, s.Centroid);
eccentricity = cat(1, s.Eccentricity);
area = cat(1, s.Area);
perimeter = cat(1, s.Perimeter);
eulerNumber = cat(1, s.EulerNumber);

% Estimate roundness and handle exceptions.
roundnessThreshold = 1.2;
roundness = 4*pi*area./perimeter.^2;
% Find objects with estimated perimeter = 0 that produce Infs and NaNs in roundness
% computation.
index = find( perimeter==0 );
% Change Infs and NaNs to 0.
roundness( index ) = 0;
roundnessMask = roundness <= roundnessThreshold;
roundness = roundness .* roundnessMask;

% Exclude subcutaneous region using Euler number.
eulerCandidates = eulerNumber == 1;

% Select regions with areas bigger than 1% of the total area.
areaThreshold = boneareaThreshold;  % 0.028;  % 0.02;
areaPercentage = area / sum(area);
areaCandidates = areaPercentage > areaThreshold;
candidateIndices = find( eulerCandidates.*areaCandidates );

% Use as criterion the roundness of remaining regions.
% boneLikelihood = roundness(candidateIndices);
boneLikelihood = roundness(candidateIndices);
% boneLikelihood = eulerNumber;
[temp, index] = sort(boneLikelihood, 1, 'descend');
clear temp;
boneLabel = candidateIndices(index(1));
BoneMarrowMask =  LabelVolume==boneLabel;
BoneMarrowMask = uint16(BoneMarrowMask);
fprintf('Bone marrow region features: Centroid: %d,%d \n', ...
    centroid(boneLabel, 1), centroid(boneLabel, 2));
fprintf('Roundness: %d, Eccentricity: %d, Area Percentage: %d\n', ...
    roundness(boneLabel), eccentricity(boneLabel), areaPercentage(boneLabel));

% Remove cortical bone using morphological operations and
% region growing.
[CorticalBoneMask, BoneMarrowMask, seeds] = GenerateBoneMasks(OneLegVolume, BoneMarrowMask, ...
    NoCorticalBoneMask, standarddevCorticalBoneRemoval, ...
    rgIterationsCorticalBoneRemoval, strelSizeFactor, seedingIterationsCorticalBoneRemoval);
% CorticalBoneMask = 1 - NoCorticalBoneMask;
CorticalBoneMask= uint16(CorticalBoneMask);

% Remove the identified bone.
NoBoneMask = uint16(OriginalOneLegVolume>0) - uint16(BoneMarrowMask) - uint16(CorticalBoneMask);
bonemarrowLabel = 3;
BoneLabelVolume = CorticalBoneMask + uint16(bonemarrowLabel * BoneMarrowMask);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

