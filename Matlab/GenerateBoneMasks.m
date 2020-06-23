function [CorticalBoneMask, CorrectedBoneMarrowMask, seeds] = ...
    GenerateBoneMasks(OneLegVolume, BoneMarrowMask,...
    NoCorticalBoneMask, std, rgIterations, ...
    strelSizeFactor, seedingIterations)
% Apply seeded region growing to produce the cortical bone region.
% The seeds come from basic thresholding.

strelSize = round(sqrt(numel(OneLegVolume(:,:,1))) / strelSizeFactor);
% DilatedBoneMarrowMask = BoneMarrowMask;
count = 1;
xIndex = []; yIndex = [];
seeds = [];
nSeeds = 0;
seedDifference = 0;
% seedingIterations = 2;
% se = strel('disk', strelSize);

infoString = sprintf('Segmenting cortical bone using: \n');
fprintf( infoString );
infoString = sprintf('Seeding iterations: %d, RG Iterations: %d, RG std: %d\n', ...
    seedingIterations, rgIterations, std);
fprintf( infoString );
infoString = sprintf('Structure element radius: %d\n', ...
    strelSize);
fprintf( infoString );

% while ((isempty(seeds) || seedDifference>0) || count <= seedingIterations )
while ( (isempty(seeds) || seedDifference>0) && count <= seedingIterations )
    % Mark a zone around the bone marrow.
    seeds = [];
    DilatedBoneMarrowMask = matitk('FBD', [count*strelSize, 1], single(BoneMarrowMask));
    DilatedBoneMarrowMask = uint16(DilatedBoneMarrowMask);
    TissueCandidates = DilatedBoneMarrowMask - BoneMarrowMask;
    % Multiply with the mask of cortical voxels.
    CorticalMaskCandidates = TissueCandidates .* uint16(1 - NoCorticalBoneMask );
%     MuscleMaskCandidates = TissueCandidates .* uint16( MuscleMask );
    I = find(CorticalMaskCandidates);
    [xIndex, yIndex, zIndex] = ind2sub(size(CorticalMaskCandidates), I);
	seedDifference = length(xIndex) - nSeeds;
	nSeeds = length(xIndex);
    for j=1:nSeeds
        if (zIndex(j) < size(BoneMarrowMask, 3)-1)
             seeds = [seeds, xIndex(j), yIndex(j), zIndex(j)];
        end
    end
	count = count + 1;
end

infoString = sprintf('Seeding iterations: %d\n', ...
    count-1);
fprintf( infoString );

% Apply region gowing to the cortical bone region.
% OneLegVolume = matitk('FCA',[5 0.0625 3], single(OneLegVolume));
% OneLegVolume = matitk('FGA',[1 9], single(OneLegVolume));
% OneLegVolume = matitk('FMEDIAN',[7 7 7], single(OneLegVolume));
% CorticalBoneMask = matitk('SCC', [std, rgIterations, 1], single(OneLegVolume), single([]), seeds);
% CorticalBoneMask = uint16(CorticalBoneMask);
% CorticalBoneMask = uint16((CorticalBoneMask + CorticalMaskCandidates) > 0);
% If this grows rtoo big, multiply with the cortical mask candidates.

CorticalBoneMask = uint16(CorticalMaskCandidates > 0);
% ResidualMuscleMask = uint16(MuscleMaskCandidates > 0);

% Apply closing to the bone mask.
seSize = 3.0;
BoneMask = CorticalBoneMask + BoneMarrowMask;
DilatedBoneMask = matitk('FBD', [seSize 1.0], uint8(BoneMask));
ClosedBoneMask = uint16(matitk('FBE', [seSize 1.0], uint8(DilatedBoneMask)));
CorticalBoneMask = uint16(ClosedBoneMask - BoneMarrowMask);

CorrectedBoneMarrowMask = uint16(BoneMarrowMask);
% CorrectedBoneMarrowMask = DilatedBoneMarrowMask - ...
%     CorticalBoneMask;
% CorrectedBoneMarrowMask = DilatedBoneMarrowMask - ...
%     CorticalBoneMask - ...
%     ResidualMuscleMask;

end