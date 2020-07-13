function InternalSATSurfaceMask = SegmentInternalSubcutaneousSurfaceLevelSets(OneLegVolume, ExperimentInfo)
% Use level sets for segmentation of subcutaneous area.

debug = 0;

% Apply padding.
originalSize = size(OneLegVolume);
PaddedOneLegVolume = zeros ([originalSize(1) originalSize(2) originalSize(3)+2 ]);
PaddedOneLegVolume(:,:,1) = OneLegVolume(:,:,1);
PaddedOneLegVolume(:,:,2:originalSize(3)+1)  = OneLegVolume;
PaddedOneLegVolume(:,:,originalSize(3)+2) = OneLegVolume(:,:,originalSize(3));

% Gradient magnitude.
% GM = matitk('FGM', [], single(OneLegVolume));
GMRG = matitk('FGMRG', [ExperimentInfo.snakeGaussianSigma], single(PaddedOneLegVolume));

% Edge potential.
% PI = matitk('FSN', [0 1 -75 180], single(GM));
%PI = matitk('FSN', [0 1 -35 80], single(GMRG));
PI = matitk('FSN', [0 1 -3 10], single(GMRG));

% Fast marching to generate initial thigh mask.
% MS = matitk('SFM', [200], single(OneLegVolume));
MS = matitk('SFM', [200], single(PaddedOneLegVolume), single([]), round(size(PaddedOneLegVolume)/2));
MS = MS <= 5;

% Generate signed distance map.
DM = matitk('FDM', [], single(MS));
DM2 = matitk('FDM', [], single(1-MS));
DMF = DM.^0.5 - DM2.^0.5;

% Apply geodesic active contours and shape detection level sets.
initialDistance = 60; % 60;
% GAC
% [propagationScaling,CurvatureScaling, AdvectionScaling, MaximumRMSError, MaxIteration]
% InternalSATSurfaceMask = matitk('SGAC', [-1.75 0.75 2.5 0.02 550], single(PI), single(DMF+initialDistance));
% InternalSATSurfaceMask = matitk('SGAC', [-1.0 0.1 1.0 0.02 750], single(PI), single(DMF+initialDistance));
InternalSATSurfaceMask = matitk('SGAC', [-1.0 0.1 1.0 0.02 350], single(PI), single(DMF+initialDistance));
InternalSATSurfaceMask = InternalSATSurfaceMask > 0;
% SDLS
% InternalSATSurfaceMask = matitk('SSDLS', [-2.75 1.25 0.02 750], single(PI), single(DMF+initialDistance));
if debug
    i = round(originalSize/2)+1;
    figure, imagesc(single(InternalSATSurfaceMask(:,:,i)).*single(PaddedOneLegVolume(:,:,i))), axis image, colorbar
end

% Remove padding from internal surface mask.
InternalSATSurfaceMask = InternalSATSurfaceMask(:,:,2:originalSize(3)+1);
InternalSATSurfaceMask = uint8(InternalSATSurfaceMask);