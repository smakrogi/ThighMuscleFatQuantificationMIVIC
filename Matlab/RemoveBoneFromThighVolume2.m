function [NoBoneMask, BoneMask] = ...
    RemoveBoneFromThighVolume2...
    (TissueLabelVolume, fatLabel, areaThreshold)
% Remove bone marrow.

% Algorithm parameters.
% areaThreshold = 20;

% Use the tissue labels to extract fat+bone (label==2).
FatMask = int16(TissueLabelVolume == fatLabel);
BoneMask = FatMask;

% Apply connected component analysis.
% 1. find connected components
midSlice = round(size(TissueLabelVolume, 3)/2);

% Apply labeling in 3D.
ForegroundMask = TissueLabelVolume > 0;
LabelVolume = bwlabeln(ForegroundMask);

% Compute the region properties of the middle slice.
s = regionprops(LabelVolume(:, :, midSlice), 'Area', 'Centroid', 'Eccentricity',...
    'Perimeter', 'EulerNumber');
eccentricity = cat(1, s.Eccentricity);
area = cat(1, s.Area);
perimeter = cat(1, s.Perimeter);
eulerNumber = cat(1, s.EulerNumber);
roundness = 4*pi*area./perimeter.^2;
boneLikelihood = roundness;
boneLikelihood = eulerNumber;
[~,index] = sort(boneLikelihood, 1, 'descend');
boneLabel = index(1);
BoneMask =  LabelVolume==boneLabel;

% for i=1:size(FatMask, 3)
%     
%     %     if i==midSlice
%     %         figure, subplot(131), imagesc(OneLegVolume(:,:,i)),
%     %         colormap(gray), title('Original'), colorbar, axis image
%     %     end
%     
%     %
%     %FatMask(:,:,i) = imfill(FatMask(:,:,i),'holes');
%     FatSliceMask  = bwlabel(FatMask(:,:,i), 4);
%     
%     % Region properties: area, centroid, perimeter, eccentricity.
%     s = regionprops(FatSliceMask, 'Area', 'Centroid', 'Eccentricity',...
%         'Perimeter');
%     eccentricity = cat(1, s.Eccentricity);
%     area = cat(1, s.Area);
%     perimeter = cat(1, s.Perimeter);
%     roundness = 4*pi*area./perimeter.^2;
%     % Find regions bigger than one voxel and compute bone likelihood
%     % by multiplying roundness and area.
%     I = find(area>areaThreshold);
%     boneLikelihood = roundness(I).*area(I);
%     boneLikelihood = roundness(I);
%     [~,index] = sort(boneLikelihood, 1, 'descend');
%     
%     % Select top two.
%     % Select the left one from the top two.
%     % determine as bone the one that approahces a circular shape.
%     boneLabel = I(index(1));
%     BoneMask(:,:,i) = FatSliceMask==boneLabel;
% end

    % Remove the identified bone.
    NoBoneMask = TissueLabelVolume - uint16(BoneMask * fatLabel);
    
end

