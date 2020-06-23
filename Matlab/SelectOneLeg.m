function [OutputImage, objectMap] = SelectOneLeg(InputImage, options)

% Apply threshold to binarize.
% threshold = OtsuThreshold(InputImage);
threshold = options.fbthreshold * max(InputImage(:));
BinaryImage = InputImage > threshold;

% % Apply morphological erosion to separate possibly contiguous thighs.
% strelSizeFactor = 64;
% strelSize = round(sqrt(prod(size(InputImage))) / strelSizeFactor);
% se = strel('disk', strelSize);
% BinaryImage = imerode(BinaryImage, se);

% Connected components.
LabelImage  = bwlabel(BinaryImage, 4);

% Region properties: size, centroid.
s = regionprops( LabelImage, 'Area', 'Centroid' );
regionSize = cat(1, s.Area);
[temp,index] = sort(regionSize, 1, 'descend');
clear temp;

% Select top two.
% Select the left or right one from the top two.
planeCoords = cat(1, s.Centroid);
planeCoords = planeCoords(index(1:2),:);

if( strcmp(options.leg, 'left') )
    [temp, oneLeg] = min( planeCoords(:,2) );   % before transpose operations: 1
    clear temp;
else
    [temp, oneLeg] = max( planeCoords(:,2) );  % before transpose operations: 1
    clear temp;
end

objectMap = LabelImage==index(oneLeg);

% Get the boundary coordinates and generate a mask.
% Find contour coordinates of all objects (the main one and holes) and ...
% pick the first one that corresponds to the exterior boundary.
[temp, objectMap, nObjects] = bwboundaries(objectMap);
clear temp;
% Thresholding to merge main object with holes.
objectMask = objectMap > 0;

% % Morphological closing to remove small holes.
% objectMask = imclose(objectMask, se);

OutputImage = InputImage .* int16(objectMask);

end