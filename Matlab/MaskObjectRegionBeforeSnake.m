function [ObjectMask] = ...
    MaskObjectRegionBeforeSnake(inputImage, ...
    intensityRange, ...
    seSize, ...
    debug, ...
    varargin)
% Foreground/background segmentation.

% Define subcutaneous adipose initial contour
% through morphological operations.
% strelSizeFactor = 512;
% strelSize = round(sqrt(numel(inputImage(:,:,1))) / strelSizeFactor);
se = strel('disk', seSize);

% Thresholding.
if (nargin>4)
    % Multiply auto threshold by factor.
    autoThreshold = graythresh(inputImage);
    thresholdFactor = varargin{1};
    threshold = autoThreshold * thresholdFactor;
else
    threshold = 0.2; 
end
thresholdedImage = inputImage>(threshold);

% Display separated trunk.
if( debug )
    figure(2),imagesc(thresholdedImage),axis image, colormap gray,colorbar, ...
        title('Thresholded')
    % saveas(gcf, 'thresholded', 'png')
end
fprintf('Foreground/background threshold: % 2.3f\n', threshold);


% Morphological operations--opening.
openedImage = imopen(thresholdedImage, se);

% Find connected components.
labelImage = bwlabel( openedImage, 4 );

if( debug )
    figure(3),imagesc(labelImage),axis image, colormap gray,colorbar, ...
        title('Connected components')
end

a = regionprops(labelImage, 'area');
areas = cat(1, a.Area);

% Select largest region.
[temp, maxLabel] = max( areas );
clear temp;
largestRegionMask = labelImage==maxLabel;

if( debug )
    figure(4),imagesc(largestRegionMask),axis image, colormap gray,colorbar, ...
        title('Largest Region')
end

OpenedLargestRegionMask = largestRegionMask;
if( debug )
    [temp,temp2,V] = find(inputImage(:).*OpenedLargestRegionMask(:));
    clear temp, clear temp2;
    figure(2), subplot(121),imagesc(intensityRange*inputImage.*OpenedLargestRegionMask),...
        axis image, colormap gray,colorbar, title('Opened Largest Region')
    subplot(122), hist(intensityRange*V, 100), axis square, ...
        grid on, title('Foreground Intensity Histogram')
    saveas(gcf, 'pre_processed', 'png')
end

% Find contour coordinates of all objects (the main one and holes) and ...
% pick the first one that corresponds to the exterior boundary.
[pointStructure, objectMap, nObjects] = bwboundaries(OpenedLargestRegionMask);
% y = pointStructure{1}(:,1);
% x = pointStructure{1}(:,2);

% Thresholding to merge main object with holes.
ObjectMask = objectMap > 0;
% Erosion to initialize the snake inside the subcutaneous
% region.

end