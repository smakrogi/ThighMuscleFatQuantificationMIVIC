function [ObjectMask, x, y] = SegmentInternalSubcutaneousSurface2D(originalImage, options)
% Delineation of internal subutaneous adipose surface.

% Multi-scale scheme.
debug = options.debug; %1;
nScales = options.nScales; %3; %2;

% Top hat parameters.
% useTopHat = options.useTopHat; %0;
% strelSizeFactor = options.strelSizeFactor; %8;

% Morphological params.
seSize = options.seSize;%10;

% Initial mask threshold.
threshold = options.fbthreshold; %100 / 1500;

% Snake algorithm.
nIterations = options.nIterations; %240; %110; 
itStep = options.itStep; %5;

originalImage = double(originalImage);
nIterationsOriginal = nIterations;

% Hierarchical scheme.
for itScale=1:nScales
    
    % Generate new scale.
    scaleFactor = itScale-nScales;
    scaledImage = imresize( originalImage, 2^scaleFactor );
    nIterations = round( nIterationsOriginal / 2^(itScale-1) );
    
    intensityRange = max(scaledImage(:))-min(scaledImage(:));
    scaledImage = (scaledImage-min(scaledImage(:))) / intensityRange;
    
    se = strel('disk', round(seSize * 2^scaleFactor));
    
    if( debug )
        figure(1),imagesc(scaledImage*intensityRange),axis image, colormap gray,colorbar, ...
            title('Scaled Image')
        saveas(gcf, 'snake_input', 'png')
    end
    
    % Define subcutaneous adipose initial contour
    % through morphological operations.
    if( itScale == 1 )
        
        % Thresholding.
        [ObjectMask] = ...
            MaskObjectRegionBeforeSnake(scaledImage, ...
            intensityRange, ...
            0, ...
            debug, ...
            threshold);
        ObjectMask = imerode(ObjectMask, se);
    else
        ObjectMask = imresize( ObjectMask,  size(scaledImage));
        ObjectMask = ObjectMask > 0;
        ObjectMask = imdilate(ObjectMask, se);
    end
    
    % Generate final object mask using snakes.
    % Implement snake force computation and evolution.
    [ObjectMask, x, y] = SnakeImplementation2D(scaledImage, ...
        ObjectMask, ...
        options, ...
        scaleFactor);
    
end

% Display final segmentation.
if( debug )
    figure(7), imagesc(originalImage); axis('image', 'off'); colormap(gray)
    snakedisp(x,y,'r');
    title(['Final result,  iter = ' num2str(nIterations*itStep)]);
    saveas(gcf, 'snake_final', 'png')
end

% Change datatype to save memory.
ObjectMask = uint8(ObjectMask);

end

