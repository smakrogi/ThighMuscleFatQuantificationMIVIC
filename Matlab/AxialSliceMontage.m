function AxialSliceMontage(InputVolume, varargin)
% Make montage of the axial slices for input and output data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 2
    sliceStep = varargin{1};
    prefix = 'original_';
    GenerateFigure(InputVolume, sliceStep, prefix);
elseif nargin == 3
    SegmentationMask = varargin{1};
    sliceStep = varargin{2};
    prefix = 'original_';
    GenerateFigure(InputVolume, sliceStep, prefix);
    prefix = 'segmentation_';
    GenerateFigure(SegmentationMask, sliceStep, prefix);
end



% for j = 1:sliceStep:maskSize(3)
%     axialMask = [axialMask InputVolume(:,:,j)];
% end
% figure,imagesc(axialMask), colormap(gray), axis off,  colorbar
% saveas(gcf, [sprintf('original_montage.png')]);

end


function GenerateFigure(SegmentationMask, sliceStep, prefix)

%sliceStep = 4;
maskSize = size(SegmentationMask);
axialMask = [];
% borderLength = 5;

% BoundingBox = FindAxialBoundingBox(SegmentationMask);
% if(BoundingBox.xmin > borderLength) BoundingBox.xmin = ...
%         BoundingBox.xmin - borderLength; end
% if(BoundingBox.xmax < maskSize(1) - borderLength) BoundingBox.xmax = ...
%         BoundingBox.xmax + borderLength; end
% if(BoundingBox.ymin > borderLength) BoundingBox.ymin = ...
%         BoundingBox.ymin - borderLength; end
% if(BoundingBox.ymax < maskSize(2) - borderLength) BoundingBox.ymax = ...
%         BoundingBox.ymax + borderLength; end
% 
% InputVolume = InputVolume(BoundingBox.xmin:BoundingBox.xmax, ...
%     BoundingBox.ymin:BoundingBox.ymax, :);
% SegmentationMask = SegmentationMask(BoundingBox.xmin:BoundingBox.xmax, ...
%     BoundingBox.ymin:BoundingBox.ymax, :, :);
       
if numel(maskSize) == 4
    nClusters = maskSize(4);
else
    nClusters = 1;
end

for i = 1:nClusters
    axialMask = [];
    for j = 1:sliceStep:maskSize(3)
        axialMask = [axialMask SegmentationMask(:,:,j,i)];
    end
    figure,imagesc(axialMask), colormap(gray), axis image, axis off,  colorbar
    saveas(gcf, [prefix sprintf('_montage_%d.png', i)]);
end

end



