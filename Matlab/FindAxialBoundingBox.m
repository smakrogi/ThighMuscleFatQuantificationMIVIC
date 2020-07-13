function [BoundingBox] = FindAxialBoundingBox(ReferenceVolume)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frame length in pixels.
frameLength = 8; % 5

% Crop matrix size.
BoundingBox.xmin = size(ReferenceVolume,1);
BoundingBox.xmax = 1;
BoundingBox.ymin = size(ReferenceVolume,2);
BoundingBox.ymax = 1;

% Compute bounding box.
for i=1:size(ReferenceVolume,3)
    [I, J] = find(ReferenceVolume(:,:,i));
    BoundingBox.xmin = min(BoundingBox.xmin, min(I));
    BoundingBox.xmax = max(BoundingBox.xmax, max(I));
    BoundingBox.ymin = min(BoundingBox.ymin, min(J));
    BoundingBox.ymax = max(BoundingBox.ymax, max(J));
end


% Generate a frame around the ROI.
if(BoundingBox.xmin > frameLength)
    BoundingBox.xmin = BoundingBox.xmin - frameLength;
end
if(BoundingBox.xmax <= size(ReferenceVolume,1) - frameLength)
    BoundingBox.xmax = BoundingBox.xmax + frameLength;
end
if(BoundingBox.ymin > frameLength)
    BoundingBox.ymin = BoundingBox.ymin - frameLength;
end
if(BoundingBox.ymax <= size(ReferenceVolume,2) - frameLength)
    BoundingBox.ymax = BoundingBox.ymax + frameLength;
end


% % Crop matrix size.
% BoundingBox.xmin = size(ReferenceVolume,1);
% BoundingBox.xmax = 1;
% BoundingBox.ymin = size(ReferenceVolume,2);
% BoundingBox.ymax = 1;
% 
% for i=1:size(ReferenceVolume,3)
%     [I, J] = find(ReferenceVolume(:,:,i));
%     BoundingBox.xmin = min(BoundingBox.xmin, min(I));
%     BoundingBox.xmax = max(BoundingBox.xmax, max(I));
%     BoundingBox.ymin = min(BoundingBox.ymin, min(J));
%     BoundingBox.ymax = max(BoundingBox.ymax, max(J));
% end

end