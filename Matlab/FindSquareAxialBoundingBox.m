function [BoundingBox] = FindSquareAxialBoundingBox(ReferenceVolume)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[BoundingBox] = FindAxialBoundingBox(ReferenceVolume);

% Make the bounding box square in the axial plane
% in order to run matlab matrix ops.
height = BoundingBox.xmax - BoundingBox.xmin;
width = BoundingBox.ymax - BoundingBox.ymin;

squareEdge  = max(height, width);

BoundingBox.xmax = BoundingBox.xmin + squareEdge;
BoundingBox.ymax = BoundingBox.ymin + squareEdge;

BoundingBox.xmax = min( BoundingBox.xmax, size(ReferenceVolume,1));
BoundingBox.ymax = min( BoundingBox.ymax, size(ReferenceVolume,2));

end