function ObjectMask = ...
    GenerateObjectMaskFromContourPoints(x, y, imageSize)
% syntax: ObjectMask = ...
%     GenerateObjectMaskFromContourPoints(x, y, imageSize);


% Generate object mask.
    [x,y] = snakeinterp(x,y,1,0.5);
    objectContourMask = zeros(imageSize);
    for itPoints=1:length(x)
        objectContourMask(round(y(itPoints)), round(x(itPoints))) = 1;
    end
    objectLabels = bwlabel(1-objectContourMask, 4);
    ObjectMask = objectLabels == 2;
    
end