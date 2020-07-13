function [ObjectMask, x, y] = SnakeImplementation2D(scaledImage, ObjectMask, options, scaleFactor)
% ObjectMask = SnakeImplementation2D(scaledImage, ObjectMask, options,
% scaleFactor);


debug = options.debug;
% Gaussian smoothing sigma.
gaussianSigma = options.gaussianSigma;%0.5;
% Initial mask threshold.
threshold = options.fbthreshold;%100 / 1500;
% Canny edge threshold.
cannyThresholds =  options.cannyThresholds;

% Parzen edge detector.
parzenKernelWidth = options.parzenKernelWidth;%9;
parzenBandwidth = options.parzenBandwidth;%0.25;

% Snake algorithm.
% alpha:  elasticity parameter
% beta:   rigidity parameter
% gamma:  viscosity parameter
% kappa:  external force weight
% kappap: pressure force weight
% px,py:  external force field
nIterations = options.nIterations;%240; %110;
alpha = options.alpha;%0.05;
beta = options.beta;%1.0;
gamma = options.gamma;%1;
kappa = options.kappa;%3;
itStep = options.itStep;%5;
kappap = options.kappap;%-0.08;
snakeExternalForce = options.snakeExternalForce;
gvfInput = options.gvfInput;
gvfIterations = options.gvfIterations;
gvfMu = options.gvfMu;

% Produce parametric contour representation.
[pointStructure, objectMap, nObjects] = bwboundaries(ObjectMask);
if( nObjects > 1 )
    midPixelCoord = round(size(objectMap) / 2) + 1;
    objectLabel = objectMap( midPixelCoord(1), midPixelCoord(2) );
    y = pointStructure{objectLabel}(:,1);
    x = pointStructure{objectLabel}(:,2);
else
    y = pointStructure{1}(:,1);
    x = pointStructure{1}(:,2);
end

% Snake segmentation set-up.
%     snakeExternalForce = 'gradient';  % 'gradient', 'dtransform', 'gvf'.
%     gvfInput = 'parzen'; %'Parzen','Canny','Gradient'
infoString = sprintf('Using %s external force...\n', snakeExternalForce);
fprintf( infoString );
scaledImage0 =  gaussianBlur(scaledImage,(gaussianSigma*(2^scaleFactor))); % gaussianBlur(scaledImage,(gaussianSigma))
switch(lower(snakeExternalForce))
    case 'dtransform'
        % Distance map force.
        DT = bwdist(1-ObjectMask);
        figure,imagesc(DT),axis image, colormap gray,colorbar
        saveas(gcf, 'distance_map', 'png')
        [px,py] = gradient(-DT);
    case 'gradient'
        % Gradient-based force.
        % % note: snake potential is the negative of edge map
        %  disp(' Compute the traditional external force ...');
        [px,py] = gradient(scaledImage0);
        %             px = px .* ObjectMask;
        %             py = py .* ObjectMask;
    case 'gvf'
        % GVF approach
        infoString = sprintf('Using %s as input for GVF...\n', gvfInput);
        fprintf( infoString );
        switch(lower(gvfInput))
            case 'canny'
                % Using Canny.
                disp(' Computing Canny Edge Map for GVF ...');
                [edgeScaledImage, cannyThresholds] = edge( scaledImage0, 'canny', cannyThresholds);
                %             if( itScale == 1 )
                edgeScaledImage = edgeScaledImage .* ObjectMask;
            case 'gradient'
                % Using gradient.
                disp(' Computing Gradient Edge Map for GVF ...');
                [gx,gy] = gradient(scaledImage0);
                edgeScaledImage = sqrt(gx.*gx+gy.*gy);
                edgeScaledImage = edgeScaledImage / max(edgeScaledImage(:));
                edgeScaledImage = edgeScaledImage .* ObjectMask;
            case 'parzen'
                % Using Parzen kernels.
                disp(' Computing Parzen Edge Map for GVF ...');
                edgeScaledImage = ParzenEdgeDetection2D(scaledImage0, ...
                    parzenKernelWidth, parzenBandwidth);
                edgeScaledImage = 1 - edgeScaledImage;
                edgeScaledImage = edgeScaledImage .* ObjectMask;
        end
        %
        if( debug )
            figure,imagesc(edgeScaledImage), axis image, colormap gray, colorbar
        end
        
        % disp(' Press any key to continue with GVF and snake');
        % pause;
        % % Compute the GVF of the edge map
        disp(' Compute GVF ...');
        [u,v] = GVF(edgeScaledImage, gvfMu, gvfIterations);
        disp(' Normalizing the GVF external force ...');
        mag = sqrt(u.*u+v.*v);
        px = u./(mag+1e-10); py = v./(mag+1e-10);
        %             px = px .* ObjectMask;
        %             py = py .* ObjectMask;
end
if( debug )
    figure, quiver(imresize(px,1/8), imresize(py,1/8)), axis image
    saveas(gcf, 'snake_external_force', 'png')
end

% Apply snake deformation.
[x, y] = DeformSnake2D(x, y, px, py, options);

% Generate object mask.
ObjectMask = GenerateObjectMaskFromContourPoints(x, y, ...
    size(scaledImage));

end

