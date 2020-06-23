function options = SetPrinceSnakeOptionsThigh(varargin)
% Set snake algorithm options.

if(isempty(varargin))
    % Debug flag.
    options.debug = 0;    
    options.leg = 'left';
    options.nScales = 3; %3;    
    % Initial mask threshold for leg selection.
    options.fbthreshold = 300 / 1500;  %300/1500
    % Structural element size for definition of initial mask.
    options.seSize = 4;  % 4, 10
    % Gaussian smoothing sigma.
    options.gaussianSigma = 1.0;
    % Top hat parameters.
    options.useTopHat = 0;
    options.strelSizeFactor = 16;
    options.snakeExternalForce =  'gradient';  % 'gradient', 'dtransform', 'parzen', 'gvf';
    options.gvfInput = 'Canny'; %'Parzen','Canny','Gradient';
else
    options.debug = varargin{1}.debug;    
    options.leg = varargin{1}.ProcessedLeg;
    options.nScales = varargin{1}.nScales{varargin{2}};
    options.fbthreshold =  varargin{1}.LegSelectionThreshold{varargin{2}};
    options.seSize = varargin{1}.snakemaskseSize;  % 4, 10
    options.gaussianSigma = varargin{1}.snakeGaussianSigma;
    options.useTopHat = varargin{1}.tophatTransform;
    options.strelSizeFactor = varargin{1}.strelSizeFactorTopHat;
    options.snakeExternalForce = varargin{1}.snakeExternalForce;
    options.gvfInput = varargin{1}.gvfInput;
end

% Canny edge threshold.
options.cannyFactor = 0.15;
options.cannyThresholds = [0.03 0.3];

% Parzen edge detector.
options.parzenKernelWidth = 3;
options.parzenBandwidth = 0.1;

% Snake algorithm.
% alpha:  elasticity parameter
% beta:   rigidity parameter
% gamma:  viscosity parameter
% kappa:  external force weight
% kappap: pressure force weight
% px,py:  external force field
options.nIterations = 220; %110, 160; 
options.alpha = 0.05; %0.05; % greater values,  harder to stretch
options.beta = 1;  %1.0; % greater values, more rigid deformation
options.gamma = 1; %1;  % greater values, slower deformation but more viscous.
options.kappa = 2.0; %1.0;  % greater values, stronger force around the edges.
options.itStep = 5;
options.kappap = -0.05;
options.gvfIterations = 50;
options.gvfMu = 0.2;
