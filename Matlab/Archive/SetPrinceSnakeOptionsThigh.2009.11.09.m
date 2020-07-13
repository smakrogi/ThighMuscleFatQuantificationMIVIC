function options = SetPrinceSnakeOptionsThigh(varargin)

% Multi-scale scheme.
options.debug = 1;
options.nScales = 3; %2;

% Top hat parameters.
options.useTopHat = 0;
options.strelSizeFactor = 8;

% Gaussian smoothing sigma.
options.gaussianSigma = 0.5;

% Initial mask threshold.
options.threshold = 100 / 1500;

% Canny edge threshold.
options.cannyFactor = 0.15;
options.seSize = 10;

% Parzen edge detector.
options.parzenKernelWidth = 9;
options.parzenBandwidth = 0.25;

% Snake algorithm.
% alpha:  elasticity parameter
% beta:   rigidity parameter
% gamma:  viscosity parameter
% kappa:  external force weight
% kappap: pressure force weight
% px,py:  external force field
options.nIterations = 160; %110; 
options.alpha = 0.05;
options.beta = 1.0;
options.gamma = 1;
options.kappa = 2.75;
options.itStep = 5;
options.kappap = -0.05;
