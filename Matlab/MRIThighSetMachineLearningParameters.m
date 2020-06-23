function MachineLearningParams = MRIThighSetMachineLearningParameters(varargin)
%% Set machine learning and statistical method parameters.
%     MachineLearningParams.feature_extraction_algorithm = varargin{1}{1};
%     MachineLearningParams.classifier_type = varargin{1}{2};
%     MachineLearningParams.error_method = varargin{1}{3};
%     MachineLearningParams.validation_entity = varargin{1}{4};
%     MachineLearningParams.redraws= varargin{1}{6};
%     MachineLearningParams.gammalibsvm = varargin{1}{7};
%     MachineLearningParams.Clibsvm = varargin{1}{8};
%     example: params = {'', 'libsvm','stratified-cross-validation', 'pixel', 10, 1, 2};

% Debug mode.
MachineLearningParams.debug = false;

% Image resampling factor.
MachineLearningParams.rescaleFactor = 1; % previously: 1/4;

% Feature extraction.
MachineLearningParams.feature_extraction_algorithm = '';
MachineLearningParams.outputDimensionality = 12; % previously: 14;
% {'pca','hdr','genetic_culling','information_based_selection',
% 'sequential_feature_selection','exhaustive_feature_selection'
% 'fisher_distance', 'mrmr'}

% HDR.
MachineLearningParams.correlationThreshold = 0.95; % before: 0.95, 0.975.

% Fisher.
MachineLearningParams.fisherDistanceThreshold = 0.35; % previously: 0.35.

% Classifier.
MachineLearningParams.classifier_type = 'libsvm';
% {'ls', 'em', 'svm', 'libsvm', 'ada_boost',
% 'parzen','multivariate_splines'}

% Validation.
MachineLearningParams.redraws = 10;  %previously: 2;
MachineLearningParams.percent = 90.0;   %previously: 75,7.5;
MachineLearningParams.error_method = 'resubstitution';  %previously: 'holdout';
% {'cross-validation','stratified-cross-validation', 'resubstitution','holdout','libsvmroc'}
% 'resubstitution' % uses all samples for training and then testing (train final model).
% 'holdout' % picks randomly a percentage for training, rest for testing.
% 'cross-validation' % m-fold cross-validation, iteratively uses (100/m)% of samples for testing,
% rest for training.
% 'stratified-cross-validation' preserves the ratio of samples from each
% class in each fold.

MachineLearningParams.validation_entity = 'pixel';
% {'subject', 'pixel'} :  preserve completeness of each subject in
% validation stage, or use single pixels.

% SVM.
MachineLearningParams.SVM.gamma = 5.0;

% LIBSVM
MachineLearningParams.SVM.kernellibsvm = 2;
% 0: linear, 1: polynomial, 2: RBF, 3: sigmoid, 4: user-specified

% gamma
% i = -2:0.5:6; % before: i = -2:0.125:1; i = -5:2:7;
% MachineLearningParams.SVM.gammalibsvm = 2.^i;
MachineLearningParams.SVM.gammalibsvm = 5.325e-4; % 1.414; % before: 1, 1.414, 1.189

% C
% j = -2:0.5:10; % before: j = 0:0.125:2; j = -5:2:15;
% MachineLearningParams.SVM.Clibsvm = 2.^j;
MachineLearningParams.SVM.Clibsvm = 2.2234e7; % 1.414; % before: 2.828, 2.378, 2.378

% Weighting parameters.
MachineLearningParams.SVM.w0libsvm = 1;
MachineLearningParams.SVM.w1libsvm = 1; % before: 1,5,8,11,16

% Parzen.
MachineLearningParams.Parzen.H_Normalizer = 5*10^-3; % previously: 5*10^4

% Manual parameter input.
if ~isempty(varargin)>0
    MachineLearningParams.feature_extraction_algorithm = varargin{1}{1};
    MachineLearningParams.classifier_type = varargin{1}{2};
    MachineLearningParams.error_method = varargin{1}{3};
    MachineLearningParams.validation_entity = varargin{1}{4};
    MachineLearningParams.redraws = varargin{1}{5};
    MachineLearningParams.SVM.gammalibsvm = varargin{1}{6};
    MachineLearningParams.SVM.Clibsvm = varargin{1}{7};
end

% Labels.
MachineLearningParams.ClassNames = {'Low Muscle Quality', 'High Muscle Quality'};

end
