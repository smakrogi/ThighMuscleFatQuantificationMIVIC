function ROC_Results = MRIThighPatternAnalysisCallerROC( filename, varargin )
% syntax: ROC_Results = MRIThighPatternAnalysisCallerROC( filename, classification_arguments );
% Run statistics and discriminant analysis using ROC to choose MQ threshold.
% Reads-in the csv file produced by MRIThighStatisticalTests and applies
% machine learning techniques.

%% Read first row with column headers.
fid = fopen(filename);
if( fid == -1 )
    error(['ERROR: Can not open ' filename '.']);
end
line = fgetl(fid);

%% Read subject IDs from first column.
i = 0;
while ischar(line)
    if i==0, 
        firstLine = line;
    else
        SubjectIDs{i} = strread(line, '%s', 1, 'delimiter', ',');
        SubjectIDs{i} = SubjectIDs{i}{1};
    end
    line = fgetl(fid);
    i=i+1;
end
Array_SubjectIDs = [];
for i=1:length(SubjectIDs)
    Array_SubjectIDs = strvcat(Array_SubjectIDs, SubjectIDs{i});
end

% Extract variable names.
Variable_Names = strread(firstLine, '%s', 'delimiter', ',');
Variable_Names = Variable_Names(2:end-9);
Array_Variable_Names = [];
for i=1:length(Variable_Names)
    Array_Variable_Names = strvcat(Array_Variable_Names, Variable_Names{i});
end

%% Close file.
fclose(fid);

%% Retrieve targets and patterns.
XLS_Data = csvread(filename, 1, 1);
% targets = XLS_Data(:, end-1);
MQ = XLS_Data(:, end-3);

%% Set parameters for search tange and step.
repeats = 1; % repeat experiments so we can average results in the end. Previous: 20. 
iterations = 100; % number of samples in ROC curve
min_MQ = min(MQ);
max_MQ = max(MQ);
step_MQ = (max_MQ-min_MQ)/iterations;
experiment_number = 1; 

%% Run ROC experiments.
for i=1:repeats
    for j=1:iterations
        
        %% Put all together in data structure.
        % The variable 'DataStruct' holds the patterns, targets and subject ids and
        % is produced by MRIThighStatisticalTests.m.
        % DataStruct members:
        % patterns: DxN matrix
        % targets: 1xN vector
        % subjectids: 1xN vector.
        % featurenames: Dx1 cell.
        DataStruct.patterns = XLS_Data(:, 1:end-10);
        DataStruct.patterns = DataStruct.patterns';
        threshold_MQ = min_MQ+(j-1) * step_MQ;
        targets = uint8(MQ >= threshold_MQ);
        DataStruct.targets = targets';
        DataStruct.subjectids = Array_SubjectIDs;
        DataStruct.featurenames = Array_Variable_Names;
        
        
        %% Set useful parameters.
        if(~isempty(varargin))
            MachineLearningParams = MRIThighSetMachineLearningParameters(varargin{1});
        else
            MachineLearningParams = MRIThighSetMachineLearningParameters;
        end
        
        % Run this using try/catch to handle errors coming from single class labels in
        % the test set.
        try
            
            Statistics = StatisticalPatternAnalysis( DataStruct, MachineLearningParams );
            
        catch ME
            
            % Display error and continue.
            Statistics = [];
            fprintf( 'Error description: %s, %s\n', ME.identifier, ME.message );
            
        end
        
        % Store statistical performance output to variable.
        if(~isempty(Statistics))
            ROC_Results(experiment_number).Statistics = Statistics;
            ROC_Results(experiment_number).threshold_MQ = threshold_MQ;
            ROC_Results(experiment_number).N_p = sum(targets);
            ROC_Results(experiment_number).N_n = length(targets) ...
                - ROC_Results(experiment_number).N_p;
            experiment_number = experiment_number + 1;
        end
        
        pause(1);
    end
end

%% Create ROC.
% Copy TPR, TNR and MQ values from ROC_results struct to separate variables.
for i=1:length(ROC_Results),
    TP(i) = ROC_Results(i).Statistics.ClassifValidation.TR(2);
    TN(i) = ROC_Results(i).Statistics.ClassifValidation.TR(1);
    threshold_MQ(i) = ROC_Results(i).threshold_MQ;
    ACC(i) = ROC_Results(i).Statistics.ClassifValidation.performance;
    N_p(i) = ROC_Results(i).N_p;
    N_n(i) = ROC_Results(i).N_n;    
    ACC_alt(i) = ((TP(i) * N_p(i)) + (TN(i) * N_n(i))) / (N_p(i) + N_n(i));
end

% Rearrange elements to facilitate mean and variance calculations by MQ
% threshold.
n_rows = length(ROC_Results) / repeats;
n_columns = repeats;
TP2 = reshape(TP', [n_rows, n_columns])';
TN2 = reshape(TN', [n_rows, n_columns])';
threshold_MQ2 = reshape(threshold_MQ', [n_rows, n_columns])';
ACC2 = reshape(ACC', [n_rows, n_columns])';
ACC_alt2 = reshape(ACC_alt', [n_rows, n_columns])';

% Compute means.
if n_columns==1
    TP_mean = TP2;
    TN_mean = TN2;
    threshold_MQ_mean = threshold_MQ2;
    ACC_mean = ACC2;
    ACC_alt_mean = ACC_alt2;
else
    TP_mean = mean(TP2);
    TN_mean = mean(TN2);
    threshold_MQ_mean = mean(threshold_MQ2);
    ACC_mean = mean(ACC2);
    ACC_alt_mean = mean(ACC_alt2);
end

% Sort and plot TP, FP pairs.
FP_mean = 1 - TN_mean;
[values, index] = sort(FP_mean, 'ascend');
figure, plot(FP_mean(index), TP_mean(index), 'linewidth', 2)
xlabel('False Positive Rate', 'fontsize', 12);
ylabel('True Positive Rate', 'fontsize', 12);
grid on;
for i=1:length(FP_mean),
    hold on, text(FP_mean(index(i)), TP_mean(index(i)), ...
        num2str(threshold_MQ_mean(index(i)), '%0.4g'));
end
title('ROC', 'fontsize', 12);

% Display TPR, TNR.
figure, plot(threshold_MQ_mean, TP_mean, 'linewidth', 2)
hold on, plot(threshold_MQ_mean, TN_mean, 'g', 'linewidth', 2)
legend('True Positive Rate', 'True Negative Rate', 'location', 'east');
grid on
xlabel('Muscle Quality', 'fontsize', 12)
ylabel('Likelihood', 'fontsize', 12)

% Display ACC.
figure, plot(threshold_MQ_mean, ACC_mean, 'linewidth', 2)
xlabel('Muscle Quality', 'fontsize', 12)
ylabel('Classification Accuracy', 'fontsize', 12)
grid on;

% Display ACC.
figure, plot(threshold_MQ_mean, ACC_alt_mean, 'linewidth', 2)
xlabel('Muscle Quality', 'fontsize', 12)
ylabel('Classification Accuracy', 'fontsize', 12)
grid on;

end
