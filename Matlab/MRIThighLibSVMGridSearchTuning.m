function [libsvmACC, libsvmgamma_n, libsvmC_n] = MRIThighLibSVMGridSearchTuning(filename, class_ids)
% syntax: [libsvmACC, libsvmgamma_n, libsvmC_n] = MRIThighLibSVMGridSearchTuning(filename, class_ids);
% % Example for libsvm training.
% filename = 'ThighMRIFeatures_20130510T141413_All_2_classes.csv';
% filename = 'ThighMRIFeatures_20130510T125956_Women_2_classes.csv';

% Read csv file.
XLS_Data = csvread(filename, 1, 1);
targets = XLS_Data(:, end-1);
XLS_Data = XLS_Data(:, 1:end-10);

% Scale features.
[scaled_patterns, pattern_means, pattern_ranges] = ScaleFeatures( XLS_Data', [], []);
I = find(targets == class_ids(1));
J = find(targets == class_ids(2));
libsvm_targets = zeros(size(targets));
libsvm_targets(I) = -1;
libsvm_targets(J) = 1;

% Set range, step, number of scales and number of repeats of grid search.
nScales = 1;
repeats = 5;
gamma_range = 24;
C_range = 24;
step = 1/2;

% Initialize counters and other variables.
count_i=0; count_j=0;
clear model, clear string, clear libsvmgamma_n, clear libsvmC_n
gamma_max = 0;
C_max = 0;
infoString = [];
tempACC = zeros(repeats, 1);

% Main search loop.
for k = 1:nScales

    % Determine step and parameter ranges.
    step = step / (2^(k-1));
    C_range = C_range / (2^(k-1));
    gamma_range = gamma_range / (2^(k-1));

    % For each experiment repeat.
    for m = 1:repeats

        % Apply random permutation.
        permuted_indices = randperm(length(libsvm_targets));
        libsvm_targets = libsvm_targets(permuted_indices);
        scaled_patterns = scaled_patterns(:, permuted_indices);
        
        %     for i = -2:2^-1:1
        for i = gamma_max-gamma_range:step:gamma_max+gamma_range
            count_i=count_i+1;
            %         for j = -2:2^-1:3
            for j = C_max-C_range:step:C_max+C_range
                count_j=count_j+1;
                libsvmgamma_n(count_i, count_j) = 2^i;
                libsvmgamma = num2str(2^i);
                libsvmC_n(count_i, count_j) = 2^j;
                libsvmC = num2str(2^j);
                string{count_i}{count_j} = ['-t 2 -e 1e-8 -v 10', ' -c ', libsvmC, ' -g ' libsvmgamma];
                libsvmACC(count_i, count_j, m) = ...
                    svmtrain(double(libsvm_targets), ...
                    double(scaled_patterns'), ...
                    string{count_i}{count_j});
            end
            count_j=0;
        end
        count_i=0;
    end
    
    % Average all repeats.
    libsvmACC = mean(libsvmACC, 3);
    
    % Find maximum ACC and arguments.
    [V, gamma_max_index] = max(libsvmACC);
    [maxACC, C_max_index] = max(V);
    C_max = libsvmC_n(gamma_max_index(C_max_index), C_max_index);
    gamma_max = libsvmgamma_n(gamma_max_index(C_max_index), C_max_index);

    % figure, surf(libsvmgamma_n, libsvmC_n, model), xlabel('libsvmgamma'), ylabel('C')
    figure(1), surf(log(libsvmgamma_n), log(libsvmC_n), libsvmACC), xlabel('libsvmgamma'), ylabel('libsvmC')
    infoString = [ infoString, sprintf('Max ACC = %d, gamma = %d, C = %d\n', maxACC, gamma_max, C_max)];
    gamma_max = log2(gamma_max);
    C_max = log2(C_max);
    
    pause( 5 );
end

% % for i = -3:0.0625:1                                                                
% % for j = 0:0.0625:3                                                                  
% '-t 2 -e 1e-8 -v 5 -c 2.9537 -g 0.52214'

fprintf(infoString);

end