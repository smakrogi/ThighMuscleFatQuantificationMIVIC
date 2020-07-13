function Stats = BLSA_Thigh_MRI_Statistics(DataTable, Ind_GroupA, Ind_GroupB)

% Add required paths.
addpath('teg_repeated_measures_ANOVA');
addpath('mixed_between_within_anova');

% Set number of features and number of groups.
nFeatures = size(DataTable, 2);
% nGroups = 3;

% Choose samples using criterion

% Label rows according to group number.
GroupA = DataTable(Ind_GroupA, :);
GroupB = DataTable(Ind_GroupB, :);

% Calculate statistics within each group and between groups.

% Iterative display barplots organized by group and dates and
% compute p-values.
count = 1;
for i=1:nFeatures
    % Measurements.
    DataGroupA = GroupA(:,i).Variables;
    DataGroupB = GroupB(:,i).Variables;
    FeatureNames{i} = DataTable(:,i).Properties.VariableNames{1};
    
    if isnumeric(DataGroupA) || islogical(DataGroupA)
        
        % Expected values.
        Stats(count).Name = FeatureNames{i};
        Stats(count).MeanA = mean(DataGroupA);
        Stats(count).StDevA = std(DataGroupA);
        
        Stats(count).MeanB = mean(DataGroupB);
        Stats(count).StDevB = std(DataGroupB);
        
        % Two-sample t-tests for cross-sectional stats.
        [Stats(count).H_val, Stats(count).P_val] = ...
            ttest2(DataGroupA, DataGroupB, 'Vartype', 'unequal');
        
        % Sample size power.
        if (Stats(count).StDevA > 0 && Stats(count).StDevB > 0 && Stats(count).MeanA ~= Stats(count).MeanB)
            Stats(count).Nout = sampsizepwr('t2',[Stats(count).MeanA, Stats(count).StDevA], ...
                Stats(count).MeanB);
            Stats(count).Power = sampsizepwr('t2',[Stats(count).MeanA, Stats(count).StDevA], ...
                Stats(count).MeanB, [], 5);
        else
            Stats(count).Nout = NaN;
            Stats(count).Power = NaN;
        end
        
        % Barplots.
        InfoString2 = sprintf('p: %.4f, h: %d', ...
            Stats(count).P_val, Stats(count).H_val);
        BarValues = [mean(DataGroupA), mean(DataGroupB)];
        ErrorValues = [std(DataGroupA), std(DataGroupB)];
        figure, barwitherr(ErrorValues, BarValues); colormap winter
        title(FeatureNames{i}, 'interpreter', 'none', 'fontsize', 18);
        set(gca,'XTickLabel',{'Group A','Group B'}, 'fontweight','b', 'fontsize', 18);
        text(double(0.5), max(BarValues(:))+max(ErrorValues(:)), ...
            InfoString2, ...
            'fontweight','b', 'fontsize', 18);
        axis normal;
        
        % Write plot to png image file.
        axis normal;
        saveas(gcf, [FeatureNames{i}, '_BarCharts.png']);
        
        
        count = count + 1;
    end
end


% Write longitudinal t-tests to csv file, pay attention to row order and their correspondence to groups.
Table = struct2table(Stats);
writetable(Table, 'BLSA_t_tests.csv');
disp(Table);

end