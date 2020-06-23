function DataStruct = MRIThighStatisticalTests(filename, selected_gender, n_classes)
% syntax: DataStruct = MRIThighStatisticalTests(filename, selected_gender, n_classes);
% selected_gender: 0 for women or 1 for men
% n_classes: 2 or 3.
% This function reads-in the feature matrix and related class information from the same input csv file,
% computes area fractions and intensity ratios, applies eigen-decomposition
% to tissues and computes statistical tests.
% It then writes the new data matrix to a csv file and the p-values to
% another csv file.


% Read excel file.


% Read first row with column headers.
fid = fopen(filename);
if( fid == -1 )
    error(['ERROR: Can not open ' filename '.']);
end
line = fgetl(fid);


% Read subject IDs from first column.
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


% Extract variable names.
Variable_Names = strread(firstLine, '%s', 'delimiter', ',');
Variable_Names = Variable_Names(5:end);
variable_names_length = length(Variable_Names);


% Close file.
fclose(fid);


% Read numerical data.
XLS_Data = csvread(filename, 1, 4);
[N, D] = size(XLS_Data);


% Retrieve BLSA fields.
% Diabetes_Mellitus = XLS_Data(1:N ,D);
% Diabetes_Mellitus_String = Variable_Names{variable_names_length};
Muscle_Quality = XLS_Data(1:N ,D-1);
BMI = XLS_Data(1:N ,D - 2);
Age = XLS_Data(1:N ,D - 4);
Gender = XLS_Data(1:N ,D - 5);


% Separate older BLSA data and variable names from the MRI ones.
BLSA_Data = XLS_Data(:,D-7:D);
BLSA_Variable_Names = Variable_Names(variable_names_length-7:variable_names_length);


% Remove BLSA field names and values from original tables..
XLS_Data = XLS_Data(:,1:D-8);
Variable_Names = Variable_Names(1:variable_names_length-8);
[N, D] = size(XLS_Data);


% Select data according to gender selection (0 for women or 1 for men).
if ~isempty(selected_gender)
    Indices = find(Gender==selected_gender);
    XLS_Data = XLS_Data(Indices,:);
    [N, D] = size(XLS_Data);
%     Diabetes_Mellitus = Diabetes_Mellitus(Indices,:);
    Muscle_Quality = Muscle_Quality(Indices,:);
    BMI = BMI(Indices,:);
    Age = Age(Indices,:);
    Gender = Gender(Indices,:);
    SubjectIDs = SubjectIDs(Indices);
    BLSA_Data = BLSA_Data(Indices,:);
end


% Determine classes for statistics and boxplots.
Muscle_Quality_Class = ...
    DetermineMuscleQualityClasses(Muscle_Quality, n_classes);


% Compute area fractions.
[XLS_Data, Variable_Names] = ComputeAreaFractionFeatures(XLS_Data, ...
    Variable_Names);


% Compute intensity fractions.
[XLS_Data, Variable_Names] = ComputeIntensityBasedFeatures(XLS_Data, ...
    Variable_Names);


% Co-variance matrix eigen-decomposition.
if D>20
    [XLS_Data, Variable_Names] = ComputeCovarianceMatrixFeatures(XLS_Data, ...
        Variable_Names);
end


% Add the BLSA dB fields to the end.
XLS_Data = [XLS_Data, BLSA_Data];
Variable_Names = [Variable_Names; BLSA_Variable_Names];

% Group-wise statistical tests.
GroupStats = ComputeStatisticalTestsAndDisplayPlots(XLS_Data, ...
    Muscle_Quality_Class, ...
    Variable_Names);

% Create data structure for further analysis.
DataStruct.GroupStats = GroupStats;
DataStruct.DataMatrix = XLS_Data;
DataStruct.VariableNames = Variable_Names;
DataStruct.SubjectIDs = SubjectIDs;
DataStruct.Gender = Gender;
DataStruct.Age = Age;
DataStruct.MuscleQualityClass = Muscle_Quality_Class;

% Write stats to csv file.
WriteGroupStatistics(GroupStats, Variable_Names);
timestamp = WriteFeatures(XLS_Data, Variable_Names, SubjectIDs, Muscle_Quality_Class);

% Write data to mat file.
save(['ThighMRIFeatures_', timestamp, '.mat'], 'DataStruct');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Muscle_Quality_Class = DetermineMuscleQualityClasses(Muscle_Quality, nClasses)

if nClasses==2
    threshold = 63.1; % M:63, W:62
    Muscle_Quality_Class = Muscle_Quality > threshold;
elseif nClasses==3
    threshold1 = 53; % M:54, W:50
    threshold2 = 72.7; % M:75, W:72
    Muscle_Quality_Class = (Muscle_Quality > threshold1) + ...
        (Muscle_Quality > threshold2);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timestamp = WriteFeatures(XLS_Data, Variable_Names, SubjectIDs, Muscle_Quality_Class)

timestamp = datestr(now, 30);
fid=fopen(['ThighMRIFeatures_', timestamp, '.csv'], 'wb', 'l');

% Write subject IDs.
fprintf(fid, '%s, ', 'Subject ID');

% Write other fields.
for i=1:length(Variable_Names)
    fprintf(fid, '%s, ', Variable_Names{i});
end
fprintf(fid, '%s, ', 'MQ Class');
fprintf(fid, '\n');

for i=1:size(XLS_Data, 1)
        fprintf(fid, '%s, ', SubjectIDs{i});
    for j=1:size(XLS_Data, 2)
        fprintf(fid, '%.6f, ', XLS_Data(i,j));
    end
    fprintf(fid, '%.6f, ', Muscle_Quality_Class(i));
    fprintf(fid, '\n');
end

fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteGroupStatistics(GroupStats, Variable_Names)

timestamp = datestr(now, 30);
fid=fopen(['ThighMRIGroupStatistics_', timestamp, '.csv'], 'wb', 'l');

nClasses = size(GroupStats.ttest2, 2);
if nClasses == 2
    fprintf(fid, ' , t-test2, , ks-test2, , fisher distance, , mean, mean, , stdev, stdev\n');
    fprintf(fid, ' , G0~G1, , G0~G1, , G0~G1, , G0, G1, , G0, G1\n');
elseif nClasses == 3
    fprintf(fid, ' , t-test2, , , ks-test2, , , fisher distance, , , mean, mean, mean, , stdev, stdev, stdev\n');
    fprintf(fid, ' , G0~G1, G0~G2, G1~G2, G0~G1, G0~G2, G1~G2, G0~G1, G0~G2, G1~G2, G0, G1, G2, , G0, G1, G2\n');
end

for i=1:length(Variable_Names)
    if nClasses == 2
        fprintf(fid, '%s, %3.6f, , ', Variable_Names{i}, ...
            GroupStats.ttest2(i,1));
        fprintf(fid, '%3.6f, , ', GroupStats.kstest2(i,1));
        fprintf(fid, '%3.6f, , ', GroupStats.fisherratio(i,1));
        fprintf(fid, '%3.6f, %3.6f, ,', GroupStats.mean(i,:));
        fprintf(fid, '%3.6f, %3.6f\n', GroupStats.stdev(i,:));
    elseif nClasses == 3
        fprintf(fid, '%s, %3.6f, %3.6f, %3.6f, ', Variable_Names{i}, ...
            GroupStats.ttest2(i,:));
        fprintf(fid, '%3.6f, %3.6f, %3.6f, ', GroupStats.kstest2(i,:));
        fprintf(fid, '%3.6f, %3.6f, %3.6f, ', GroupStats.fisherratio(i,:));
        fprintf(fid, '%3.6f, %3.6f, %3.6f, ', GroupStats.mean(i,:));
        fprintf(fid, '%3.6f, %3.6f, %3.6f\n', GroupStats.stdev(i,:));
    end
end

fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Statistics = ComputeStatisticalTestsAndDisplayPlots(DataMatrix, ...
    Targets, ...
    Variable_Names)

[N, D] = size(DataMatrix);

nClasses = length(unique(Targets));
index0 = find(Targets==0);
index1 = find(Targets==1);
index2 = find(Targets==2);

Statistics.ttest2 = zeros(D, length(unique(Targets)));
Statistics.kstest2 = zeros(D, length(unique(Targets)));
Statistics.fisherratio = zeros(D, length(unique(Targets)));

for i=1:D
    if ~isnan(DataMatrix(1:N,i))

        % Create filename string without dots and spaces.
        filename_string = strrep(Variable_Names{i}, '.', '');
        filename_string = strrep( filename_string, ' ', '');
        filename_string = strrep( filename_string, '\n', '');
        
        % Means and stdevs.
        Statistics.mean(i, 1) = mean(DataMatrix(index0, i));
        Statistics.mean(i, 2) = mean(DataMatrix(index1, i));
        if nClasses==3
            Statistics.Mean(i, 3) = mean(DataMatrix(index2, i));
        end
        
        Statistics.stdev(i, 1) = std(DataMatrix(index0, i));
        Statistics.stdev(i, 2) = std(DataMatrix(index1, i));
        if nClasses==3
            Statistics.stdev(i, 3) = std(DataMatrix(index2, i));
        end
        
        % Two-sample Student's t tests
        [H, Statistics.ttest2(i, 1)] = ttest2(DataMatrix(index0, i), ...
            DataMatrix(index1, i));
        if nClasses==3
            [H, Statistics.ttest2(i, 2)] = ttest2(DataMatrix(index0, i), ...
                DataMatrix(index2, i));
            [H, Statistics.ttest2(i, 3)] = ttest2(DataMatrix(index1, i), ...
                DataMatrix(index2, i));
        end

        %     Statistics.rho(i,:,:) = corr(DataStruct.patterns(i,index0), ...
        %         DataStruct.patterns(i,index1), 'type', 'Spearman');

        % Two-sample Kolmogorov-Smirnov test
        [H, Statistics.kstest2(i, 1)] = kstest2(DataMatrix(index0, i), ...
            DataMatrix(index1, i));
        if nClasses==3
            [H, Statistics.kstest2(i, 2)] = kstest2(DataMatrix(index0, i), ...
                DataMatrix(index2, i));
            [H, Statistics.kstest2(i, 3)] = kstest2(DataMatrix(index1, i), ...
                DataMatrix(index2, i));
        end

        % Fisher ratio.
        Statistics.fisherratio(i, 1) = abs(mean(DataMatrix(index0, i)) - ...
            mean(DataMatrix(index1, i))) ...
            / (std(DataMatrix(index0, i)) + ...
            std(DataMatrix(index1, i)));
        if nClasses==3
            Statistics.fisherratio(i, 2) = abs(mean(DataMatrix(index0, i)) - ...
                mean(DataMatrix(index2, i))) ...
                / (std(DataMatrix(index0, i)) + ...
                std(DataMatrix(index2, i)));
            Statistics.fisherratio(i, 3) = abs(mean(DataMatrix(index1, i)) - ...
                mean(DataMatrix(index2, i))) ...
                / (std(DataMatrix(index1, i)) + ...
                std(DataMatrix(index2, i)));
        end
        
        % Build histograms.
        if nClasses==2
                [n0,x0] = hist(DataMatrix(index0, i), 10);
                [n1,x1] = hist(DataMatrix(index1, i), 10);
                figure(1)
                bar(x0, n0, 'r')
                hold on
                bar(x1, n1, 'g')
                title([Variable_Names{i}, ' Histogram'], 'FontSize', 18)
                legend('G0', 'G1')
                hold off
        elseif nClasses==3
                [n0,x0] = hist(DataMatrix(index0, i), 10);
                [n1,x1] = hist(DataMatrix(index1, i), 10);
                [n2,x2] = hist(DataMatrix(index2, i), 10);                
                figure(1)
                bar(x0, n0, 'r')
                hold on
                bar(x1, n1, 'g')
                hold on
                bar(x2, n2, 'b')
                title(['Histogram of ', Variable_Names{i}], 'FontSize', 18)
                legend('G0', 'G1', 'G2')
                hold off
        end
        saveas(gcf, [filename_string, 'Hist.png']);
        
        % Display box plots for each variable grouped by quality.
        % Save figures in png format.
        figure(2), boxplot(DataMatrix(1:N,i)', Targets, 'whisker', 10), ...
            title(Variable_Names{i}, 'FontSize', 18),
        saveas(gcf, [filename_string, 'Box.png']);
        
                
                % Scatter plots vs. Age
                % figure, plot(Age, SAT_AreaFraction, '.'), ...
                %     title('SAT Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'SATAreaFractionvsAge.png');
                % figure, plot(Age, MuscleAreaFraction, '.'), ...
                %     title('Muscle Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'MuscleAreaFractionvsAge.png');
                % figure, plot(Age, InterMFAreaFraction, '.'), ...
                %     title('InterMF Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'InterMF Area FractionvsAge.png');
                % figure, plot(Age, IntraMFAreaFraction, '.'), ...
                %     title('IntraMF Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'IntraMFAreaFractionvsAge.png');
                % figure, plot(Age, InterMFAreaFraction+IntraMFAreaFraction, '.'), ...
                %     title('InterMF+IntraMF Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'InterIntraMFAreaFractionvsAge.png');

    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NewXLS_Data, NewVariable_Names] = ComputeAreaFractionFeatures(XLS_Data, ...
    Variable_Names)
% Compute area fractions and display box plots.
% Save figures in png format.
[N, D] = size(XLS_Data);

Muscle_Area = XLS_Data(1:N,1);
SAT_Area = XLS_Data(1:N,2);
InterMF_Area = XLS_Data(1:N,3);
IntraMF_Area = XLS_Data(1:N,4);
Bone_Area = XLS_Data(1:N,5);

SAT_AreaFraction = SAT_Area ./ (Muscle_Area + InterMF_Area + IntraMF_Area + SAT_Area);
% MuscleAreaFraction = Muscle_Area ./ (Muscle_Area + InterMF_Area + IntraMF_Area + SAT_Area);
% InterMFAreaFraction = InterMF_Area ./ (Muscle_Area + InterMF_Area + IntraMF_Area + SAT_Area);
% IntraMFAreaFraction = IntraMF_Area ./ (Muscle_Area + InterMF_Area + IntraMF_Area + SAT_Area);
MuscleAreaFraction = Muscle_Area ./ ( Muscle_Area + InterMF_Area + IntraMF_Area );
if sum(IntraMF_Area) == 0
    InterMFAreaFraction = InterMF_Area ./ ( Muscle_Area + InterMF_Area + IntraMF_Area );
    IntraMFAreaFraction = IntraMF_Area ./ ( Muscle_Area + InterMF_Area + IntraMF_Area );
else
    InterMFAreaFraction = InterMF_Area ./ ( InterMF_Area + IntraMF_Area );
    IntraMFAreaFraction = IntraMF_Area ./ ( InterMF_Area + IntraMF_Area );
end

XLS_Data = [XLS_Data, SAT_AreaFraction, MuscleAreaFraction, ...
    InterMFAreaFraction, IntraMFAreaFraction];
Variable_Names = [Variable_Names; 'SAT Area Fraction'; 'Muscle Area Fraction'; ...
    'InterMF Area Fraction'; 'IntraMF Area Fraction'];

% % Box plots
% figure, boxplot(SAT_AreaFraction, Muscle_Quality_Class, 'whisker', 10), ...
%     title('SAT Area Fraction', 'FontSize', 18), saveas(gcf, 'SATAreaFraction.png');
% figure, boxplot(MuscleAreaFraction, Muscle_Quality_Class, 'whisker', 10), ...
%     title('Muscle Area Fraction', 'FontSize', 18), saveas(gcf, 'MuscleAreaFraction.png');
% figure, boxplot(InterMFAreaFraction, Muscle_Quality_Class, 'whisker', 10), ...
%     title('InterMF Area Fraction', 'FontSize', 18), saveas(gcf, 'InterMFAreaFraction.png');
% figure, boxplot(IntraMFAreaFraction, Muscle_Quality_Class, 'whisker', 10), ...
%     title('IntraMF Area Fraction', 'FontSize', 18), saveas(gcf, 'IntraMFAreaFraction.png');
% figure, boxplot(InterMFAreaFraction+IntraMFAreaFraction, Muscle_Quality_Class, ...
%     'whisker', 10), title('InterMF+IntraMF Area Fraction', 'FontSize', 18), ...
%     saveas(gcf, 'InterIntraMFAreaFraction.png');

% Scatter plots vs. Age
% figure, plot(Age, SAT_AreaFraction, '.'), ...
%     title('SAT Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'SATAreaFractionvsAge.png');
% figure, plot(Age, MuscleAreaFraction, '.'), ...
%     title('Muscle Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'MuscleAreaFractionvsAge.png');
% figure, plot(Age, InterMFAreaFraction, '.'), ...
%     title('InterMF Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'InterMF Area FractionvsAge.png');
% figure, plot(Age, IntraMFAreaFraction, '.'), ...
%     title('IntraMF Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'IntraMFAreaFractionvsAge.png');
% figure, plot(Age, InterMFAreaFraction+IntraMFAreaFraction, '.'), ...
%     title('InterMF+IntraMF Area Fraction vs. Age', 'FontSize', 18), saveas(gcf, 'InterIntraMFAreaFractionvsAge.png');

NewXLS_Data = XLS_Data;
NewVariable_Names = Variable_Names;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NewXLS_Data, NewVariable_Names] = ComputeIntensityBasedFeatures(XLS_Data, ...
    Variable_Names)

% Compute intensity fractions and display box plots.
% Save figures in png format.
[N, D] = size(XLS_Data);

% NS
Muscle_NS_Intensity = XLS_Data(1:N,6);
SAT_NS_Intensity = XLS_Data(1:N,7);
InterMF_NS_Intensity = XLS_Data(1:N,8);
IntraMF_NS_Intensity = XLS_Data(1:N,9);

ScaledIntensityFeatures.ScaledMuscleNSIntensity = Muscle_NS_Intensity ./ SAT_NS_Intensity;
ScaledIntensityFeatures.ScaledInterMFNSIntensity = InterMF_NS_Intensity ./ SAT_NS_Intensity;
ScaledIntensityFeatures.ScaledIntraMFNSIntensity = IntraMF_NS_Intensity ./ SAT_NS_Intensity;
% ScaledIntensityFeatures.ScaledMuscleNSIntensity = Muscle_NS_Intensity ./ (SAT_NS_Intensity + Muscle_NS_Intensity);
% ScaledIntensityFeatures.ScaledInterMFNSIntensity = InterMF_NS_Intensity ./ (SAT_NS_Intensity + Muscle_NS_Intensity);
% ScaledIntensityFeatures.ScaledIntraMFNSIntensity = IntraMF_NS_Intensity ./ (SAT_NS_Intensity + Muscle_NS_Intensity);

% figure, boxplot(ScaledIntensityFeatures.ScaledMuscleNSIntensity, Muscle_Quality, 'whisker', 10), ...
%     title('Muscle NS Intensity Fraction', 'FontSize', 18), saveas(gcf, 'ScaledMuscleNSIntensity.png');
% figure, boxplot(ScaledIntensityFeatures.ScaledInterMFNSIntensity, Muscle_Quality, 'whisker', 10), ...
%     title('InterMF NS Intensity Fraction', 'FontSize', 18), saveas(gcf, 'ScaledInterMFNSIntensity.png');
% if(~isnan(ScaledIntensityFeatures.ScaledIntraMFNSIntensity)), figure, boxplot(ScaledIntensityFeatures.ScaledIntraMFNSIntensity, Muscle_Quality, 'whisker', 10), ...
%         title('IntraMF NS Intensity Fraction', 'FontSize', 18), saveas(gcf, 'ScaledIntraMFNSIntensity.png'); end

% FS
Muscle_FS_Intensity = XLS_Data(1:N,11);
SAT_FS_Intensity = XLS_Data(1:N,12);
InterMF_FS_Intensity = XLS_Data(1:N,13);
IntraMF_FS_Intensity = XLS_Data(1:N,14);

ScaledIntensityFeatures.ScaledMuscleFSIntensity = Muscle_FS_Intensity ./ SAT_FS_Intensity;
ScaledIntensityFeatures.ScaledInterMFFSIntensity = InterMF_FS_Intensity ./ SAT_FS_Intensity;
ScaledIntensityFeatures.ScaledIntraMFFSIntensity = IntraMF_FS_Intensity ./ SAT_FS_Intensity;
% ScaledIntensityFeatures.ScaledMuscleFSIntensity = Muscle_FS_Intensity ./ (SAT_FS_Intensity + Muscle_FS_Intensity);
% ScaledIntensityFeatures.ScaledInterMFFSIntensity = InterMF_FS_Intensity ./ (SAT_FS_Intensity + Muscle_FS_Intensity);
% ScaledIntensityFeatures.ScaledIntraMFFSIntensity = IntraMF_FS_Intensity ./ (SAT_FS_Intensity + Muscle_FS_Intensity);

% figure, boxplot(ScaledIntensityFeatures.ScaledMuscleFSIntensity, Muscle_Quality, 'whisker', 10), ...
%     title('Muscle FS Intensity Fraction', 'FontSize', 18), saveas(gcf, 'ScaledMuscleFSIntensity.png');
% figure, boxplot(ScaledIntensityFeatures.ScaledInterMFFSIntensity, Muscle_Quality, 'whisker', 10), ...
%     title('InterMF FS Intensity Fraction', 'FontSize', 18), saveas(gcf, 'ScaledInterMFFSIntensity.png');
% if(~isnan(ScaledIntensityFeatures.ScaledIntraMFFSIntensity)), figure, boxplot(ScaledIntensityFeatures.ScaledIntraMFFSIntensity, Muscle_Quality, 'whisker', 10), ...
%         title('IntraMF FS Intensity Fraction', 'FontSize', 18), saveas(gcf, 'ScaledIntraMFFSIntensity.png'); end

% WS
Muscle_WS_Intensity = XLS_Data(1:N,16);
SAT_WS_Intensity = XLS_Data(1:N,17);
InterMF_WS_Intensity = XLS_Data(1:N,18);
IntraMF_WS_Intensity = XLS_Data(1:N,19);

ScaledIntensityFeatures.ScaledMuscleWSIntensity = Muscle_WS_Intensity ./ SAT_WS_Intensity;
ScaledIntensityFeatures.ScaledInterMFWSIntensity = InterMF_WS_Intensity ./ SAT_WS_Intensity;
ScaledIntensityFeatures.ScaledIntraMFWSIntensity = IntraMF_WS_Intensity ./ SAT_WS_Intensity;
% ScaledIntensityFeatures.ScaledMuscleWSIntensity = Muscle_WS_Intensity ./ (SAT_WS_Intensity + Muscle_WS_Intensity);
% ScaledIntensityFeatures.ScaledInterMFWSIntensity = InterMF_WS_Intensity ./ (SAT_WS_Intensity + Muscle_WS_Intensity);
% ScaledIntensityFeatures.ScaledIntraMFWSIntensity = IntraMF_WS_Intensity ./ (SAT_WS_Intensity + Muscle_WS_Intensity);

% figure, boxplot(ScaledIntensityFeatures.ScaledMuscleWSIntensity, Muscle_Quality, 'whisker', 10), ...
%     title('Muscle WS Intensity Fraction', 'FontSize', 18), saveas(gcf, 'ScaledMuscleWSIntensity.png');
% figure, boxplot(ScaledIntensityFeatures.ScaledInterMFWSIntensity, Muscle_Quality, 'whisker', 10), ...
%     title('InterMF WS Intensity Fraction', 'FontSize', 18), saveas(gcf, 'ScaledInterMFWSIntensity.png');
% if(~isnan(ScaledIntensityFeatures.ScaledIntraMFWSIntensity)), figure, boxplot(ScaledIntensityFeatures.ScaledIntraMFWSIntensity, Muscle_Quality, 'whisker', 10), ...
%         title('IntraMF WS Intensity Fraction', 'FontSize', 18), saveas(gcf, 'ScaledIntraMFWSIntensity.png'); end


% Append scaled intensity features to our data matrix.
XLS_Data = [XLS_Data, ...
    ScaledIntensityFeatures.ScaledMuscleNSIntensity, ...
    ScaledIntensityFeatures.ScaledInterMFNSIntensity, ...
    ScaledIntensityFeatures.ScaledIntraMFNSIntensity, ...
    ScaledIntensityFeatures.ScaledMuscleFSIntensity, ...
    ScaledIntensityFeatures.ScaledInterMFFSIntensity, ...
    ScaledIntensityFeatures.ScaledIntraMFFSIntensity, ...
    ScaledIntensityFeatures.ScaledMuscleWSIntensity, ...
    ScaledIntensityFeatures.ScaledInterMFWSIntensity, ...
    ScaledIntensityFeatures.ScaledIntraMFWSIntensity];

Variable_Names = [Variable_Names; ...
    'Scaled Muscle NS Intensity'; 'Scaled InterMF NS Intensity'; 'Scaled IntraMF NS Intensity';
    'Scaled Muscle FS Intensity'; 'Scaled InterMF FS Intensity'; 'Scaled IntraMF FS Intensity';
    'Scaled Muscle WS Intensity'; 'Scaled InterMF WS Intensity'; 'Scaled IntraMF WS Intensity'];

NewXLS_Data = XLS_Data;
NewVariable_Names = Variable_Names;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NewXLS_Data, NewVariable_Names] = ComputeCovarianceMatrixFeatures(XLS_Data, ...
    Variable_Names)

N = size(XLS_Data, 1);
MuscleCovarianceMatrixDeterminant = zeros(N, 1);
MuscleCovarianceMatrixEigenVector11 = zeros(N, 1);
MuscleCovarianceMatrixEigenVector12 = zeros(N, 1);
MuscleCovarianceMatrixEigenValue1 = zeros(N, 1);
MuscleCovarianceMatrixEigenValue2 = zeros(N, 1);

InterMFCovarianceMatrixDeterminant = zeros(N, 1);
InterMFCovarianceMatrixEigenVector11 = zeros(N, 1);
InterMFCovarianceMatrixEigenVector12 = zeros(N, 1);
InterMFCovarianceMatrixEigenValue1 = zeros(N, 1);
InterMFCovarianceMatrixEigenValue2 = zeros(N, 1);

IntraMFCovarianceMatrixDeterminant = zeros(N, 1);
IntraMFCovarianceMatrixEigenVector11 = zeros(N, 1);
IntraMFCovarianceMatrixEigenVector12 = zeros(N, 1);
IntraMFCovarianceMatrixEigenValue1 = zeros(N, 1);
IntraMFCovarianceMatrixEigenValue2 = zeros(N, 1);

% Compute determinants.
for i=1:N

    Muscle_Covariance_Matrix = zeros(2,2);
    Muscle_Covariance_Matrix(1,1) = XLS_Data(i, 23);
    Muscle_Covariance_Matrix(1,2) = XLS_Data(i, 24);
    Muscle_Covariance_Matrix(2,2) = XLS_Data(i, 25);
    Muscle_Covariance_Matrix(2,1) = Muscle_Covariance_Matrix(1,2);

    Muscle_Covariance_Matrix_Measures = CovarianceMatrixAlgebra(Muscle_Covariance_Matrix);
    MuscleCovarianceMatrixEigenVector11(i) = Muscle_Covariance_Matrix_Measures.eigenVector11;
    MuscleCovarianceMatrixEigenVector12(i) = Muscle_Covariance_Matrix_Measures.eigenVector12;
    MuscleCovarianceMatrixEigenValue1(i) = Muscle_Covariance_Matrix_Measures.eigenValue1;
    MuscleCovarianceMatrixEigenValue2(i) = Muscle_Covariance_Matrix_Measures.eigenValue2;
    MuscleCovarianceMatrixDeterminant(i) = Muscle_Covariance_Matrix_Measures.determinant;


    InterMF_Covariance_Matrix = zeros(2,2);
    InterMF_Covariance_Matrix(1,1) = XLS_Data(i, 28);
    InterMF_Covariance_Matrix(1,2) = XLS_Data(i, 29);
    InterMF_Covariance_Matrix(2,2) = XLS_Data(i, 30);
    InterMF_Covariance_Matrix(2,1) = InterMF_Covariance_Matrix(1,2);

    InterMF_Covariance_Matrix_Measures = CovarianceMatrixAlgebra(InterMF_Covariance_Matrix);
    InterMFCovarianceMatrixEigenVector11(i) = InterMF_Covariance_Matrix_Measures.eigenVector11;
    InterMFCovarianceMatrixEigenVector12(i) = InterMF_Covariance_Matrix_Measures.eigenVector12;
    InterMFCovarianceMatrixEigenValue1(i) = InterMF_Covariance_Matrix_Measures.eigenValue1;
    InterMFCovarianceMatrixEigenValue2(i) = InterMF_Covariance_Matrix_Measures.eigenValue2;
    InterMFCovarianceMatrixDeterminant(i) = InterMF_Covariance_Matrix_Measures.determinant;


    IntraMF_Covariance_Matrix = zeros(2,2);
    IntraMF_Covariance_Matrix(1,1) = XLS_Data(i, 33);
    IntraMF_Covariance_Matrix(1,2) = XLS_Data(i, 34);
    IntraMF_Covariance_Matrix(2,2) = XLS_Data(i, 35);
    IntraMF_Covariance_Matrix(2,1) = IntraMF_Covariance_Matrix(1,2);

    IntraMF_Covariance_Matrix_Measures = CovarianceMatrixAlgebra(IntraMF_Covariance_Matrix);
    IntraMFCovarianceMatrixEigenVector11(i) = IntraMF_Covariance_Matrix_Measures.eigenVector11;
    IntraMFCovarianceMatrixEigenVector12(i) = IntraMF_Covariance_Matrix_Measures.eigenVector12;
    IntraMFCovarianceMatrixEigenValue1(i) = IntraMF_Covariance_Matrix_Measures.eigenValue1;
    IntraMFCovarianceMatrixEigenValue2(i) = IntraMF_Covariance_Matrix_Measures.eigenValue2;
    IntraMFCovarianceMatrixDeterminant(i) = IntraMF_Covariance_Matrix_Measures.determinant;

end


CovMatFeatures.MuscleCovarianceMatrixEigenVector11 = MuscleCovarianceMatrixEigenVector11;
CovMatFeatures.MuscleCovarianceMatrixEigenVector12 = MuscleCovarianceMatrixEigenVector12;
CovMatFeatures.MuscleCovarianceMatrixEigenValue1 = MuscleCovarianceMatrixEigenValue1;
CovMatFeatures.MuscleCovarianceMatrixEigenValue2 = MuscleCovarianceMatrixEigenValue2;
CovMatFeatures.MuscleCovarianceMatrixDeterminant = MuscleCovarianceMatrixDeterminant;


CovMatFeatures.InterMFCovarianceMatrixEigenVector11 = InterMFCovarianceMatrixEigenVector11;
CovMatFeatures.InterMFCovarianceMatrixEigenVector12 = InterMFCovarianceMatrixEigenVector12;
CovMatFeatures.InterMFCovarianceMatrixEigenValue1 = InterMFCovarianceMatrixEigenValue1;
CovMatFeatures.InterMFCovarianceMatrixEigenValue2 = InterMFCovarianceMatrixEigenValue2;
CovMatFeatures.InterMFCovarianceMatrixDeterminant = InterMFCovarianceMatrixDeterminant;


CovMatFeatures.IntraMFCovarianceMatrixEigenVector11 = IntraMFCovarianceMatrixEigenVector11;
CovMatFeatures.IntraMFCovarianceMatrixEigenVector12 = IntraMFCovarianceMatrixEigenVector12;
CovMatFeatures.IntraMFCovarianceMatrixEigenValue1 = IntraMFCovarianceMatrixEigenValue1;
CovMatFeatures.IntraMFCovarianceMatrixEigenValue2 = IntraMFCovarianceMatrixEigenValue2;
CovMatFeatures.IntraMFCovarianceMatrixDeterminant = IntraMFCovarianceMatrixDeterminant;

% Box plots
% figure, boxplot(CovMatFeatures.MuscleCovarianceMatrixDeterminant, Muscle_Quality_Class, 'whisker', 10), ...
%     title('Muscle Cov. Matrix Det.', 'FontSize', 18), saveas(gcf, 'MuscleCovarianceMatrixDeterminant.png');
% figure, boxplot(CovMatFeatures.InterMFCovarianceMatrixDeterminant, Muscle_Quality_Class, 'whisker', 10), ...
%     title('InterMF Cov. Matrix Det.', 'FontSize', 18), saveas(gcf, 'InterMFCovarianceMatrixDeterminant.png');
% if(CovMatFeatures.IntraMFCovarianceMatrixDeterminant~=0), figure, boxplot(CovMatFeatures.IntraMFCovarianceMatrixDeterminant, Muscle_Quality_Class, 'whisker', 10), ...
%         title('IntraMF Cov. Matrix Det.', 'FontSize', 18), saveas(gcf, 'IntraMFCovarianceMatrixDeterminant.png'); end
% Scatter plots vs. Age
%     figure, plot(Age, CovMatFeatures.MuscleCovarianceMatrixDeterminant, '.'), ...
%         title('Muscle Cov. Matrix Det. vs. Age', 'FontSize', 18), saveas(gcf, 'MuscleCovarianceMatrixDeterminantvsAge.png');
%     figure, plot(Age, CovMatFeatures.InterMFCovarianceMatrixDeterminant, '.'), ...
%         title('InterMF Cov. Matrix Det. vs. Age', 'FontSize', 18), saveas(gcf, 'InterMFCovarianceMatrixDeterminantvsAge.png');
%     if(~isnan(CovMatFeatures.IntraMFCovarianceMatrixDeterminant)),     figure, plot(Age, CovMatFeatures.IntraMFCovarianceMatrixDeterminant, '.'), ...
%         title('IntraMF Cov. Matrix Det. vs. Age', 'FontSize', 18), saveas(gcf, 'IntraMFCovarianceMatrixDeterminantvsAge.png'); end


% Append cov. matrix features to our data matrix.
XLS_Data = [XLS_Data, ...
    CovMatFeatures.MuscleCovarianceMatrixEigenVector11, ...
    CovMatFeatures.MuscleCovarianceMatrixEigenVector12, ...
    CovMatFeatures.MuscleCovarianceMatrixEigenValue1, ...
    CovMatFeatures.MuscleCovarianceMatrixEigenValue2, ...
    CovMatFeatures.MuscleCovarianceMatrixDeterminant, ...
    CovMatFeatures.InterMFCovarianceMatrixEigenVector11, ...
    CovMatFeatures.InterMFCovarianceMatrixEigenVector12, ...
    CovMatFeatures.InterMFCovarianceMatrixEigenValue1, ...
    CovMatFeatures.InterMFCovarianceMatrixEigenValue2, ...
    CovMatFeatures.InterMFCovarianceMatrixDeterminant, ...
    CovMatFeatures.IntraMFCovarianceMatrixEigenVector11, ...
    CovMatFeatures.IntraMFCovarianceMatrixEigenVector12, ...
    CovMatFeatures.IntraMFCovarianceMatrixEigenValue1, ...
    CovMatFeatures.IntraMFCovarianceMatrixEigenValue2, ...
    CovMatFeatures.IntraMFCovarianceMatrixDeterminant];

Variable_Names = [Variable_Names; ...
    'Muscle Cov. Matrix Eigenvector11'; 'Muscle Cov. Matrix Eigenvector12'; ...
    'Muscle Cov. Matrix Eigenvalue1'; 'Muscle Cov. Matrix Eigenvalue2'; 'Muscle Cov. Matrix Det.'; ...
    'InterMF Cov. Matrix Eigenvector11'; 'InterMF Cov. Matrix Eigenvector12'; ...
    'InterMF Cov. Matrix Eigenvalue1'; 'InterMF Cov. Matrix Eigenvalue2'; ...
    'InterMF Cov. Matrix Det.'; ...
    'IntraMF Cov. Matrix Eigenvector11'; 'IntraMF Cov. Matrix Eigenvector12'; ...
    'IntraMF Cov. Matrix Eigenvalue1'; 'IntraMF Cov. Matrix Eigenvalue2'; ...
    'IntraMF Cov. Matrix Det.'];

NewXLS_Data = XLS_Data;
NewVariable_Names = Variable_Names;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CovarianceMatrixMeasures = CovarianceMatrixAlgebra(CovarianceMatrix)

[EigenVectors, Diagonal] = eig(CovarianceMatrix);

eigenValues = [Diagonal(1,1), Diagonal(2,2)];
[SortedEigenValues, I] = sort(eigenValues);

CovarianceMatrixMeasures.eigenValue1 = SortedEigenValues(2);
CovarianceMatrixMeasures.eigenValue2 = SortedEigenValues(1);

CovarianceMatrixMeasures.eigenVector11 = EigenVectors(1,I(2));
CovarianceMatrixMeasures.eigenVector12 = EigenVectors(2,I(2));

CovarianceMatrixMeasures.determinant = det(CovarianceMatrix);

end