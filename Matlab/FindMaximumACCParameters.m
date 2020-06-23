function [maxMeanACC, argMaxMeanACC, MeanACC, ACC, ImportedData] = FindMaximumACCParameters(dirName)
% Find maximizing parameters of ACC.
% syntax: [maxACC, argMaxACC, MeanACC, ACC, ImportedData] = FindMaximumACCParameters(dirName);

% Get all csv filenames with experimental results in directory.
D = dir([dirName, '*.csv']);

% Import data from each csv file and retrieve ACC values.
for i=1:length(D); 
    system(['grep ''Params'' ', D(i).name,'> text.txt']);  
    ImportedData=importdata('text.txt', ','); 
    ACC{i} = str2num(cat(1,ImportedData.textdata{:,8})); 
end

% Calculate mean ACC for each parameter setting.
MeanACC = zeros(size(ACC{1}));

for i=1:length(ACC), 
    MeanACC = MeanACC + ACC{i}; 
end

MeanACC = MeanACC / length(ACC);

% Find arg max of ACC and display ACC as well.
[maxMeanACC, argMaxMeanACC] = max(MeanACC);

% Display plot of mean ACC vs the parameter setting.
figure, plot(MeanACC);

% List the maximizing parameters.
ImportedData.textdata{argMaxMeanACC, 1:2}

end