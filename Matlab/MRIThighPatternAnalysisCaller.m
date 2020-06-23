function Statistics = MRIThighPatternAnalysisCaller( filename, varargin )
% syntax: Statistics = MRIThighPatternAnalysisCaller( filename, classification_arguments );
% Run statistics and discriminant analysis algorithms.
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

%% Set useful parameters.
if(~isempty(varargin))
    MachineLearningParams = MRIThighSetMachineLearningParameters(varargin{1});
else
    MachineLearningParams = MRIThighSetMachineLearningParameters;
end

%% Retrieve targets and patterns.
XLS_Data = csvread(filename, 1, 1);
% targets = XLS_Data(:, end-1);
MQ = XLS_Data(:, end-3);
targets = uint8(MQ >= MachineLearningParams.threshold_MQ);

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
DataStruct.targets = targets';
DataStruct.subjectids = Array_SubjectIDs;
DataStruct.featurenames = Array_Variable_Names;

% Run this using try/catch to handle errors coming from single class labels in
% the test set.
try
    
    Statistics = StatisticalPatternAnalysis( DataStruct, MachineLearningParams );
    
catch ME

    % Display error and continue.
    Statistics = [];
    fprintf( 'Error description: %s, %s\n', ME.identifier, ME.message );
        
end

end
