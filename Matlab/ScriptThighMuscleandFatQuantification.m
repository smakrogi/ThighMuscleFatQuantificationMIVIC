% ScriptThighMuscleandFatQuantification

% Read dataset info.
ExperimentFiles = { 'FST13D_List_Nii_Preprocessed_20120713_165626.txt', ...
    'SubjectIDs_List_Nii_Preprocessed_20120713_165626.txt', ...
    'T13D_List_Nii_Preprocessed_20120713_165626.txt', ...
    'WST13D_List_Nii_Preprocessed_20120713_165626.txt'};
PathName = '/home/makrogianniss/harbor3t_home/matlab/MRIThighAnalysis/MuscleQualityTests';
ExperimentInfo = SetDataPathandFilenamesThighNew(ExperimentFiles, ...
    PathName);

% Process all data in the lists.
nParticipants = length(cat(1,ExperimentInfo.SubjectID));
for (subjectNumber=1:nParticipants)
    [stats, InternalSATSurfaceMask, MuscleMask] = ThighQuantificationCL(ExperimentFiles, PathName, ...
        subjectNumber);
    close all
end