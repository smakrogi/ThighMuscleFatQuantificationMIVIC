function [stats, InternalSATSurfaceMask, MuscleMask] = ThighMuscleandFatQuantification(configFile, subjectNumber, rescaleFactor, nClusters)
% syntax: [stats, InternalSATSurfaceMask, MuscleMask] =
% ThighMuscleandFatQuantification(configFile, subjectNumber, rescaleFactor, nClusters);
% rescaleFactor: factor for rescaling of axial plane in clustering.
% nClusters: number of clusters.
% 1. Parametric deformable model segmentation.
% 2. Separation of VAT and SAT using clustering.
% 3. Computation of muscle and fat statistics.

tStart = tic;

% Read dataset info.
% SetDataPathandFilenames2;
% SetDataPathandFilenamesLinux;
% SetDataPathandFilenamesThigh;
eval(configFile);

% Parametric deformable model segmentation of the
% internal subcutaneous surface.
[InternalSATSurfaceMask, OneLegVolume] = ...
    SegmentThighInternalSATSurfaceVolume...
    (subjectNumber, ExperimentInfo);

% Morphological operations followed by
% central clustering to define all fat voxels.
[TissueLabelVolume, fatLabel] = ProcessThighImages...
    (subjectNumber, ExperimentInfo, rescaleFactor, nClusters, ...
    [ExperimentInfo.SubjectID{subjectNumber} '_InternalSATSurfaceMask.hdr'], ...
    [ExperimentInfo.SubjectID{subjectNumber} '_OneLegVolume.hdr']);
% disp('Clustering completed, press any key to continue.');
% pause;
% close all;

% Calulcate muscle, SAT, and IMFAT tissue volumes.
[stats, MuscleMask] = CalculateThighMuscleandAdiposeVolumes...
    (subjectNumber, TissueLabelVolume, OneLegVolume, ...
    InternalSATSurfaceMask, ExperimentInfo, fatLabel);

tElapsed=toc(tStart);

infoString = sprintf('\nElapsed time #%d(sec) \n', tElapsed);
fprintf(infoString);

end




