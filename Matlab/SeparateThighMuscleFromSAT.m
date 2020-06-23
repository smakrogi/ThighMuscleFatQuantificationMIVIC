function [stats, InternalSATSurfaceMask, TissueLabelVolume2] = ...
    SeparateThighMuscleFromSAT(configFile, subjectNumber, rescaleFactor, nClusters)
% Obsolete function.
% syntax: [stats, InternalSATSurfaceMask, MuscleMask] =
% SeparateThighMuscleFromSAT(configFile, subjectNumber, rescaleFactor, nClusters);
% rescaleFactor: factor for rescaling of axial plane in clustering.
% nClusters: number of clusters.
% 1. Parametric deformable model segmentation.
% 2. Separation of VAT and SAT using clustering.
% 3. Computation of muscle and fat statistics.
% 3T NIA/NIH, by S.K. Makrogiannis. 
tStart = tic;

% Read dataset info.
% SetDataPathandFilenames2;
% SetDataPathandFilenamesLinux;
% SetDataPathandFilenamesThigh;
eval(configFile);

lengthN = length(subjectNumber);
for i=subjectNumber(1):subjectNumber(lengthN)
    
    % Parametric deformable model segmentation of the
    % internal subcutaneous surface.
    infoString = ['\nSubject: ' ,   ExperimentInfo.SubjectID{i}, '\n'];
    fprintf(infoString);
    [InternalSATSurfaceMask, OneLegVolume, BoundingBox, NoBoneMask] = ...
        SegmentThighInternalSATSurfaceVolume...
        (i, ExperimentInfo);
    
    % Morphological operations followed by
    % central clustering to define all fat voxels.
    [TissueLabelVolume, fatLabel] = ProcessThighImages...
        (i, ExperimentInfo, rescaleFactor, nClusters, ...
        [ExperimentInfo.SubjectID{i} '_InternalSATSurfaceMask.hdr'], ...
        BoundingBox, NoBoneMask);
    % disp('Clustering completed, press any key to continue.');
    % pause;
    % close all;
    
    % Calulcate muscle, SAT, and IMFAT tissue volumes.
    [stats, TissueLabelVolume2, infoString] = CalculateThighMuscleandAdiposeVolumes...
        (i, TissueLabelVolume, OneLegVolume, ...
        InternalSATSurfaceMask, NoBoneMask, ExperimentInfo, BoundingBox, fatLabel);
    
    tElapsed=toc(tStart);
    
    infoString = sprintf('\nElapsed time #%d(sec) \n', tElapsed);
    fprintf(infoString);
    
end

end

