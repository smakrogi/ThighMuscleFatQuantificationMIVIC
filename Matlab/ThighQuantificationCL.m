function [stats, TissueLabelVolume2] = ...
    ThighQuantificationCL(ExperimentFiles, PathName, subjectNumber)
% syntax: [stats, TissueLabelVolume2] =
% ThighMuscleandFatQuantification(ExperimentFiles, subjectNumber);
% 1. Parametric deformable model segmentation.
% 2. Separation of muscle and IMAT using clustering.
% 3. Computation of muscle and fat statistics.
% 3T NIA/NIH, by S.K. Makrogiannis. 

% Read dataset info.
ExperimentInfo = SetDataPathandFilenamesThighNew(ExperimentFiles, ...
    PathName);

lengthN = length(subjectNumber);
for i=subjectNumber(1):subjectNumber(lengthN)

    % Try/catch block for exception handling.
    try
        
        % Start counting time.
        tStart = tic;
        
        % Parametric deformable model segmentation of the
        % internal subcutaneous surface.
        infoString = ['\nSubject: ' ,   ExperimentInfo.SubjectID{i}, '\n'];
        fprintf(infoString);
  
        %     [InternalSATSurfaceMask, OneLegVolume, BoundingBox, NoBoneMask] = ...
        %         SegmentThighInternalSATSurfaceVolume...
        %         (i, ExperimentInfo);
        
        % Morphological operations followed by
        % central clustering to define all fat voxels.
        [TissueLabelVolume, tissueLabels, DistributionModels] = ProcessThighImages...
            (i, ExperimentInfo, ...
            [ExperimentInfo.SubjectID{i} '_InternalSATSurfaceMask.hdr'], ...
            [ExperimentInfo.SubjectID{i} '_OneLegVolume.hdr'], ...
            [ExperimentInfo.SubjectID{i} '_BoneLabelVolume.hdr']);
        % disp('Clustering completed, press any key to continue.');
        % pause;
        % close all;
        clear TissueLabelVolume;
        
        % Calulcate muscle, SAT, and IMFAT tissue volumes.
        imagefileprefix = [ExperimentInfo.SubjectID{i}, '_', ...
            ExperimentInfo.DistanceLearningMethod, '_', ...
            ExperimentInfo.ClusteringMethod];
        [stats, TissueLabelVolume2, infoString] = CalculateThighMuscleandAdiposeVolumes...
            (i, ExperimentInfo, tissueLabels, DistributionModels, ...
            [imagefileprefix '_TissueLabelVolume.hdr'], ...
            [ExperimentInfo.SubjectID{i} '_OneLegVolume.hdr'], ...
            [ExperimentInfo.SubjectID{i} '_InternalSATSurfaceMask.hdr'], ...
            [ExperimentInfo.SubjectID{i} '_BoneLabelVolume.hdr']);
        fprintf(infoString);
        
        tElapsed=toc(tStart);
        infoString = sprintf('\nElapsed time %.2f(sec) \n', tElapsed);
        fprintf(infoString);
        
    catch ME
        
        % Display error and continue.
        fprintf( 'Error description: %s\n', ME.message );
        fprintf( 'Continuing with next image in batch mode.\n');
        
    end
    
end

end

