% MATLAB
%
% Files
%   AnalyzeThighDataInFeatureSpace                     - syntax: [tissueIdx, OutDataMatrix] = ...
%   CalculateThighMuscleandAdiposeVolumes              - syntax: [stats, TissueLabelVolume2, infoString] = CalculateThighMuscleandAdiposeVolumes...
%   createfigure_fat                                   - CREATEFIGURE(X1,Y1,X2,Y2)
%   createfigure_muscle                                - CREATEFIGURE(X1,Y1,X2,Y2)
%   FindMaximumACCParameters                           - Find maximizing parameters of ACC.
%   GenerateBoneMasks                                  - Apply seeded region growing to produce the cortical bone region.
%   MRIThighLibSVMGridSearchTuning                     - syntax: [libsvmACC, libsvmgamma_n, libsvmC_n] = MRIThighLibSVMGridSearchTuning(filename, class_ids);
%   MRIThighPatternAnalysisCaller                      - syntax: Statistics = MRIThighPatternAnalysisCaller( filename, classification_arguments );
%   MRIThighPatternAnalysisCallerROC                   - syntax: ROC_Results = MRIThighPatternAnalysisCallerROC( filename, classification_arguments );
%   MRIThighQuantificationGUI                          - THIGHQUANTIFICATION M-file for ThighQuantification.fig
%   MRIThighSetMachineLearningParameters               - Set machine learning and statistical method parameters.
%   MRIThighStatisticalTests                           - syntax: DataStruct = MRIThighStatisticalTests(filename, selected_gender, n_classes);
%   ProcessThighImages                                 - syntax: [TissueLabelVolumeOutput, fatLabel] = ProcessThighImages...
%   RemoveBoneFromThighVolume                          - Remove the cortical bone.
%   RemoveBoneFromThighVolume2                         - Remove bone marrow.
%   RemoveBoneFromThighVolumeBeforeSATExtraction       - 
%   ScriptThighMuscleandFatQuantification              - ScriptThighMuscleandFatQuantification
%   SegmentThighInternalSATSurfaceVolume               - SegmentInternalSubcutaneousSurfaceVolume
%   SelectOneLeg                                       - Apply threshold to binarize.
%   SeparateThighMuscleFromSAT                         - Obsolete function.
%   SetDataPathandFilenamesThigh                       - SetDataPathandFilenames
%   SetDataPathandFilenamesThighNew                    - SetDataPathandFilenamesThighNew
%   SetDataPathandFilenamesThighPreprocessed           - SetDataPathandFilenames
%   SetDataPathandFilenamesThighPreprocessedRegistered - SetDataPathandFilenamesThighPreprocessedRegistered
%   SetPrinceSnakeOptionsThigh                         - Set snake algorithm options.
%   SetPrinceSnakeOptionsThighParzenGVF                - Set snake algorithm options.
%   subclust                                           - Locates data cluster centers using subtractive clustering. 
%   ThighQuantification                                - M-file for ThighQuantification.fig
%   ThighQuantificationCL                              - syntax: [stats, TissueLabelVolume2] =
%   ThighQuantificationParameters                      - M-file for ThighQuantificationParameters.fig
%   TopHatTransformMultiSlice                          - 
%   TransposeVolume                                    - Used after reading and writing routines.
%   WarpToNonSuppressedVolume                          - syntax: WarpedVolume = WarpToNonSuppressedVolume(inputvolumeFilename,
