function WarpedVolume = WarpToNonSuppressedVolume(inputvolumeFilename, targetvolumeFilename)
% syntax: WarpedVolume = WarpToNonSuppressedVolume(inputvolumeFilename,
% targetvolumeFilename);

% ~/Software/bin/BRAINSStandAlone_x86_64/bin/BRAINSFit --movingVolume 20100126_100802BLSA788101_FatSuppressedNoSAT.hdr  --fixedVolume 20100126_100802BLSA788101_NonSuppressedNoSAT.hdr --outputVolume testFS.hdr --transformType BSpline --costFunctionConvergenceFactor 1e10
% Program parameters.
warpingProgram = '~/Software/bin/BRAINSStandAlone_x86_64/bin/BRAINSFit ';
parameterString = '  --transformType Rigid,BSpline --maxBSplineDisplacement 5 ';
parameterString = [parameterString ' --costFunctionConvergenceFactor 1e10 --numberOfHistogramBins 50'];
parameterString = [parameterString ' --initializeTransformMode useMomentsAlign'];

% Create filename of warped image.
[pathstr, name, ext] = fileparts(inputvolumeFilename);
warpedvolumeFilename = [name, '_WarpedToNS.hdr'];
% Creat command string.
commandLine = [warpingProgram, ' --movingVolume ', inputvolumeFilename, ...
    ' --fixedVolume ' targetvolumeFilename, ' --outputVolume ', warpedvolumeFilename, parameterString];

% New path to gcc libraries.
oldPath = getenv('LD_LIBRARY_PATH');
newPath = ['/lib64:/usr/lib64:', oldPath];
setenv('LD_LIBRARY_PATH', newPath);

% Run command.
infoString = 'Elastic registration to non suppressed volume... \n';
infoString = [infoString, commandLine, '\n'];
fprintf(infoString);
system(commandLine);

% Revert to old path.
setenv('LD_LIBRARY_PATH', oldPath);

% Load warped volume to memory and cast to uint16.
WarpedVolume = load_untouch_nii(warpedvolumeFilename);
WarpedVolume = WarpedVolume.img;
% WarpedVolume= uint16( WarpedVolume );

end
