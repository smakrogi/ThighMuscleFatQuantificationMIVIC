function varargout = ThighQuantification(varargin)
% THIGHQUANTIFICATION M-file for ThighQuantification.fig
%      THIGHQUANTIFICATION, by itself, creates a new THIGHQUANTIFICATION or raises the existing
%      singleton*.
%
%      H = THIGHQUANTIFICATION returns the handle to a new THIGHQUANTIFICATION or the handle to
%      the existing singleton*.
%
%      THIGHQUANTIFICATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THIGHQUANTIFICATION.M with the given input arguments.
%
%      THIGHQUANTIFICATION('Property','Value',...) creates a new THIGHQUANTIFICATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThighQuantification_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThighQuantification_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThighQuantification

% Last Modified by GUIDE v2.5 25-Jun-2012 20:10:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThighQuantification_OpeningFcn, ...
                   'gui_OutputFcn',  @ThighQuantification_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ThighQuantification is made visible.
function ThighQuantification_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThighQuantification (see VARARGIN)

% Choose default command line output for ThighQuantification
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ThighQuantification wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ThighQuantification_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_ReadSubjects.
function pushbutton_ReadSubjects_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ReadSubjects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Open a dialog box for file reading.

% Read data from current figure.
FigData = guidata(gcf);

[FileNames,PathName]=uigetfile({'*.txt','*.txt files'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Pick a configuration file','MultiSelect','on');

% Put together path and configuration filename.
% configFile = [PathName, FileName];
% [temp, name, ext] = fileparts(FileName) ;
% clear temp;
% eval(name);
ExperimentInfo = SetDataPathandFilenamesThighNew(FileNames,PathName);
FigData.FileNames = FileNames;
FigData.PathName = PathName;
FigData.ExperimentInfo = ExperimentInfo;
FigData.slice_Number = 1;

% Display subject list..
set(handles.listboxParticipants,'String',ExperimentInfo.SubjectID)

% Copy the updated Figdata structure to the gcf data.
guidata(gcf,FigData);

% --- Executes on button press in pushbutton_RemoveSAT.
function pushbutton_RemoveSAT_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_RemoveSAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data from figure.
FigData = guidata(gcf);

% Clear the final result figure.
cla(handles.axesNoSubcutaneousFat, 'reset');
cla(handles.axesBar, 'reset');
cla(handles.axesBarH, 'reset');
cla(handles.axesSegmentationMask);
cla(handles.axesColorbar);
set(handles.text_Quantification, 'String', ' ');

% Read the selected participant.
formatIdentifier = get(handles.listboxParticipants,'Value');
formatIdentifierString = get(handles.listboxParticipants,'String');
volumeFormat = cell2mat(formatIdentifierString(formatIdentifier));

% Run the first step in the process.
FigData.subjectNumber = formatIdentifier;

% ExperimentFiles.FileNames = FigData.FileNames;
% ExperimentFiles.PathName = FigData.PathName;
%     [stats, InternalSATSurfaceMask, TissueLabelVolume2] = ...
% ThighQuantificationCL(ExperimentFiles, FigData.subjectNumber);

infoString = ['\nSubject: ' ,   FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '\n'];
fprintf(infoString);
[FigData.InternalSATSurfaceMask, FigData.OneLegVolume, ...
    FigData.BoundingBox, FigData.NoBoneMask] = ...
    SegmentThighInternalSATSurfaceVolume...
    (FigData.subjectNumber, FigData.ExperimentInfo);

FigData.number_of_Slices = size(FigData.OneLegVolume, 3);
FigData.slice_Number = round(FigData.number_of_Slices / 2);
showImage(handles.axesInput, ...
    FigData.OneLegVolume, ...
    FigData.slice_Number, 'input');
showImage(handles.axesNoSubcutaneousFat, ...
    FigData.InternalSATSurfaceMask, ...
    FigData.slice_Number, 'label');

SetSliderBar(FigData, handles);
cla(handles.axesSegmentationMask);

% Copy the updated Figdata structure to the gcf data.
guidata(gcf,FigData);


% --- Executes on button press in pushbuttonDetectIMFat.
function pushbuttonDetectIMFat_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDetectIMFat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data from figure.
FigData = guidata(gcf);

% Central clustering to define all fat voxels.
[FigData.TissueLabelVolume, FigData.tissueLabels, FigData.DistributionModels] = ProcessThighImages...
    (FigData.subjectNumber, FigData.ExperimentInfo, ...
    [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber} '_InternalSATSurfaceMask.hdr'], ...
    [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '_OneLegVolume.hdr'], ...
    [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '_BoneLabelVolume.hdr'] );

% showImage(handles.axesNoSubcutaneousFat, FigData.TissueLabelVolume, ...
%     FigData.slice_Number, 'label');
% PlotImageSamplesAndSave(handles.axesNoSubcutaneousFat, FigData);
LoadAndPlotImageSamples(handles, FigData);
cla(handles.axesSegmentationMask);

% Copy the updated Figdata structure to the gcf data.
guidata(gcf,FigData);


% --- Executes on button press in pushbuttonQuantify.
function pushbuttonQuantify_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonQuantify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data from figure.
FigData = guidata(gcf);

% % Calulcate muscle, SAT, and IMFAT tissue volumes.
imagefileprefix = [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '_', ...
    FigData.ExperimentInfo.DistanceLearningMethod, ...
    '_', FigData.ExperimentInfo.ClusteringMethod];
[temp, FigData.TissueLabelVolume2, FigData.infoString] = CalculateThighMuscleandAdiposeVolumes...
    (FigData.subjectNumber, FigData.ExperimentInfo, FigData.tissueLabels, FigData.DistributionModels, ...
    [imagefileprefix '_TissueLabelVolume.hdr'], ...
    [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '_OneLegVolume.hdr'], ...
    [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber} '_InternalSATSurfaceMask.hdr'], ...
    [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '_BoneLabelVolume.hdr'] );
clear temp;

% Display metadata.
set(handles.text_Quantification,'String', FigData.infoString)

showImage(handles.axesSegmentationMask, FigData.TissueLabelVolume2, ...
     FigData.slice_Number, 'label');
displayColorbar(handles.axesColorbar);

% Copy the updated Figdata structure to the gcf data.
guidata(gcf,FigData);



% --- Executes on button press in pushbuttonRunBatchProcess.
function pushbuttonRunBatchProcess_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRunBatchProcess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 

% ExperimentFiles.FileNames = FigData.FileNames;
% ExperimentFiles.PathName = FigData.PathName;
% nParticipants = length(cat(1,handles.ExperimentInfo.SubjectID));
% for (subjectNumber=1:nParticipants)
%     [stats, InternalSATSurfaceMask, TissueLabelVolume2] = ...
% ThighQuantificationCL(ExperimentFiles, subjectNumber);
% %     close all
% end

FigData = guidata(gcf);

% Clear the final result figure.
cla(handles.axesNoSubcutaneousFat);
cla(handles.axesBar);
cla(handles.axesBarH);
cla(handles.axesSegmentationMask);
set(handles.text_Quantification, 'String', ' ');

% Start from the selected dataset.
formatIdentifier = get(handles.listboxParticipants,'Value');
firstSubject = formatIdentifier;
nParticipants = length(cat(1,handles.ExperimentInfo.SubjectID));
subjectIT = firstSubject;

while (get(handles.pushbuttonRunBatchProcess,'Value')==1.0 && ...
        subjectIT <= nParticipants)
    % for (subjectIT=firstSubject:nParticipants)
    
    % Try/catch block for exception handling.
    try
        % Start time counting.
        tStart = tic;
        
        FigData.subjectNumber = subjectIT;
        
        % % Clear the final result figure.
        % cla(handles.axesNoSubcutaneousFat);
        % cla(handles.axesSegmentationMask);
        % set(handles.text_Quantification, 'String', ' ');
        
        % 1st step.
        infoString = ['\nSubject: ' ,   FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '\n'];
        fprintf(infoString);
        
        [FigData.InternalSATSurfaceMask, FigData.OneLegVolume, ...
            FigData.BoundingBox, FigData.NoBoneMask] = ...
            SegmentThighInternalSATSurfaceVolume...
            (FigData.subjectNumber, FigData.ExperimentInfo);
        
        FigData.number_of_Slices = size(FigData.OneLegVolume, 3);
        FigData.slice_Number = round(FigData.number_of_Slices / 2);
        showImage(handles.axesInput, FigData.OneLegVolume, ...
            FigData.slice_Number, 'input');
        SetSliderBar(FigData, handles);
        
        showImage(handles.axesNoSubcutaneousFat, ...
            FigData.InternalSATSurfaceMask, ...
            FigData.slice_Number, 'label');
        % colorbar;
        
        % Copy the updated Figdata structure to the gcf data.
        guidata(gcf,FigData);
        
        %     % 2nd step.
        pushbuttonDetectIMFat_Callback(hObject, eventdata, handles);

        %     % 3rd step.
        pushbuttonQuantify_Callback(hObject, eventdata, handles);
        
        % Display elapsed time.
        tElapsed=toc(tStart);
        infoString = sprintf('\nElapsed time: %.2f(sec) \n', tElapsed);
        fprintf(infoString);
        
%         % Close all figures left from 'command line' mode.
%         close all;
        
    %Basic exception handling.
    catch ME
        
        % Display error and continue.
        fprintf( 'Error description: %s\n', ME.message );
        fprintf( 'Continuing with next image in batch mode.\n');
        
    end
    
    subjectIT = subjectIT + 1;
    
end

% Toggle off batch processing.
set(handles.pushbuttonRunBatchProcess,'Value', 0.0)

% Copy the updated Figdata structure to the gcf data.
guidata(gcf,FigData);



% --- Executes on selection change in listboxParticipants.
function listboxParticipants_Callback(hObject, eventdata, handles)
% hObject    handle to listboxParticipants (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxParticipants contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxParticipants
% Load gui data.
FigData = guidata(gcf);

% Clear the final result figure.
cla(handles.axesInput);
cla(handles.axesNoSubcutaneousFat, 'reset');
cla(handles.axesBar, 'reset');
cla(handles.axesBarH, 'reset');
cla(handles.axesSegmentationMask);
cla(handles.axesColorbar);
set(handles.text_Quantification, 'String', ' ');

% Read the selected participant.
formatIdentifier = get(handles.listboxParticipants,'Value');

% Run the first step in the process.
FigData.subjectNumber = formatIdentifier;


% Input volume.
onelegvolumeFilename = [ FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, ...
    '_OneLegVolume.hdr'];
fid = fopen( onelegvolumeFilename );
if fid~=-1
    FigData.OneLegVolume = load_untouch_nii( onelegvolumeFilename );
    FigData.OneLegVolume = FigData.OneLegVolume.img;
    % FigData.OneLegVolume = TransposeVolume(FigData.OneLegVolume);
    % Update slider.
    FigData.number_of_Slices = size(FigData.OneLegVolume, 3);
    FigData.slice_Number = round(FigData.number_of_Slices / 2);
    SetSliderBar(FigData, handles);
    % Update slider.
    FigData.number_of_Slices = size(FigData.OneLegVolume, 3);
    FigData.slice_Number = round(FigData.number_of_Slices / 2);
    SetSliderBar(FigData, handles);
    showImage(handles.axesInput, ...
        FigData.OneLegVolume, ...
        FigData.slice_Number, 'input');
    fclose( fid );
end


% IMFAT Cluster map.
imagefileprefix = [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '_', ...
    FigData.ExperimentInfo.DistanceLearningMethod, '_', ...
    FigData.ExperimentInfo.ClusteringMethod];

fid = fopen( [imagefileprefix,  '_TissueLabelVolume.hdr']);
if fid~=-1
    FigData.TissueLabelVolume =  load_untouch_nii([imagefileprefix, '_TissueLabelVolume.hdr']);
    FigData.TissueLabelVolume = FigData.TissueLabelVolume.img;
%     FigData.TissueLabelVolume = TransposeVolume(FigData.TissueLabelVolume);
    % Update slider.
    FigData.number_of_Slices = size(FigData.TissueLabelVolume, 3);
    FigData.slice_Number = round(FigData.number_of_Slices / 2);
    SetSliderBar(FigData, handles);
    % Update slider.
    FigData.number_of_Slices = size(FigData.TissueLabelVolume, 3);
    FigData.slice_Number = round(FigData.number_of_Slices / 2);
    SetSliderBar(FigData, handles);
    showImage(handles.axesNoSubcutaneousFat, ...
        FigData.TissueLabelVolume, ...
        FigData.slice_Number, 'label');
    fclose( fid );
end


% Scatter plot of WS and FS intensities.
FigData.DistributionModels.Samples = ...
    LoadAndPlotImageSamples(handles, ...
    FigData);


% Final tissue labels.
fid = fopen( [imagefileprefix,  '_TissueLabelVolume2.hdr']);
if fid~=-1
    FigData.TissueLabelVolume2 = load_untouch_nii([imagefileprefix '_TissueLabelVolume2.hdr']);
    FigData.TissueLabelVolume2 = FigData.TissueLabelVolume2.img;
%     FigData.TissueLabelVolume2 = TransposeVolume(FigData.TissueLabelVolume2);
    % Update slider.
    FigData.number_of_Slices = size(FigData.TissueLabelVolume2, 3);
    FigData.slice_Number = round(FigData.number_of_Slices / 2);
    SetSliderBar(FigData, handles);
    % Update slider.
    FigData.number_of_Slices = size(FigData.TissueLabelVolume2, 3);
    FigData.slice_Number = round(FigData.number_of_Slices / 2);
    SetSliderBar(FigData, handles);
    showImage(handles.axesSegmentationMask, ...
        FigData.TissueLabelVolume2, ...
        FigData.slice_Number, 'label');
    fclose( fid );
    displayColorbar(handles.axesColorbar);
end

% Quantification results.
fid=fopen([imagefileprefix '_ThighQuantification.csv']);

if fid~=-1
%     infoString = fscanf(fid, '%s');
    infoString = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'Delimiter', ',' , 'CollectOutput',1);
    infoString2 = [ infoString{1}(1,:); infoString{1}(2,:) ];
    % Display metadata.
    set(handles.text_Quantification,'String', infoString2)
    fclose(fid);
end

% Copy the updated Figdata structure to the gcf data.
guidata(gcf,FigData);



% --- Executes during object creation, after setting all properties.
function listboxParticipants_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxParticipants (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
ThighQuantificationParameters;
 

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderSliceNumber_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSliceNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

switch get(handles.sliderSliceNumber, 'Enable')
    case 'on'
        sliceNumber = round(get(handles.sliderSliceNumber, 'Value'));
    case {'inactive', 'off'}
        sliceNumber = 1;
end

if(strcmp(get(handles.axesInput, 'Visible'), 'on'))
    showImage(handles.axesInput, handles.OneLegVolume, ...
        sliceNumber, 'input');
end
% if(strcmp(get(handles.axesNoSubcutaneousFat, 'Visible'), 'on'))
% showImage(handles.axesNoSubcutaneousFat, handles.InternalSATSurfaceMask, ...
%     sliceNumber);
% end
if(strcmp(get(handles.axesNoSubcutaneousFat, 'Visible'), 'on'))
showImage(handles.axesNoSubcutaneousFat, handles.TissueLabelVolume, ...
    sliceNumber, 'label');
end
if(strcmp(get(handles.axesSegmentationMask, 'Visible'), 'on'))
showImage(handles.axesSegmentationMask, handles.TissueLabelVolume2, ...
    sliceNumber, 'label');
end


% --- Executes during object creation, after setting all properties.
function sliderSliceNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSliceNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% Displays image and intensity values.
function showImage(handleFigure, Volume, sliceNumber, colormap)

switch lower(colormap)
    case 'input'
        RGB = Volume(:, :, sliceNumber);
    case 'likelihood'
        RGB = ind2rgb(Volume(:, :, sliceNumber), cool);
    case 'label'
        RGB = ind2rgb(Volume(:, :, sliceNumber), cool(6));
end
axes(handleFigure), imagesc(RGB), axis image;
% impixelinfo(handleFigure);


% Scatter plot of pixel intensities.
function PlotImageSamples(handles, FigData)

legendLabels = cell(1, FigData.ExperimentInfo.nClusters);
purple = [119/255 73/255 152/255];
plot_labels = {'k', 'r', 'g', purple, 'm', 'y'};
plotMarkerSize = 2;
plotMarkerShape = '+';
fs_vector = [];
ws_vector = [];

for count=1:FigData.ExperimentInfo.nClusters
    
    % WS and FS scatterplot
    axes(handles.axesNoSubcutaneousFat), axis square
    plot(FigData.DistributionModels.Samples{count}.FS, ...
        FigData.DistributionModels.Samples{count}.WS, ...
        plotMarkerShape, ...
        'MarkerEdgeColor', plot_labels{count}, 'MarkerFaceColor', plot_labels{count}, ...
        'MarkerSize', plotMarkerSize); h1 = gca;
    hold on, xlabel('fat suppressed'); ylabel('water suppressed'); grid on, axis square;

    if count==FigData.tissueLabels.muscle legendLabels{count} = 'muscle';
    elseif count==FigData.tissueLabels.fat legendLabels{count} = 'interMF';
    else legendLabels {count} = 'intraMF';
    end

    fs_vector = [fs_vector; FigData.DistributionModels.Samples{count}.FS];
    ws_vector = [ws_vector; FigData.DistributionModels.Samples{count}.WS];

end

% Build and plot marginal histograms.
nbins = 32;
max_fs = max(fs_vector); min_fs = min(fs_vector);
fs_bin_width = (max_fs-min_fs) / nbins;
[nx,cx]= hist(fs_vector, ...
    (min_fs:fs_bin_width:max_fs));

max_ws = max(ws_vector); min_ws = min(ws_vector);
ws_bin_width = (max_ws-min_ws) / nbins;
[ny,cy]= hist(ws_vector, (min_ws:ws_bin_width:max_ws));

axes(handles.axesNoSubcutaneousFat) 
legend(legendLabels, 'Location', 'NorthEast')
axis([min_fs max_fs min_ws max_ws ])
hold off;

for count=1:FigData.ExperimentInfo.nClusters
    [nx,cx]= hist(FigData.DistributionModels.Samples{count}.FS, ...
        (min_fs:fs_bin_width:max_fs));

    axes(handles.axesBar);
    bar(cx,-nx, 0.8, plot_labels{count}); axis('off'); hold on;

    [ny,cy]= hist(FigData.DistributionModels.Samples{count}.WS, ...
        (min_ws:ws_bin_width:max_ws));

    axes(handles.axesBarH);
    barh(cy,-ny, 0.8, plot_labels{count}); axis('off'); hold on;
end

axes(handles.axesBar); 
xlim([min_fs max_fs])
hold off;

axes(handles.axesBarH); 
ylim([min_ws max_ws])
hold off;
    
        
% Load supplemtentary data.
function SupplementaryData = LoadAndPlotImageSamples(handles, FigData)

imagefileprefix = [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '_', ...
    FigData.ExperimentInfo.DistanceLearningMethod, '_', ...
    FigData.ExperimentInfo.ClusteringMethod];
SupplementaryData = [];

fid = fopen([imagefileprefix,  FigData.ExperimentInfo.SupplementaryDataFilenameSuffix]);
if fid~=-1
    fclose( fid );
    load([imagefileprefix,  FigData.ExperimentInfo.SupplementaryDataFilenameSuffix], ...
        'SupplementaryData');
    FigData.DistributionModels.Samples = SupplementaryData.DataMatrix;
    FigData.tissueLabels = SupplementaryData.tissueLabels;
    PlotImageSamples(handles, FigData)
end


% function PlotImageSamplesAndSave(handleFigure, FigData)
% 
% PlotImageSamples(handleFigure, FigData);
% 
% imagefileprefix = [FigData.ExperimentInfo.SubjectID{FigData.subjectNumber}, '_', ...
%     FigData.ExperimentInfo.DistanceLearningMethod, '_', ...
%     FigData.ExperimentInfo.ClusteringMethod];
% 
% fid = fopen([imagefileprefix,  FigData.ExperimentInfo.SupplementaryDataFilenameSuffix]);
% if fid~=-1
%     fclose( fid );
%     load([imagefileprefix,  FigData.ExperimentInfo.SupplementaryDataFilenameSuffix], ...
%         'SupplementaryData');
% end
% 
% SupplementaryData.DataMatrix = FigData.DistributionModels.Samples;
% SupplementaryData.tissueLabels = FigData.tissueLabels;
% 
% save([imagefileprefix,  FigData.ExperimentInfo.SupplementaryDataFilenameSuffix], ...
%     'SupplementaryData');


% --- Executes on button press in pushbuttonAllsteps.
function pushbuttonAllsteps_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAllsteps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tStart = tic;

% Run all 3 steps sequentially.
pushbutton_RemoveSAT_Callback(hObject, eventdata, handles);
pushbuttonDetectIMFat_Callback(hObject, eventdata, handles);
pushbuttonQuantify_Callback(hObject, eventdata, handles);

tElapsed=toc(tStart);

infoString = sprintf('\nElapsed time %.2f(sec) \n', tElapsed);
fprintf(infoString);

    
% % Colorbar of segmentation labels.
% function displayColorbar(handleFigure)
% % Display the colorbar.
% hold(gca,'off');
% A = zeros(100,10); A(1:50, :) = 4; A(51:100, :) = 2;
% a=linspace(3,0,100);
% A = repmat(a', [1,10]);
% axes(handleFigure), imagesc(round(A)), axis image, axis off;
% text(12, 15, 'Muscle'), text(12, 40, 'IMFAT'), text(12, 60, 'SAT');, text(12, 85, 'Air'); 
% hold(gca,'off');


function displayColorbar(handleFigure)
% Display the colorbar.

A = zeros(100,10); A(1:50, :) = 4; A(51:100, :) = 2;
a = linspace(5.4,0.5,100);
A = round(repmat(a', [1,10]));
RGB = ind2rgb(A, cool(6));
axes(handleFigure), imagesc(RGB), axis image, axis off;
text(12, 10, 'Muscle'), text(12, 30, 'IntraMF'), text(12, 50, 'InterMF'), text(12, 70, 'SAT'); text(12, 90, 'Air'); 

% Set the slider bar
function SetSliderBar(FigData, handles)
% Set the slider bar values.
if FigData.number_of_Slices > 1
    set(handles.sliderSliceNumber, 'Enable', 'on');
    set(handles.sliderSliceNumber, 'max', FigData.number_of_Slices);
    set(handles.sliderSliceNumber, 'min', 1);
    set(handles.sliderSliceNumber, 'Value', FigData.slice_Number);
    slider_step1 = 1 / (FigData.number_of_Slices-1);
%     slider_step1 = 1;
    set(handles.sliderSliceNumber, 'SliderStep', [slider_step1 slider_step1]);
else
    set(handles.sliderSliceNumber, 'Value', slice_Number);
    set(handles.sliderSliceNumber, 'HandleVisibility', 'off');
    set(handles.sliderSliceNumber, 'Enable', 'off');
end


% --- Executes during object creation, after setting all properties.
function axesLogo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axesLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axesLogo

[X,Y] = meshgrid(-12:0.2:12);
R = sqrt(X.^2 + Y.^2) + eps;
Z = sin(R)./R;
% mesh(X,Y,Z), axis off,
surf(X,Y,Z,'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
text(15,10,'BLSA'), axis off, axis tight;
