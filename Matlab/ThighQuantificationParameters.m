function varargout = ThighQuantificationParameters(varargin)
% THIGHQUANTIFICATIONPARAMETERS M-file for ThighQuantificationParameters.fig
%      THIGHQUANTIFICATIONPARAMETERS, by itself, creates a new THIGHQUANTIFICATIONPARAMETERS or raises the existing
%      singleton*.
%
%      H = THIGHQUANTIFICATIONPARAMETERS returns the handle to a new THIGHQUANTIFICATIONPARAMETERS or the handle to
%      the existing singleton*.
%
%      THIGHQUANTIFICATIONPARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THIGHQUANTIFICATIONPARAMETERS.M with the given input arguments.
%
%      THIGHQUANTIFICATIONPARAMETERS('Property','Value',...) creates a new THIGHQUANTIFICATIONPARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThighQuantificationParameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThighQuantificationParameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThighQuantificationParameters

% Last Modified by GUIDE v2.5 08-Apr-2013 10:09:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThighQuantificationParameters_OpeningFcn, ...
                   'gui_OutputFcn',  @ThighQuantificationParameters_OutputFcn, ...
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


% --- Executes just before ThighQuantificationParameters is made visible.
function ThighQuantificationParameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThighQuantificationParameters (see VARARGIN)

% Choose default command line output for ThighQuantificationParameters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ThighQuantificationParameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ThighQuantificationParameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listboxDR.
function listboxDR_Callback(hObject, eventdata, handles)
% hObject    handle to listboxDR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxDR contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxDR


% --- Executes during object creation, after setting all properties.
function listboxDR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxDR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxClusteringMethod.
function listboxClusteringMethod_Callback(hObject, eventdata, handles)
% hObject    handle to listboxClusteringMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxClusteringMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxClusteringMethod


% --- Executes during object creation, after setting all properties.
function listboxClusteringMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxClusteringMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFirstProcessedSlice_Callback(hObject, eventdata, handles)
% hObject    handle to editFirstProcessedSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFirstProcessedSlice as text
%        str2double(get(hObject,'String')) returns contents of editFirstProcessedSlice as a double


% --- Executes during object creation, after setting all properties.
function editFirstProcessedSlice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFirstProcessedSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editLastProcessedSlice_Callback(hObject, eventdata, handles)
% hObject    handle to editLastProcessedSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLastProcessedSlice as text
%        str2double(get(hObject,'String')) returns contents of editLastProcessedSlice as a double


% --- Executes during object creation, after setting all properties.
function editLastProcessedSlice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLastProcessedSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNClusters_Callback(hObject, eventdata, handles)
% hObject    handle to editNClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNClusters as text
%        str2double(get(hObject,'String')) returns contents of editNClusters as a double


% --- Executes during object creation, after setting all properties.
function editNClusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNScales_Callback(hObject, eventdata, handles)
% hObject    handle to editNScales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNScales as text
%        str2double(get(hObject,'String')) returns contents of editNScales as a double


% --- Executes during object creation, after setting all properties.
function editNScales_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNScales (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editErosionMaskSize_Callback(hObject, eventdata, handles)
% hObject    handle to editErosionMaskSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editErosionMaskSize as text
%        str2double(get(hObject,'String')) returns contents of editErosionMaskSize as a double


% --- Executes during object creation, after setting all properties.
function editErosionMaskSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editErosionMaskSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxSnakeExternalForces.
function listboxSnakeExternalForces_Callback(hObject, eventdata, handles)
% hObject    handle to listboxSnakeExternalForces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxSnakeExternalForces contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxSnakeExternalForces


% --- Executes during object creation, after setting all properties.
function listboxSnakeExternalForces_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxSnakeExternalForces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Pass parameters to the father structure.
ParentFigData = guidata(gcbo);
% ParentFigData.
% ExperimentInfo.snakemaskseSize =
% ExperimentInfo.DistanceLearningMethods =
% ExperimentInfo.OutputDimensionality
% ExperimentInfo.ClusteringMethod
% ExperimentInfo.FirstUseableSlice
% ExperimentInfo.LastUseableSlice
% ExperimentInfo.ProcessedLeg
% nClusters
% ExperimentInfo.nScale

close( gcf );


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Pass parameters to the father structure.
MainAppHandle = findobj( 'name', 'ThighQuantification' );
ParentFigData = guidata( MainAppHandle );
% subjectNumber = ParentFigData.subjectNumber;
subjectNumber = get(ParentFigData.listboxParticipants,'Value');
% ExperimentInfo.FirstUseableSlice
formatIdentifierString = get(handles.editFirstProcessedSlice, 'String');
ParentFigData.ExperimentInfo.FirstUseableSlice = str2num( formatIdentifierString );
% ExperimentInfo.LastUseableSlice
formatIdentifierString = get(handles.editLastProcessedSlice, 'String');
ParentFigData.ExperimentInfo.LastUseableSlice = str2num( formatIdentifierString );
% ExperimentInfo.commandlineMode
formatIdentifierString = get(handles.checkboxCommandLineMode, 'Value');
ParentFigData.ExperimentInfo.commandlineMode = formatIdentifierString;
% nClusters
formatIdentifierString = get(handles.editNClusters, 'String');
ParentFigData.ExperimentInfo.nClusters = str2num( formatIdentifierString );
% ExperimentInfo.nScale
formatIdentifierString = get(handles.editNScales, 'String');
ParentFigData.ExperimentInfo.nScales{subjectNumber} = str2num( formatIdentifierString );
% ExperimentInfo.ProcessedLeg
formatIdentifier = get(handles.listboxSelectLeg, 'Value');
formatIdentifierString = get(handles.listboxSelectLeg,'String');
ParentFigData.ExperimentInfo.ProcessedLeg = cell2mat(formatIdentifierString(formatIdentifier));
% volumeFormat = cell2mat(formatIdentifierString(formatIdentifier));
% ExperimentInfo.DistanceLearningMethods =
formatIdentifier = get(handles.listboxDR,'Value');
ParentFigData.ExperimentInfo.DistanceLearningMethod = ...
    ParentFigData.ExperimentInfo.DistanceLearningMethods{formatIdentifier};
% ExperimentInfo.ClusteringMethod
formatIdentifier = get(handles.listboxClusteringMethod,'Value');
ParentFigData.ExperimentInfo.ClusteringMethod = ...
    ParentFigData.ExperimentInfo.ClusteringMethods{formatIdentifier};
% ExperimentInfo.commandlineMode
formatIdentifierString = get(handles.checkboxWarpWSandFStoNS, 'Value');
ParentFigData.ExperimentInfo.warpWSandFStoNS = formatIdentifierString;
% ExperimentInfo.SATExtractionAlgorithm
formatIdentifier = get(handles.listboxSATExtractionAlgorithms, 'Value');
ParentFigData.ExperimentInfo.SATExtractionAlgorithm = ...
    ParentFigData.ExperimentInfo.SATExtractionAlgorithms{formatIdentifier};
% ExperimentInfo.snakeExternalForce
formatIdentifier = get(handles.listboxSnakeExternalForces,'Value');
ParentFigData.ExperimentInfo.snakeExternalForce = ...
    ParentFigData.ExperimentInfo.snakeExternalForces{formatIdentifier};
% ExperimentInfo.gvfInput
formatIdentifier = get(handles.listboxGVFInput,'Value');
ParentFigData.ExperimentInfo.gvfInput = ...
    ParentFigData.ExperimentInfo.gvfInputs{formatIdentifier};
% ExperimentInfo.snakeGaussianSigma =
formatIdentifierString = get(handles.editSnakeGaussianSigma, 'String');
ParentFigData.ExperimentInfo.snakeGaussianSigma = str2num( formatIdentifierString );
% ExperimentInfo.snakemaskseSize =
formatIdentifierString = get(handles.editErosionMaskSize, 'String');
ParentFigData.ExperimentInfo.snakemaskseSize = str2num( formatIdentifierString );
% ExperimentInfo.editBoneAreaThreshold =
formatIdentifierString = get(handles.editBoneAreaThreshold, 'String');
ParentFigData.ExperimentInfo.boneAreaThreshold = str2num( formatIdentifierString );
% ExperimentInfo.OutputDimensionality
formatIdentifierString = get(handles.editNDimensions, 'String');
ParentFigData.ExperimentInfo.OutputDimensionality = ...
    str2num( formatIdentifierString );
% ExperimentInfo.nComponents = 
formatIdentifierString = get(handles.editNGMMComponents, 'String');
ParentFigData.ExperimentInfo.GMMnDistributions = ...
    str2num( formatIdentifierString );
% ExperimentInfo.rescaleFactor
formatIdentifierString = get(handles.editRescaleFactor, 'String');
ParentFigData.ExperimentInfo.rescaleFactor = ...
    str2num( formatIdentifierString );
% ExperimentInfo.LegSelectionThreshold
formatIdentifierString = get(handles.editLegSelectionThreshold, 'String');
ParentFigData.ExperimentInfo.LegSelectionThreshold{subjectNumber} = ...
    str2num( formatIdentifierString );
% Copy the updated Figdata structure to the gcf data.
guidata(MainAppHandle, ParentFigData);




% --- Executes on selection change in listboxSelectLeg.
function listboxSelectLeg_Callback(hObject, eventdata, handles)
% hObject    handle to listboxSelectLeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxSelectLeg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxSelectLeg


% --- Executes during object creation, after setting all properties.
function listboxSelectLeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxSelectLeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNDimensions_Callback(hObject, eventdata, handles)
% hObject    handle to editNDimensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNDimensions as text
%        str2double(get(hObject,'String')) returns contents of editNDimensions as a double


% --- Executes during object creation, after setting all properties.
function editNDimensions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNDimensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNGMMComponents_Callback(hObject, eventdata, handles)
% hObject    handle to editNGMMComponents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNGMMComponents as text
%        str2double(get(hObject,'String')) returns contents of editNGMMComponents as a double


% --- Executes during object creation, after setting all properties.
function editNGMMComponents_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNGMMComponents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editRescaleFactor_Callback(hObject, eventdata, handles)
% hObject    handle to editRescaleFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRescaleFactor as text
%        str2double(get(hObject,'String')) returns contents of editRescaleFactor as a double


% --- Executes during object creation, after setting all properties.
function editRescaleFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRescaleFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pushbutton1_Callback(hObject, eventdata, handles);
close( gcf );


% --- Executes on selection change in listboxSATExtractionAlgorithms.
function listboxSATExtractionAlgorithms_Callback(hObject, eventdata, handles)
% hObject    handle to listboxSATExtractionAlgorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxSATExtractionAlgorithms contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxSATExtractionAlgorithms


% --- Executes during object creation, after setting all properties.
function listboxSATExtractionAlgorithms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxSATExtractionAlgorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxGVFInput.
function listboxGVFInput_Callback(hObject, eventdata, handles)
% hObject    handle to listboxGVFInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxGVFInput contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxGVFInput


% --- Executes during object creation, after setting all properties.
function listboxGVFInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxGVFInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editLegSelectionThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editLegSelectionThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLegSelectionThreshold as text
%        str2double(get(hObject,'String')) returns contents of editLegSelectionThreshold as a double


% --- Executes during object creation, after setting all properties.
function editLegSelectionThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLegSelectionThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSnakeGaussianSigma_Callback(hObject, eventdata, handles)
% hObject    handle to editSnakeGaussianSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSnakeGaussianSigma as text
%        str2double(get(hObject,'String')) returns contents of editSnakeGaussianSigma as a double


% --- Executes during object creation, after setting all properties.
function editSnakeGaussianSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSnakeGaussianSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editBoneAreaThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editBoneAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBoneAreaThreshold as text
%        str2double(get(hObject,'String')) returns contents of editBoneAreaThreshold as a double


% --- Executes during object creation, after setting all properties.
function editBoneAreaThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBoneAreaThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxCommandLineMode.
function checkboxCommandLineMode_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCommandLineMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCommandLineMode


% --- Executes on button press in checkboxWarpWSandFStoNS.
function checkboxWarpWSandFStoNS_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxWarpWSandFStoNS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxWarpWSandFStoNS
