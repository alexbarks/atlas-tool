function varargout = heat_maps(varargin)
% HEAT_MAPS MATLAB code for heat_maps.fig
%      HEAT_MAPS, by itself, creates a new HEAT_MAPS or raises the existing
%      singleton*.
%
%      H = HEAT_MAPS returns the handle to a new HEAT_MAPS or the handle to
%      the existing singleton*.
%
%      HEAT_MAPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HEAT_MAPS.M with the given input arguments.
%
%      HEAT_MAPS('Property','Value',...) creates a new HEAT_MAPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before heat_maps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to heat_maps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help heat_maps

% Last Modified by GUIDE v2.5 26-Aug-2016 11:33:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @heat_maps_OpeningFcn, ...
                   'gui_OutputFcn',  @heat_maps_OutputFcn, ...
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


% --- Executes just before heat_maps is made visible.
function heat_maps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to heat_maps (see VARARGIN)

% Choose default command line output for heat_maps
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes heat_maps wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = heat_maps_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in WSS_calc.
function WSS_calc_Callback(hObject, eventdata, handles)
% hObject    handle to WSS_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDir = pwd;
choice = questdlg('Do you want to process an individual patient or batch several?', ...
	'WSS calculation', ...
	'Individual','Batch','Individual');
% Handle response
switch choice
    case 'Individual'
        
        [FileName,MrstructPath,FilterIndex] = uigetfile(currDir,'Select the mag_struct file');
        try
            cd(MrstructPath);
        catch
            warndlg('File not found!');
            return;
        end
        [FileName,MimicsSegPath,FilterIndex] = uigetfile([MrstructPath '\..\*'],'Select the Mimics aorta or PA text file, or the already computed _grayvalues_mask_struct mat file');
        try
            cd(MimicsSegPath);
        catch
            warndlg('File not found!');
            return;
        end
        cd(currDir)
        WssFraction = get(handles.WssFraction,'Value');
        WssThresh = get(handles.WssThresh,'Value');
        plotFlag = get(handles.checkbox_plot,'Value');
        saveFlag = get(handles.checkbox_saveHist,'Value');
        WSS_syst = get(handles.WSS_sysTime,'Value');
        WSS_syst_avg = get(handles.WSS_syst_avg,'Value');
        WSS_allTimes = get(handles.WSS_allTimes,'Value');
        intracranial = get(handles.checkbox_intracranial,'Value');
        % TimeFlag: 0 if only peak systole; 1 if average over 5 systolic phases; 2 if
        % all phases
        if WSS_syst == 1
            TimeFlag = 0;
        elseif WSS_syst_avg == 1
            TimeFlag = 1;
        elseif WSS_allTimes == 1
            TimeFlag = 2;
        end
        wss_ensight_Flag = get(handles.checkbox_ensight,'Value');
        % % If WSS calculation for systole only force to 0 >> taken into account in
        % mimics_to_Wss
        % if (WSS_allTimes == 0 && wss_ensight_Flag == 1)
        %     wss_ensight_Flag = 0;
        % end
        hematocritFlag = get(handles.checkbox_hematocrit,'Value');
        multipleMasksFlag = get(handles.checkbox_masks,'Value');
        if(isequal(FileName(end-3:end),'.txt'))
            mimicsFileFlag = 1;
        elseif(isequal(FileName(end-3:end),'.mat'))
            mimicsFileFlag = 0;
            MimicsSegPath=[MimicsSegPath FileName];
        end
        mimics_to_Wss([MrstructPath],[MimicsSegPath],WssFraction,WssThresh,plotFlag,saveFlag,TimeFlag,wss_ensight_Flag,hematocritFlag,multipleMasksFlag,mimicsFileFlag,FileName,intracranial)
        h = msgbox('Wall shear stress calculation done!');
   
    case 'Batch'
        folder_name = uigetdir(currDir,'Select the folder containing all patient individual folders');
        WSS_syst = get(handles.WSS_sysTime,'Value');
        WSS_syst_avg = get(handles.WSS_syst_avg,'Value');
        WSS_allTimes = get(handles.WSS_allTimes,'Value');
        intracranial = get(handles.checkbox_intracranial,'Value');
        % TimeFlag: 0 if only peak systole; 1 if average over 5 systolic phases; 2 if
        % all phases
        if WSS_syst == 1
            TimeFlag = 0;
        elseif WSS_syst_avg == 1
            TimeFlag = 1;
        elseif WSS_allTimes == 1
            TimeFlag = 2;
        end
        wss_batch(folder_name, TimeFlag,intracranial);
        h = msgbox('Wall shear stress calculation done!');
end


% --- Executes on button press in checkbox_plot.
function checkbox_plot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot


% --- Executes on slider movement.
function WssFraction_Callback(hObject, eventdata, handles)
% hObject    handle to WssFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

wss_frac = get(hObject,'Value');
set(handles.wss_value, 'String', wss_frac);


% --- Executes during object creation, after setting all properties.
function WssFraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WssFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function wss_value_Callback(hObject, eventdata, handles)
% hObject    handle to wss_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wss_value as text
%        str2double(get(hObject,'String')) returns contents of wss_value as a double

wss_frac=str2num(get(hObject,'String'));
set(handles.WssFraction, 'Value', wss_frac);


% --- Executes during object creation, after setting all properties.
function wss_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wss_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function WssThresh_Callback(hObject, eventdata, handles)
% hObject    handle to WssThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

wss_threshold = get(hObject,'Value');
set(handles.wssThreshold_value, 'String', wss_threshold);


% --- Executes during object creation, after setting all properties.
function WssThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WssThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function wssThreshold_value_Callback(hObject, eventdata, handles)
% hObject    handle to wssThreshold_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wssThreshold_value as text
%        str2double(get(hObject,'String')) returns contents of wssThreshold_value as a double

wss_threshold=str2num(get(hObject,'String'));
set(handles.WssThresh, 'Value', wss_threshold);


% --- Executes during object creation, after setting all properties.
function wssThreshold_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wssThreshold_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_saveHist.
function checkbox_saveHist_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_saveHist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_saveHist


% --- Executes on button press in checkbox_ensight.
function checkbox_ensight_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ensight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ensight

ensight = get(hObject,'Value');
WSS_sysTime = get(handles.WSS_sysTime, 'Value');
if (ensight == 1 && WSS_sysTime == 1)
    set(hObject,'Value',0);
    warndlg('WSS calculation at every time step is required')
end


% --- Executes on button press in checkbox_hematocrit.
function checkbox_hematocrit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_hematocrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_hematocrit


% --- Executes on button press in checkbox_masks.
function checkbox_masks_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_masks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_masks


% --- Executes on button press in heatMaps_creation.
function heatMaps_creation_Callback(hObject, eventdata, handles)
% hObject    handle to heatMaps_creation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currDir = pwd;
% atlas_folder = '\\10.61.223.37\data_imaging\cv_mri\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls';
peak_systolicFlag = get(handles.signArea_calc,'Value');
if peak_systolicFlag == 0   % average over 5 systolic time points
    atlas_folder = '\\10.61.223.37\data_imaging\cv_mri\Aorta-4D_Flow\Results\Pim\PimsProjectWrapUp\7_AgeMatching\data\control_atlases\5time_averaged\all_controls\atlas.mat';
elseif peak_systolicFlag == 1   % peak systolic time point only
    atlas_folder = '\\10.61.223.37\data_imaging\cv_mri\Aorta-4D_Flow\Results\Pim\PimsProjectWrapUp\7_AgeMatching\data\control_atlases\peakSystolic\all_controls\atlas.mat';
end
[FileName,AtlasPath,FilterIndex] = uigetfile(atlas_folder,'Select the atlas.mat file');  % donner chemin par defaut
try
    cd(AtlasPath);
catch
    warndlg('File not found!');
    return;
end
[FileName,FilePath,FilterIndex] = uigetfile(currDir,'Select the mag_struct file');
try
    cd(FilePath);
catch
    warndlg('File not found!');
    return;
end
cd(currDir)
plotFlag = get(handles.RE_calc,'Value');
calculateIE_Flag = get(handles.IE_calc,'Value');
calculate_area_of_higherlowerFlag = get(handles.velVolume_wssArea_calc,'Value');
images_for_surgeryFlag = get(handles.peakSyst_comp,'Value');
heat_map_traffic_light_scalars_affine_registration(AtlasPath,FilePath(1:end-1),plotFlag,calculateIE_Flag,calculate_area_of_higherlowerFlag,peak_systolicFlag,images_for_surgeryFlag)
h = msgbox('Heat maps creation done');


% --- Executes on button press in RE_calc.
function RE_calc_Callback(hObject, eventdata, handles)
% hObject    handle to RE_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RE_calc


% --- Executes on button press in IE_calc.
function IE_calc_Callback(hObject, eventdata, handles)
% hObject    handle to IE_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IE_calc


% --- Executes on button press in velVolume_wssArea_calc.
function velVolume_wssArea_calc_Callback(hObject, eventdata, handles)
% hObject    handle to velVolume_wssArea_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velVolume_wssArea_calc


% --- Executes on button press in signArea_calc.
function signArea_calc_Callback(hObject, eventdata, handles)
% hObject    handle to signArea_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of signArea_calc


% --- Executes on button press in peakSyst_comp.
function peakSyst_comp_Callback(hObject, eventdata, handles)
% hObject    handle to peakSyst_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of peakSyst_comp


% --- Executes on button press in WSS_visu.
function WSS_visu_Callback(hObject, eventdata, handles)
% hObject    handle to WSS_visu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
WSS_syst = get(handles.WSS_sysTime,'Value');
WSS_syst_avg = get(handles.WSS_syst_avg,'Value');
WSS_allTimes = get(handles.WSS_allTimes,'Value');
% TimeFlag: 0 if only peak systole; 1 if average over 5 systolic phases; 2 if
% all phases
if WSS_syst == 1
    TimeFlag = 0;
elseif WSS_syst_avg == 1
    TimeFlag = 1;
elseif WSS_allTimes == 1
    TimeFlag = 2;
end
view_WSS_vectors(TimeFlag);


% --- Executes on button press in checkbox_intracranial.
function checkbox_intracranial_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_intracranial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_intracranial
