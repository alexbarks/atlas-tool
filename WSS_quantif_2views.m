function varargout = WSS_quantif_2views(varargin)
% WSS_QUANTIF_2VIEWS MATLAB code for WSS_quantif_2views.fig
%      WSS_QUANTIF_2VIEWS, by itself, creates a new WSS_QUANTIF_2VIEWS or raises the existing
%      singleton*.
%
%      H = WSS_QUANTIF_2VIEWS returns the handle to a new WSS_QUANTIF_2VIEWS or the handle to
%      the existing singleton*.
%
%      WSS_QUANTIF_2VIEWS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WSS_QUANTIF_2VIEWS.M with the given input arguments.
%
%      WSS_QUANTIF_2VIEWS('Property','Value',...) creates a new WSS_QUANTIF_2VIEWS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WSS_quantif_2views_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WSS_quantif_2views_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WSS_quantif_2views

% Last Modified by GUIDE v2.5 15-Nov-2016 16:25:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WSS_quantif_2views_OpeningFcn, ...
    'gui_OutputFcn',  @WSS_quantif_2views_OutputFcn, ...
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


% --- Executes just before WSS_quantif_2views is made visible.
function WSS_quantif_2views_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WSS_quantif_2views (see VARARGIN)

% Choose default command line output for WSS_quantif_2views
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WSS_quantif_2views wait for user response (see UIRESUME)
% uiwait(handles.figure1);

currDir = pwd;
[FileName,MrstructPath,FilterIndex] = uigetfile(currDir,'Select the Wss_struct file of your patient');
try
    cd(MrstructPath);
catch
    warndlg('File not found!');
    return;
end

load(FileName);
extName = FileName(12:end);
load(strcat('Wss_point_cloud_',extName));
load(strcat('mask_struct_',extName));
% load Wss_struct_aorta.mat
% load Wss_point_cloud_aorta
load vel_struct
% load mask_struct_aorta

mask2 = mrstruct_mask.dataAy;
mask2_vox = mrstruct_mask.vox;

L2 = (mask2 ~= 0);
contours = zeros(size(L2));
contours(L2==0) = -1;
contours(L2==1) = 1;
[F,V] = isosurface(contours,0); % make a surface from the detected contours
V = V .* (ones(size(V,1),1) * mask2_vox(1:3));

velocity = mrStruct.dataAy; clear mrStruct
for t = 1:size(velocity,5)
    vx = squeeze(velocity(:,:,:,1,t));
    vy = squeeze(velocity(:,:,:,2,t));
    vz = squeeze(velocity(:,:,:,3,t));
    vmagn = sqrt(vx.^2 + vy.^2 + vz.^2);
    mean_velo(t) = mean(vmagn(L2));
end

h_meanVel=figure('Name','Mean velocity');
plot(1:size(velocity,5),mean_velo,'-ko','LineWidth',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',14);
ylabel('Mean velocity (m/s)');
xlabel('Time frame #');

[I,time] = find(mean_velo==max(mean_velo));

WSS_all = Wss_point_cloud; clear Wss_point_cloud

choice = questdlg('Do you want to extract WSS values calculated...', ...
    'WSS regional quantification', ...
    'at peak systole?','while averaging up to 5 systolic timesteps?','at peak systole?');
% Handle response
switch choice
    case 'at peak systole?'
        TimeFlag=0;
    case 'while averaging up to 5 systolic timesteps?'
        TimeFlag=1;
end

if TimeFlag==0
    % Peak systolic WSS
    figure(h_meanVel)
    hold on, plot(time,mean_velo(time),'-ko','LineWidth',4,...
        'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    if size(Wss_matrix,5) > 5
        WSS = Wss_matrix(:,:,:,:,time);        
    elseif size(Wss_matrix,5) == 5
        WSS = Wss_matrix(:,:,:,:,3);
    elseif size(Wss_matrix,5) == 4
        WSS = Wss_matrix(:,:,:,:,2);
    elseif size(Wss_matrix,5) == 3
        WSS = Wss_matrix(:,:,:,:,1);
    elseif size(Wss_matrix,5) == 1    % Emilie: calculated at only one time (peak systole)
        WSS = Wss_matrix(:,:,:,:,1);
    end
    wss_m = squeeze(sum(WSS.^2,4).^0.5);
elseif TimeFlag==1  % Averaged systolic WSS
    if size(Wss_matrix,5) == 1
        warndlg('WSS was previously calculated only at peak systole!');
        return;
        % Velocity averaged over x systolic time frames
    elseif time == 2    % second time frame is peak systole: averaging over 4 timesteps
        disp('AVERAGE OVER 4 TIME FRAMES!')
        data2.x_value_wss_t1 = Wss_matrix(:,:,:,1,time-1);data2.y_value_wss_t1 = Wss_matrix(:,:,:,2,time-1);data2.z_value_wss_t1 = Wss_matrix(:,:,:,3,time-1);
        data2.x_value_wss_t2 = Wss_matrix(:,:,:,1,time);  data2.y_value_wss_t2 = Wss_matrix(:,:,:,2,time);  data2.z_value_wss_t2 = Wss_matrix(:,:,:,3,time);
        data2.x_value_wss_t3 = Wss_matrix(:,:,:,1,time+1);data2.y_value_wss_t3 = Wss_matrix(:,:,:,2,time+1);data2.z_value_wss_t3 = Wss_matrix(:,:,:,3,time+1);
        data2.x_value_wss_t4 = Wss_matrix(:,:,:,1,time+2);data2.y_value_wss_t4 = Wss_matrix(:,:,:,2,time+2);data2.z_value_wss_t4 = Wss_matrix(:,:,:,3,time+2);
        data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4)./4;
        data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4)./4;
        data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4)./4;
        figure(h_meanVel)
        hold on, plot(time-1:time+2,mean_velo(time-1:time+2),'-ko','LineWidth',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    elseif time == 1 % first time frame is peak systole: averaging over 3 timesteps
        disp('AVERAGE OVER 3 TIME FRAMES!')
        data2.x_value_wss_t1 = Wss_matrix(:,:,:,1,time);  data2.y_value_wss_t1 = Wss_matrix(:,:,:,2,time);  data2.z_value_wss_t1 = Wss_matrix(:,:,:,3,time);
        data2.x_value_wss_t2 = Wss_matrix(:,:,:,1,time+1);data2.y_value_wss_t2 = Wss_matrix(:,:,:,2,time+1);data2.z_value_wss_t2 = Wss_matrix(:,:,:,3,time+1);
        data2.x_value_wss_t3 = Wss_matrix(:,:,:,1,time+2);data2.y_value_wss_t3 = Wss_matrix(:,:,:,2,time+2);data2.z_value_wss_t3 = Wss_matrix(:,:,:,3,time+2);
        data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3)./3;
        data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3)./3;
        data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3)./3;
        figure(h_meanVel)
        hold on, plot(time:time+2,mean_velo(time:time+2),'-ko','LineWidth',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    else % timestep > 2 is peak systole: averaging over 5 timesteps
        if size(Wss_matrix,5) > 5
            data2.x_value_wss_t1 = Wss_matrix(:,:,:,1,time-2);data2.y_value_wss_t1 = Wss_matrix(:,:,:,2,time-2);data2.z_value_wss_t1 = Wss_matrix(:,:,:,3,time-2);
            data2.x_value_wss_t2 = Wss_matrix(:,:,:,1,time-1);data2.y_value_wss_t2 = Wss_matrix(:,:,:,2,time-1);data2.z_value_wss_t2 = Wss_matrix(:,:,:,3,time-1);
            data2.x_value_wss_t3 = Wss_matrix(:,:,:,1,time);  data2.y_value_wss_t3 = Wss_matrix(:,:,:,2,time);  data2.z_value_wss_t3 = Wss_matrix(:,:,:,3,time);
            data2.x_value_wss_t4 = Wss_matrix(:,:,:,1,time+1);data2.y_value_wss_t4 = Wss_matrix(:,:,:,2,time+1);data2.z_value_wss_t4 = Wss_matrix(:,:,:,3,time+1);
            data2.x_value_wss_t5 = Wss_matrix(:,:,:,1,time+2);data2.y_value_wss_t5 = Wss_matrix(:,:,:,2,time+2);data2.z_value_wss_t5 = Wss_matrix(:,:,:,3,time+2);
            figure(h_meanVel)
            hold on, plot(time-2:time+2,mean_velo(time-2:time+2),'-ko','LineWidth',4,...
                'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
        elseif size(Wss_matrix,5) == 5
            data2.x_value_wss_t1 = Wss_matrix(:,:,:,1,1);data2.y_value_wss_t1 = Wss_matrix(:,:,:,2,1);data2.z_value_wss_t1 = Wss_matrix(:,:,:,3,1);
            data2.x_value_wss_t2 = Wss_matrix(:,:,:,1,2);data2.y_value_wss_t2 = Wss_matrix(:,:,:,2,2);data2.z_value_wss_t2 = Wss_matrix(:,:,:,3,2);
            data2.x_value_wss_t3 = Wss_matrix(:,:,:,1,3);data2.y_value_wss_t3 = Wss_matrix(:,:,:,2,3);data2.z_value_wss_t3 = Wss_matrix(:,:,:,3,3);
            data2.x_value_wss_t4 = Wss_matrix(:,:,:,1,4);data2.y_value_wss_t4 = Wss_matrix(:,:,:,2,4);data2.z_value_wss_t4 = Wss_matrix(:,:,:,3,4);
            data2.x_value_wss_t5 = Wss_matrix(:,:,:,1,5);data2.y_value_wss_t5 = Wss_matrix(:,:,:,2,5);data2.z_value_wss_t5 = Wss_matrix(:,:,:,3,5);
            figure(h_meanVel)
            hold on, plot(1:5,mean_velo(1:5),'-ko','LineWidth',4,...
                'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
        end
        data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4 + data2.x_value_wss_t5)./5;
        data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4 + data2.y_value_wss_t5)./5;
        data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4 + data2.z_value_wss_t5)./5;
    end
    wss_m = sqrt(data2.x_value_wss.^2 + data2.y_value_wss.^2 + data2.z_value_wss.^2);
end

% Calculate the MIPS
mip1 = squeeze(max(wss_m,[],3));
mip2 = squeeze(max(wss_m,[],2));
mip3 = squeeze(max(wss_m,[],1));

dataSize = size(wss_m);
setappdata(handles.figure1, 'wss', wss_m);
setappdata(handles.figure1, 'dataSize', dataSize);
setappdata(handles.figure1, 'TimeFlag', TimeFlag);

% plot velocity MIP - view 1
axes(handles.axes_proj1)
imagesc(mip1);
set(gca,'dataaspectRatio',mask2_vox(1:3));
colorbar('position',[0.910447761194025 0.03286384976525822 0.036182722749886875 0.9389671361502351]);
axis off
% plot velocity MIP - view 2
axes(handles.axes_proj2)
imagesc(mip2);
set(gca,'dataaspectRatio',mask2_vox([1 3 2]))
axis off
% plot velocity MIP - view 3
axes(handles.axes_proj3);
imagesc(mip3);
set(gca,'dataaspectRatio',mask2_vox([3 2 1]))
axis off
%Set Colormap to have black be the lowest value (zero)
cmap = colormap;
cmap(1,1) = 0;
cmap(1,2) = 0;
cmap(1,3) = 0;
% colormap(cmap)
set(handles.text_colorbar, 'String', 'WSS (Pa)');

% 3D plot
x = V(:,1)/mask2_vox(1);
y = V(:,2)/mask2_vox(2);
z = V(:,3)/mask2_vox(3);
axes(handles.axes_3d), patch('Faces',F,'Vertices',[x y z], ...
    'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',1);
%%%
count = 0;
angles(1) = 0;
load mag_struct
magnitude = flipdim(double(mrStruct.dataAy(:,:,:,3)),3);
magnitude(magnitude == 0) = 3;
magnitude(magnitude == 1) = 3;
magnitude(magnitude == 2) = 3;
hold on
s = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3)/2,magnitude(:,:,size(magnitude,3)/2),'EdgeColor','none');
set(s,'HandleVisibility','off','Visible','off');
colormap(handles.axes_3d,gray)
view([180 -90])
% caxis([0 64]);
aspectRatio = 1./mask2_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))
%camlight(-45,0); lighting phong
%%%%
set(handles.slider_anat3dView, 'Value',size(magnitude,3)/2);
set(handles.slider_anat3dView, 'Max',size(magnitude,3));
set(handles.slider_anat3dView, 'SliderStep',[1/(size(magnitude,3)-1) 10/(size(magnitude,3)-1)]);
%%%
axis equal;axis off; axis ij
% h_cam3d=rotate3d(handles.axes_3d);
% set(h_cam3d,'enable','on');
% set(h_cam3d,'rotateStyle','orbit');

setappdata(handles.figure1, 'ROIcount', 0);
setappdata(handles.figure1, 'magnImages', magnitude);
setappdata(handles.figure1, 'hMagnImages', s);
setappdata(handles.figure1, 'nbRot', count);
setappdata(handles.figure1, 'anglesRot', angles);
set(gcf, 'toolbar', 'figure')


% --- Outputs from this function are returned to the command line.
function varargout = WSS_quantif_2views_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_draw.
function pushbutton_draw_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check which projection
proj = get(handles.popupmenu_proj, 'Value');
axes(eval(['handles.axes_proj' num2str(proj)]))

ROIcount = getappdata(handles.figure1, 'ROIcount');
% if ROIcount >= 8
%     warndlg('Already 8 ROIs!', 'Too many ROIs');
%     return;
% end
if isappdata(handles.figure1, 'hROI')
    h_roi = getappdata(handles.figure1, 'hROI');
end
if isappdata(handles.figure1, 'posROI')
    posROI = getappdata(handles.figure1, 'posROI');
end
dataSize = getappdata(handles.figure1, 'dataSize');

% Polygon and mask
switch proj
    case 1
        h_roi{1,ROIcount+1}  = impoly;
        posROI{1,ROIcount+1} = getPosition(h_roi{1,ROIcount+1});
        xi      = posROI{1,ROIcount+1}(:,1);
        yi      = posROI{1,ROIcount+1}(:,2);
        xmax    = max(xi);
        xmin    = min(xi);
        ymax    = max(yi);
        ymin    = min(yi);
        posROI{2,ROIcount+1} = [1,ymin;dataSize(3),ymin;dataSize(3),ymax;1,ymax];
        posROI{3,ROIcount+1} = [1,xmin;dataSize(3),xmin;dataSize(3),xmax;1,xmax];
        axes(handles.axes_proj2)
        h_roi{2,ROIcount+1} = impoly(gca,posROI{2,ROIcount+1});
        axes(handles.axes_proj3)
        h_roi{3,ROIcount+1} = impoly(gca,posROI{3,ROIcount+1});
    case 2
        h_roi{2,ROIcount+1}  = impoly;
        posROI{2,ROIcount+1}   = getPosition(h_roi{2,ROIcount+1});
        zi      = posROI{2,ROIcount+1}(:,1);
        yi      = posROI{2,ROIcount+1}(:,2);
        zmax    = max(zi);
        zmin    = min(zi);
        ymax    = max(yi);
        ymin    = min(yi);
        posROI{1,ROIcount+1} = [1,ymin;dataSize(2),ymin;dataSize(2),ymax;1,ymax];
        posROI{3,ROIcount+1} = [zmin,1;zmin,dataSize(2);zmax,dataSize(2);zmax,1];
        axes(handles.axes_proj1)
        h_roi{1,ROIcount+1} = impoly(gca,posROI{1,ROIcount+1});
        axes(handles.axes_proj3)
        h_roi{3,ROIcount+1} = impoly(gca,posROI{3,ROIcount+1});
    case 3
        h_roi{3,ROIcount+1}  = impoly;
        posROI{3,ROIcount+1}     = getPosition(h_roi{3,ROIcount+1});
        zi      = posROI{3,ROIcount+1}(:,1);
        xi      = posROI{3,ROIcount+1}(:,2);
        zmax    = max(zi);
        zmin    = min(zi);
        xmax    = max(xi);
        xmin    = min(xi);
        posROI{1,ROIcount+1} = [xmin,1;xmin,dataSize(1);xmax,dataSize(1);xmax,1];
        posROI{2,ROIcount+1} = [zmin,1;zmin,dataSize(1);zmax,dataSize(1);zmax,1];
        axes(handles.axes_proj1)
        h_roi{1,ROIcount+1} = impoly(gca,posROI{1,ROIcount+1});
        axes(handles.axes_proj2)
        h_roi{2,ROIcount+1} = impoly(gca,posROI{2,ROIcount+1});
end

setappdata(handles.figure1, 'ROIcount', ROIcount+1);
setappdata(handles.figure1, 'hROI', h_roi);
setappdata(handles.figure1, 'posROI', posROI);

% Color code the ROI
if ROIcount+1 == 1
    set(handles.text_ROI1,'Visible', 'on');
elseif ROIcount+1 == 2
    setColor(h_roi{1,2}, 'r')
    setColor(h_roi{2,2}, 'r')
    setColor(h_roi{3,2}, 'r')
    set(handles.text_ROI2,'Visible', 'on');
elseif ROIcount+1 == 3
    setColor(h_roi{1,3}, 'g')
    setColor(h_roi{2,3}, 'g')
    setColor(h_roi{3,3}, 'g')
    set(handles.text_ROI3,'Visible', 'on');
elseif ROIcount+1 == 4
%     setColor(h_roi{1,4}, 'c')
%     setColor(h_roi{2,4}, 'c')
%     setColor(h_roi{3,4}, 'c')
    setColor(h_roi{1,4}, 'y')
    setColor(h_roi{2,4}, 'y')
    setColor(h_roi{3,4}, 'y')
    set(handles.text_ROI4,'Visible', 'on');
    set(handles.text_ROI4,'ForegroundColor', 'y');
elseif ROIcount+1 == 5
%     setColor(h_roi{1,5}, 'm')
%     setColor(h_roi{2,5}, 'm')
%     setColor(h_roi{3,5}, 'm')
    set(handles.text_ROI5,'Visible', 'on');
    set(handles.text_ROI5,'ForegroundColor', 'b');
elseif ROIcount+1 == 6
%     setColor(h_roi{1,6}, 'y')
%     setColor(h_roi{2,6}, 'y')
%     setColor(h_roi{3,6}, 'y')
    setColor(h_roi{1,6}, 'r')
    setColor(h_roi{2,6}, 'r')
    setColor(h_roi{3,6}, 'r')
    set(handles.text_ROI6,'Visible', 'on');
    set(handles.text_ROI6,'ForegroundColor', 'r');
elseif ROIcount+1 == 7
%     setColor(h_roi{1,7}, 'w')
%     setColor(h_roi{2,7}, 'w')
%     setColor(h_roi{3,7}, 'w')
    setColor(h_roi{1,7}, 'g')
    setColor(h_roi{2,7}, 'g')
    setColor(h_roi{3,7}, 'g')
    set(handles.text_ROI7,'Visible', 'on');
    set(handles.text_ROI7,'ForegroundColor', 'g');
elseif ROIcount+1 == 8
%     setColor(h_roi{1,8}, [.5 .5 .5])
%     setColor(h_roi{2,8}, [.5 .5 .5])
%     setColor(h_roi{3,8}, [.5 .5 .5])
    setColor(h_roi{1,8}, 'y')
    setColor(h_roi{2,8}, 'y')
    setColor(h_roi{3,8}, 'y')
    set(handles.text_ROI8,'Visible', 'on');
    set(handles.text_ROI8,'ForegroundColor', 'y');
elseif ROIcount+1 == 9
%     setColor(h_roi{1,9}, [1 153/255 0])
%     setColor(h_roi{2,9}, [1 153/255 0])
%     setColor(h_roi{3,9}, [1 153/255 0])
    set(handles.text_ROI9,'Visible', 'on');
    set(handles.text_ROI9,'ForegroundColor', 'b');
elseif ROIcount+1 == 10
%     setColor(h_roi{1,10}, [51/255 153/255 51/255])
%     setColor(h_roi{2,10}, [51/255 153/255 51/255])
%     setColor(h_roi{3,10}, [51/255 153/255 51/255])
    setColor(h_roi{1,10}, 'r')
    setColor(h_roi{2,10}, 'r')
    setColor(h_roi{3,10}, 'r')
    set(handles.text_ROI10,'Visible', 'on');
    set(handles.text_ROI10,'ForegroundColor', 'r');
elseif ROIcount+1 == 11
%     setColor(h_roi{1,11}, [102/255 0 102/255])
%     setColor(h_roi{2,11}, [102/255 0 102/255])
%     setColor(h_roi{3,11}, [102/255 0 102/255])
    setColor(h_roi{1,11}, 'g')
    setColor(h_roi{2,11}, 'g')
    setColor(h_roi{3,11}, 'g')
    set(handles.text_ROI11,'Visible', 'on');
    set(handles.text_ROI11,'ForegroundColor', 'g');
elseif ROIcount+1 == 12
%     setColor(h_roi{1,12}, [0 153/255 153/255])
%     setColor(h_roi{2,12}, [0 153/255 153/255])
%     setColor(h_roi{3,12}, [0 153/255 153/255])
    setColor(h_roi{1,12}, 'y')
    setColor(h_roi{2,12}, 'y')
    setColor(h_roi{3,12}, 'y')
    set(handles.text_ROI12,'Visible', 'on');
    set(handles.text_ROI12,'ForegroundColor', 'y');
end

% Add callbacks for when impolys are adjusted
addNewPositionCallback(h_roi{1,ROIcount+1},@(pos) redrawRoi23(pos, hObject, handles, ROIcount+1));
addNewPositionCallback(h_roi{2,ROIcount+1},@(pos) redrawRoi13(pos, hObject, handles, ROIcount+1));
addNewPositionCallback(h_roi{3,ROIcount+1},@(pos) redrawRoi12(pos, hObject, handles, ROIcount+1));


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ROIcount = getappdata(handles.figure1, 'ROIcount');
if ROIcount~=0
    warndlg('ROIs have already been drawn, please save and reload data','Too many ROIs');
    return;
end

currDir = pwd;
load(strcat(currDir,'\..','\regional_masks\ROIs'));

% Plot
for i = 1:size(posROI,2)
    axes(handles.axes_proj1)
    h_roi{1,i} = impoly(gca, posROI{1,i});
    axes(handles.axes_proj2)
    h_roi{2,i} = impoly(gca, posROI{2,i});
    axes(handles.axes_proj3)
    h_roi{3,i} = impoly(gca, posROI{3,i});
    setappdata(handles.figure1, 'posROI', posROI);
    setappdata(handles.figure1, 'hROI', h_roi);
    setappdata(handles.figure1, 'ROIcount', size(posROI,2));
    % Add Callbacks for loaded ROIs
    addNewPositionCallback(h_roi{1,i},@(pos) redrawRoi23(pos, [], handles, i));
    addNewPositionCallback(h_roi{2,i},@(pos) redrawRoi13(pos, [], handles, i));
    addNewPositionCallback(h_roi{3,i},@(pos) redrawRoi12(pos, [], handles, i));
end

% setappdata(handles.figure1, 'posROI', posROI);
% setappdata(handles.figure1, 'hROI', h_roi);
% setappdata(handles.figure1, 'ROIcount', size(posROI,2));

% Color code loaded ROIs according to colorROI fcn
if size(posROI,2) < 6
    warndlg('ROIs are missing, please draw all ROIs','Missing ROIs');
    return;
else
    set(handles.text_ROI1,'Visible', 'on');
    setColor(h_roi{1,2}, 'r')
    setColor(h_roi{2,2}, 'r')
    setColor(h_roi{3,2}, 'r')
    set(handles.text_ROI2,'Visible', 'on');
    setColor(h_roi{1,3}, 'g')
    setColor(h_roi{2,3}, 'g')
    setColor(h_roi{3,3}, 'g')
    set(handles.text_ROI3,'Visible', 'on');
%     setColor(h_roi{1,4}, 'c')
%     setColor(h_roi{2,4}, 'c')
%     setColor(h_roi{3,4}, 'c')
    setColor(h_roi{1,4}, 'y')
    setColor(h_roi{2,4}, 'y')
    setColor(h_roi{3,4}, 'y')
    set(handles.text_ROI4,'Visible', 'on');
%     setColor(h_roi{1,5}, 'm')
%     setColor(h_roi{2,5}, 'm')
%     setColor(h_roi{3,5}, 'm')
    set(handles.text_ROI5,'Visible', 'on');
%     setColor(h_roi{1,6}, 'y')
%     setColor(h_roi{2,6}, 'y')
%     setColor(h_roi{3,6}, 'y')
    setColor(h_roi{1,6}, 'r')
    setColor(h_roi{2,6}, 'r')
    setColor(h_roi{3,6}, 'r')
    set(handles.text_ROI6,'Visible', 'on');
    if size(posROI,2) == 8
        setColor(h_roi{1,7}, 'w')
        setColor(h_roi{2,7}, 'w')
        setColor(h_roi{3,7}, 'w')
        set(handles.text_ROI7,'Visible', 'on');
        setColor(h_roi{1,8}, [.5 .5 .5])
        setColor(h_roi{2,8}, [.5 .5 .5])
        setColor(h_roi{3,8}, [.5 .5 .5])
        set(handles.text_ROI8,'Visible', 'on');
    elseif size(posROI,2) == 12
        setColor(h_roi{1,7}, 'g')
        setColor(h_roi{2,7}, 'g')
        setColor(h_roi{3,7}, 'g')
        set(handles.text_ROI7,'Visible', 'on');
        setColor(h_roi{1,8}, 'y')
        setColor(h_roi{2,8}, 'y')
        setColor(h_roi{3,8}, 'y')
        set(handles.text_ROI8,'Visible', 'on');
%         setColor(h_roi{1,9}, [1 153/255 0])
%         setColor(h_roi{2,9}, [1 153/255 0])
%         setColor(h_roi{3,9}, [1 153/255 0])
        set(handles.text_ROI9,'Visible', 'on');
%         setColor(h_roi{1,10}, [51/255 153/255 51/255])
%         setColor(h_roi{2,10}, [51/255 153/255 51/255])
%         setColor(h_roi{3,10}, [51/255 153/255 51/255])
        setColor(h_roi{1,10}, 'r')
        setColor(h_roi{2,10}, 'r')
        setColor(h_roi{3,10}, 'r')
        set(handles.text_ROI10,'Visible', 'on');
%         setColor(h_roi{1,11}, [102/255 0 102/255])
%         setColor(h_roi{2,11}, [102/255 0 102/255])
%         setColor(h_roi{3,11}, [102/255 0 102/255])
        setColor(h_roi{1,11}, 'g')
        setColor(h_roi{2,11}, 'g')
        setColor(h_roi{3,11}, 'g')
        set(handles.text_ROI11,'Visible', 'on');
%         setColor(h_roi{1,12}, [0 153/255 153/255])
%         setColor(h_roi{2,12}, [0 153/255 153/255])
%         setColor(h_roi{3,12}, [0 153/255 153/255])
        setColor(h_roi{1,12}, 'y')
        setColor(h_roi{2,12}, 'y')
        setColor(h_roi{3,12}, 'y')
        set(handles.text_ROI12,'Visible', 'on');
    end
end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ROIcount = getappdata(handles.figure1, 'ROIcount');
% if ROIcount < 8
%     warndlg('Less than 8 ROIs! Please draw all 8 ROIs', 'Missing ROIs');
%     return;
% end
wss_m = getappdata(handles.figure1, 'wss');
dataSize = getappdata(handles.figure1, 'dataSize');
TimeFlag = getappdata(handles.figure1, 'TimeFlag');
h_roi = getappdata(handles.figure1, 'hROI');
posROI = getappdata(handles.figure1, 'posROI');

% Excel file lines
indices{1,1} = 'ROI #';
indices{2,1} = 'mean';
indices{3,1} = 'median';
indices{4,1} = 'max';
indices{5,1} = 'min';
indices{6,1} = 'std';
indices{7,1} = 'max5percent';
indices{8,1} = 'max2percent';
indices{9,1} = 'min5percent';
indices{10,1} = 'min2percent';

for k=1:size(h_roi,2)
    roi1 = poly2mask(posROI{1,k}(:,1), posROI{1,k}(:,2),dataSize(1),dataSize(2));
    roi2 = poly2mask(posROI{2,k}(:,1), posROI{2,k}(:,2),dataSize(1),dataSize(3));
    roi3 = poly2mask(posROI{3,k}(:,1), posROI{3,k}(:,2),dataSize(2),dataSize(3));
    for i = 1:dataSize(3)
        voi1(:,:,i) = roi1;
    end
    for i = 1:dataSize(2)
        voi2(:,i,:) = roi2;
    end
    for i = 1:dataSize(1)
        voi3(i,:,:) = roi3;
    end
    voi_all = voi1 + voi2 + voi3;
    voi_all(voi_all < 3) = 0;
    voi_all(voi_all == 3) = 1;
    wss_mask = wss_m.*voi_all;
    wss_mask(wss_mask == 0) = NaN;
    % Compute quantitative indices
    indices{1,k+1} = k;
    indices{2,k+1} = nanmean(wss_mask(:));
    indices{3,k+1} = nanmedian(wss_mask(:));
    indices{4,k+1} = nanmax(wss_mask(:));
    indices{5,k+1} = nanmin(wss_mask(:));
    indices{6,k+1} = nanstd(wss_mask(:));
    WSS_sorted = sort(wss_mask(:));
    WSS_sorted = WSS_sorted(~isnan(WSS_sorted));
    indices{7,k+1} = mean(WSS_sorted(end-5/100*ceil(length(WSS_sorted)):end));
    indices{8,k+1} = mean(WSS_sorted(end-2/100*ceil(length(WSS_sorted)):end));
    indices{9,k+1} = mean(WSS_sorted(1:5/100*ceil(length(WSS_sorted))));
    indices{10,k+1} = mean(WSS_sorted(1:2/100*ceil(length(WSS_sorted))));
end

currDir = pwd;
mkdir(strcat(currDir, '\..'),'regional_masks')
cd(strcat(currDir, '\..','\regional_masks'));
% Save ROIs
save('ROIs', 'posROI');
% Save quantitative indices in an Excel sheet
if TimeFlag==0
    xls_file = 'wss_indices_peakSyst.xls';
elseif TimeFlag==1
    xls_file = 'wss_indices_systAvg.xls';
end
xlswrite(xls_file,indices);
% Save figures
A1 = (handles.axes_proj1);
A2 = (handles.axes_proj2);
A3 = (handles.axes_proj3);
T1 = handles.text_ROI1;
T2 = handles.text_ROI2;
T3 = handles.text_ROI3;
T4 = handles.text_ROI4;
T5 = handles.text_ROI5;
T6 = handles.text_ROI6;
T7 = handles.text_ROI7;
T8 = handles.text_ROI8;
F1=figure; %('Visible', 'off');
set(gcf, 'color', [.5 .5 .5]);
set(gcf, 'InvertHardCopy', 'off');
c2 = copyobj(A1,F1);
c3 = copyobj(T1,F1);
c4 = copyobj(T2,F1);
c5 = copyobj(T3,F1);
c6 = copyobj(T4,F1);
c7 = copyobj(T5,F1);
c8 = copyobj(T6,F1);
c9 = copyobj(T7,F1);
c10 = copyobj(T8,F1);
% set(gca,'pos',[0 0 1 1]);
c1 =colorbar;
cpos = get(c2,'Position');
cpos(2) = 2;
% cpos(2) = .125;
% cpos(4) = .75;
set(c2,'Position',cpos);
cpos = get(c10,'Position');
cpos(2) = 2;
set(c10,'Position',cpos);
cpos = get(c9,'Position');
cpos(2) = 4;
set(c9,'Position',cpos);
cpos = get(c8,'Position');
cpos(2) = 6;
set(c8,'Position',cpos);
cpos = get(c7,'Position');
cpos(2) = 8;
set(c7,'Position',cpos);
cpos = get(c6,'Position');
cpos(2) = 10;
set(c6,'Position',cpos);
cpos = get(c5,'Position');
cpos(2) = 12;
set(c5,'Position',cpos);
cpos = get(c4,'Position');
cpos(2) = 14;
set(c4,'Position',cpos);
cpos = get(c3,'Position');
cpos(2) = 16;
set(c3,'Position',cpos);
set(c1,'YColor', 'white');
cmap = colormap;
cmap(1,1) = 0;
cmap(1,2) = 0;
cmap(1,3) = 0;
colormap(cmap)
% print(F1,'-djpeg','-noui', '-r250', 'projection1.jpg');
print(F1,'-djpeg', '-r250', 'projection1.jpg');
close(F1)

F2=figure; %('Visible', 'off');
set(gcf, 'color', [.5 .5 .5]);
set(gcf, 'InvertHardCopy', 'off');
c2 = copyobj(A2,F2);
c3 = copyobj(T1,F2);
c4 = copyobj(T2,F2);
c5 = copyobj(T3,F2);
c6 = copyobj(T4,F2);
c7 = copyobj(T5,F2);
c8 = copyobj(T6,F2);
c9 = copyobj(T7,F2);
c10 = copyobj(T8,F2);
c1 =colorbar;
cpos = get(c2,'Position');
cpos(2) = 2;
cpos(1) = 20;
set(c2,'Position',cpos);
cpos = get(c10,'Position');
cpos(2) = 2;
set(c10,'Position',cpos);
cpos = get(c9,'Position');
cpos(2) = 4;
set(c9,'Position',cpos);
cpos = get(c8,'Position');
cpos(2) = 6;
set(c8,'Position',cpos);
cpos = get(c7,'Position');
cpos(2) = 8;
set(c7,'Position',cpos);
cpos = get(c6,'Position');
cpos(2) = 10;
set(c6,'Position',cpos);
cpos = get(c5,'Position');
cpos(2) = 12;
set(c5,'Position',cpos);
cpos = get(c4,'Position');
cpos(2) = 14;
set(c4,'Position',cpos);
cpos = get(c3,'Position');
cpos(2) = 16;
set(c3,'Position',cpos);
set(c1,'YColor', 'white');
cmap = colormap;
cmap(1,1) = 0;
cmap(1,2) = 0;
cmap(1,3) = 0;
colormap(cmap)
print(F2,'-djpeg', '-r500', 'projection2.jpg');
close(F2)

F3=figure; %('Visible', 'off');
set(gcf, 'color', [.5 .5 .5]);
set(gcf, 'InvertHardCopy', 'off');
c2 = copyobj(A3,F3);
c3 = copyobj(T1,F3);
c4 = copyobj(T2,F3);
c5 = copyobj(T3,F3);
c6 = copyobj(T4,F3);
c7 = copyobj(T5,F3);
c8 = copyobj(T6,F3);
c9 = copyobj(T7,F3);
c10 = copyobj(T8,F3);
c1 =colorbar;
cpos = get(c2,'Position');
cpos(2) = 2;
set(c2,'Position',cpos);
cpos = get(c10,'Position');
cpos(2) = 2;
set(c10,'Position',cpos);
cpos = get(c9,'Position');
cpos(2) = 4;
set(c9,'Position',cpos);
cpos = get(c8,'Position');
cpos(2) = 6;
set(c8,'Position',cpos);
cpos = get(c7,'Position');
cpos(2) = 8;
set(c7,'Position',cpos);
cpos = get(c6,'Position');
cpos(2) = 10;
set(c6,'Position',cpos);
cpos = get(c5,'Position');
cpos(2) = 12;
set(c5,'Position',cpos);
cpos = get(c4,'Position');
cpos(2) = 14;
set(c4,'Position',cpos);
cpos = get(c3,'Position');
cpos(2) = 16;
set(c3,'Position',cpos);
set(c1,'YColor', 'white');
cmap = colormap;
cmap(1,1) = 0;
cmap(1,2) = 0;
cmap(1,3) = 0;
colormap(cmap)
print(F3,'-djpeg', '-r500', 'projection3.jpg');
close(F3)

cd(currDir)
h1 = msgbox('WSS indices were calculated and saved in the regional_masks folder');


% --- Executes on selection change in popupmenu_proj.
function popupmenu_proj_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_proj contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_proj


% --- Executes during object creation, after setting all properties.
function popupmenu_proj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_proj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function redrawRoi23(pos, hObject, handles, k)
% Redraws ROI 2 and 3 when 1 is changed

% Get handles
h_roi = getappdata(handles.figure1, 'hROI');
posROI = getappdata(handles.figure1, 'posROI');

posROI{1,k} = pos;

% Get new pos
xi      = posROI{1,k}(:,1);
yi      = posROI{1,k}(:,2);
xmax    = max(xi);
xmin    = min(xi);
ymax    = max(yi);
ymin    = min(yi);

% Check if bounds have changed
% ROI reduced
for i = 1:size(posROI{2,k}(:,2))
    if ymin > posROI{2,k}(i,2)
        posROI{2,k}(i,2) = ymin;
    end
    if ymax < posROI{2,k}(i,2)
        posROI{2,k}(i,2) = ymax;
    end
end
for i = 1:size(posROI{3,k}(:,2))
    if xmin > posROI{3,k}(i,2)
        posROI{3,k}(i,2) = xmin;
    end
    if xmax < posROI{3,k}(i,2)
        posROI{3,k}(i,2) = xmax;
    end
end

% ROI enlarged
temp = posROI{2,k};
for i = 1:2
    [mintemp coord] = min(temp(:,2)); %min
    if ymin < mintemp
        posROI{2,k}(coord,2) = ymin;
        temp(coord,2) = NaN;
    end
end
temp = posROI{2,k};
for i = 1:2
    [maxtemp coord] = max(temp(:,2)); %max
    if ymax > maxtemp
        posROI{2,k}(coord,2) = ymax;
        temp(coord,2) = NaN;
    end
end
temp = posROI{3,k};
for i = 1:2
    [mintemp coord] = min(temp(:,2)); %min
    if xmin < mintemp
        posROI{3,k}(coord,2) = xmin;
        temp(coord,2) = NaN;
    end
end
temp = posROI{3,k};
for i = 1:2
    [maxtemp coord] = max(temp(:,2)); %max
    if xmax > maxtemp
        posROI{3,k}(coord,2) = xmax;
        temp(coord,2) = NaN;
    end
end

% Redraw impolys
setPosition(h_roi{2,k},posROI{2,k});
setPosition(h_roi{3,k},posROI{3,k});

% Set handles
setappdata(handles.figure1, 'hROI', h_roi);
setappdata(handles.figure1, 'posROI', posROI);

function redrawRoi13(pos, hObject, handles, k)
% Redraws ROI 1 and 3 when 2 is changed

% Get handles
h_roi = getappdata(handles.figure1, 'hROI');
posROI = getappdata(handles.figure1, 'posROI');

posROI{2,k} = pos;

% Get new pos
zi      = posROI{2,k}(:,1);
yi      = posROI{2,k}(:,2);
zmax    = max(zi);
zmin    = min(zi);
ymax    = max(yi);
ymin    = min(yi);

% Check if bounds have changed
% ROI reduced
for i = 1:size(posROI{1,k}(:,2))
    if ymin > posROI{1,k}(i,2)
        posROI{1,k}(i,2) = ymin;
    end
    if ymax < posROI{1,k}(i,2)
        posROI{1,k}(i,2) = ymax;
    end
end
for i = 1:size(posROI{3,k}(:,1))
    if zmin > posROI{3,k}(i,1)
        posROI{3,k}(i,1) = zmin;
    end
    if zmax < posROI{3,k}(i,1)
        posROI{3,k}(i,1) = zmax;
    end
end

% ROI enlarged
temp = posROI{1,k};
for i = 1:2
    [mintemp coord] = min(temp(:,2)); %min
    if ymin < mintemp
        posROI{1,k}(coord,2) = ymin;
        temp(coord,2) = NaN;
    end
end
temp = posROI{1,k};
for i = 1:2
    [maxtemp coord] = max(temp(:,2)); %max
    if ymax > maxtemp
        posROI{1,k}(coord,2) = ymax;
        temp(coord,2) = NaN;
    end
end
temp = posROI{3,k};
for i = 1:2
    [mintemp coord] = min(temp(:,1)); %min
    if zmin < mintemp
        posROI{3,k}(coord,1) = zmin;
        temp(coord,1) = NaN;
    end
end
temp = posROI{3,k};
for i = 1:2
    [maxtemp coord] = max(temp(:,1)); %max
    if zmax > maxtemp
        posROI{3,k}(coord,1) = zmax;
        temp(coord,1) = NaN;
    end
end

% Redraw impolys
setPosition(h_roi{1,k},posROI{1,k});
setPosition(h_roi{3,k},posROI{3,k});

% Set handles
setappdata(handles.figure1, 'hROI', h_roi);
setappdata(handles.figure1, 'posROI', posROI);

function redrawRoi12(pos, hObject, handles, k)
% Redraws ROI 1 and 2 when 3 is changed

% Get handles
h_roi = getappdata(handles.figure1, 'hROI');
posROI = getappdata(handles.figure1, 'posROI');

posROI{3,k} = pos;

% Get new pos
zi      = posROI{3,k}(:,1);
xi      = posROI{3,k}(:,2);
zmax    = max(zi);
zmin    = min(zi);
xmax    = max(xi);
xmin    = min(xi);

% Check if bounds have changed
% ROI reduced
for i = 1:size(posROI{1,k}(:,1))
    if xmin > posROI{1,k}(i,1)
        posROI{1,k}(i,1) = xmin;
    end
    if xmax < posROI{1,k}(i,1)
        posROI{1,k}(i,1) = xmax;
    end
end
for i = 1:size(posROI{2,k}(:,1))
    if zmin > posROI{2,k}(i,1)
        posROI{2,k}(i,1) = zmin;
    end
    if zmax < posROI{2,k}(i,1)
        posROI{2,k}(i,1) = zmax;
    end
end

% ROI enlarged
temp = posROI{1,k};
for i = 1:2
    [mintemp coord] = min(temp(:,1)); %min
    if xmin < mintemp
        posROI{1,k}(coord,1) = xmin;
        temp(coord,1) = NaN;
    end
end
temp = posROI{1,k};
for i = 1:2
    [maxtemp coord] = max(temp(:,1)); %max
    if xmax > maxtemp
        posROI{1,k}(coord,1) = xmax;
        temp(coord,1) = NaN;
    end
end
temp = posROI{2,k};
for i = 1:2
    [mintemp coord] = min(temp(:,1)); %min
    if zmin < mintemp
        posROI{2,k}(coord,1) = zmin;
        temp(coord,1) = NaN;
    end
end
temp = posROI{2,k};
for i = 1:2
    [maxtemp coord] = max(temp(:,1)); %max
    if zmax > maxtemp
        posROI{2,k}(coord,1) = zmax;
        temp(coord,1) = NaN;
    end
end

% Redraw impolys
setPosition(h_roi{1,k},posROI{1,k});
setPosition(h_roi{2,k},posROI{2,k});

% Set handles
setappdata(handles.figure1, 'hROI', h_roi);
setappdata(handles.figure1, 'posROI', posROI);


% --- Executes on slider movement.
function slider_3dView_Callback(hObject, eventdata, handles)
% hObject    handle to slider_3dView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

count = getappdata(handles.figure1, 'nbRot');
angles = getappdata(handles.figure1, 'anglesRot');

count = count + 1;
dtheta = get(handles.slider_3dView,'Value');
dphi = 0;
if count == 1;
    dtheta2 = dtheta*3;
else
    dtheta2 = (dtheta - angles(count-1))*3;
end
axes(handles.axes_3d);
camorbit(dtheta2,dphi,'data',[0 1 0])
angles(count) = dtheta;

setappdata(handles.figure1, 'nbRot', count);
setappdata(handles.figure1, 'anglesRot', angles);


% --- Executes during object creation, after setting all properties.
function slider_3dView_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_3dView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox_anat3dView.
function checkbox_anat3dView_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_anat3dView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_anat3dView

show = round(get(handles.checkbox_anat3dView,'Value'));
s = getappdata(handles.figure1, 'hMagnImages');
if (show == 1)
    set(handles.slider_anat3dView, 'Enable', 'on');
    patchobj = findobj(s);
    set(patchobj,'HandleVisibility','on','Visible','on');
elseif (show == 0)
    set(handles.slider_anat3dView, 'Enable', 'off');
    patchobj = findobj(s);
    set(patchobj,'HandleVisibility','off','Visible','off');
end


% --- Executes on slider movement.
function slider_anat3dView_Callback(hObject, eventdata, handles)
% hObject    handle to slider_anat3dView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

s = getappdata(handles.figure1, 'hMagnImages');
magnitude = getappdata(handles.figure1, 'magnImages');
sliceobj = findobj(s);
delete(sliceobj)
slice = round(get(handles.slider_anat3dView,'Value'));
s = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]).*(size(magnitude,3)-(slice-1)),magnitude(:,:,slice),'EdgeColor','none');
setappdata(handles.figure1, 'hMagnImages', s);


% --- Executes during object creation, after setting all properties.
function slider_anat3dView_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_anat3dView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox_colorbar.
function checkbox_colorbar_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_colorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_colorbar

colbar = round(get(handles.checkbox_colorbar,'Value'));
if (colbar == 1)
    colormap(handles.axes_proj1,jet)
    cmap = colormap;
    cmap(1,1) = 0;
    cmap(1,2) = 0;
    cmap(1,3) = 0;
    colormap(cmap)
elseif (colbar == 0)
    colormap(handles.axes_proj1,gray)
end
