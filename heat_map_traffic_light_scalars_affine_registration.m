function [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,PATHNAME,plotFlag,calculateIE_Flag,calculate_area_of_higherlowerFlag,peak_systolicFlag,images_for_surgeryFlag)

%%% [heat_map,traffic_light]=heat_map_traffic_light_scalars_affine_registration(offset,plotFlag,calculateRE_Flag,calculateIE_Flag,calculate_area_of_significanceFlag,peak_systolicFlag)
%
% This function creates the traffic light figure for velocity and the heat map for WSS for individual patients by comparison with the velocity
% and WSS atlases (ensemble-averaged velocity and WSS maps).
%
% The method for the WSS heat map was published in (so please cite this paper when using this method):
%
% Characterization of abnormal wall shear stress using 4D flow MRI in human bicuspid aortopathy
% van Ooij P, Potters WV, Collins J, Carr M, Carr J, Malaisrie SC, Fedak PWM, McCarthy PM, Markl M, Barker AJ
% Ann Biomed Eng. Accepted 2014 August 8
%
% And the method for the traffic light map hopefuly will be published in (so please cite this paper when using this method):
%
% Brief report: Visualization of Aortic Outflow Patterns for Different Bicuspid Valve Fusion Types
% van Ooij P, Collins J, Carr J, Malaisrie SC, McCarthy PM, Markl M, Barker AJ
% New England Journal of Medicine, 2014
%
% First, the individual patient is co-registered to the ensemble-averaged geometry (created with 'make_geometry_point_cloud.m') and the
% velocity and WSS values of the ensemble-averaged maps (created with 'make_atlas_point_cloud_scalars_affine_registration.m') are interpolated
% to the patient geometry. A comparison between the ensemble-averaged velocity/WSS values and individual velocity/WSS values is carried out.
% Velocities of the individual patient that are higher than mean+2 standard devations (2SD), mean+1 SD and lower than mean_1SD than the
% ensemble-averaged velocities are visualized in red, yellow and green, respectively (hence "traffic light map"). WSS values of the individual
% patient that are higher than mean + 2SD and lower than mean - 2SD are visualized in red and blue, respectively (hence "heat map").
%
% 2014, Pim van Ooij, Northwestern University
%
% Input
% 1)AtlasPath        : The path where the 'atlas.mat' is located.
% 2)MrstructPath     : The path where 'mag_struct' and 'vel_struct.mat' are located.
% 3)calculateRE_Flag: When switched on the registration error will be calculated. However, this is a different RE than the one in the paper
%                     mentioned above as the RE in this function is calculated from AFFINE registration whereas the RE reported in the paper
%                     is calculated from RIGID registration in the function 'make_geometry_point_cloud.m'. See: van Ooij et al. Magn Res Med 2014
% 4)calculateIE_Flag: When switched on the interpolation error (RE, see paper mentioned above) will be calculated. Note that ROIs are needed
%                     which can be drawn manually when switched on. See: van Ooij et al. Magn Res Med 2014
% 5)calculate_velvolume_and_WSSarea_total: The volumes of the red, yellow and green volumes in the traffic light and the surface areas in the heat map are
%                     printed to the screen
% 6)calculate_area_of_significanceFlag: When switched on the area of significance will be calculated (see first paper mentioned above). Note that ROIs are
%                     needed which can be drawn manually when switched on. However, if calculateIE_Flag is switched on than you need to do this only once.
% 7)peak_systolicFlag: When switched on the comparison of the individual patient with the ensemble-averaged maps will be performed for the peak
%                      systolic time frame of the indivual patient only. Note that it makes sense that peak systolic ensemble-averaged are loaded in
%                      as well. Default is the atlas for 5 systolic timesteps averaged.
%
% Output
% 1)traffic_light:     The traffic light map for velocity is saved to this struct
% 2) heat_map:         The heat map for WSS is saved to this struct
%
% Usage
% This code is for creating velocity traffic light maps and WSS heat maps from 'data_done' structs for individual patients as created by
% Pims_postprocessing tool.
% The function for creating traffic light and heat maps from mrStructs is under construction
%
% Examples:
% [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,MrstructPath,plotFlag,calculateRE_Flag,calculateIE_Flag,calculate_velvolume_and_WSSarea_total,calculate_area_of_higherlowerFlag,peak_systolicFlag)
% [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration('','',0,0,0,0,0,0,0,0)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create, or lookup default path cache for flirt and cygwin (in c:\temp)
% find if c:\temp exists, if exist look for cache file
if (exist('c:\temp','dir')==7 && exist('c:\temp\atlas_tool.cfg','file')==2)
    %read settings
    fid = fopen('c:\temp\atlas_tool.cfg'); % need to assign 'w' if want to write to this file
    path_flirt  = fgetl(fid);
    path_cygwin = fgetl(fid);
    fclose(fid);
    clear fid
elseif (exist('c:\temp','dir')==0 || exist('c:\temp\atlas_tool.cfg','file')==0) %else create config
    % make c:temp dir, turn off warning if already exists
    wid = 'MATLAB:MKDIR:DirectoryExists';
    warning('off',wid)
    mkdir('c:\temp')
    warning('on',wid)
    % get working directory and create cfg file with path
    path_flirt = uigetdir('c:\temp','Select your working directory for flirt');
    path_cygwin = uigetdir('c:\temp','Select your working directory for cygwin');
    if ischar(path_flirt) || ischar(path_cygwin)
        i_tmp = (path_flirt=='\'); %repalce control char backslash with slash (in order to be able to write the path)
        path_flirt(i_tmp) = '/';
        i_tmp = (path_cygwin=='\'); %repalce control char backslash with slash (in order to be able to write the path)
        path_cygwin(i_tmp) = '/';
        fid = fopen('c:\temp\atlas_tool.cfg','w');
        fprintf(fid,[path_flirt '\r' path_cygwin]);
        fclose(fid);
        cd(path_flirt)
    end
    clear i_tmp fid working_path wid
end

if nargin < 3
    AtlasPath = '';
    PATHNAME = '';
end

if ~exist(AtlasPath) == 2 || isempty(AtlasPath)
   [FILENAME_atlas,AtlasPath] = uigetfile('C:\1_Chicago\Data\MIMICS\BAVcohorts\control_atlases\41_50','Load atlas.mat');
   %FILENAME_atlas = 'atlas.mat' 
   %AtlasPath = 'C:\1_Chicago\Data\MIMICS\age_matching\controls\19_30'
else
    FILENAME_atlas = 'atlas.mat';
end

if ~exist(PATHNAME) == 2 || isempty(PATHNAME)
    [PATHNAME] = uigetdir('C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\1_PT18_SE\mrstruct','Select patient folder containing mrstruct folder');
  %  PATHNAME = 'X:\cv_mri\Aorta-4D_Flow\BAV\PT271-MJ\3dpc\mrstruct'
    FILENAME1 = 'mask_struct_aorta';        % 1: Load mask
    FILENAME2 = 'vel_struct';               % 2: Load velocity
    FILENAME3 = 'Wss_point_cloud_aorta';    % 3: Load WSS
    FILENAME4 = 'mag_struct';
    MrstructPath = PATHNAME;
else
    MrstructPath = PATHNAME;
    ind_sep=findstr(MrstructPath,'\');
    PATHNAME = MrstructPath(1:ind_sep(end)-1);
    FILENAME2 = 'vel_struct';               % 2: Load velocity
    FILENAME4 = 'mag_struct';
    tempDir=pwd;
    cd(MrstructPath)
    maskStructFiles=ls('mask_struct_*');
    if size(maskStructFiles,1)>1
        [FileName,FilePath,FilterIndex] = uigetfile(MrstructPath,'There are several ''mask_struct'' files, please select the right one');
        FILENAME1 = FileName;        % 1: Load mask
    else
        FILENAME1 = maskStructFiles(1,:);        % 1: Load mask
    end
    WssPtCloudFiles=ls('Wss_point_cloud_*');
    if size(WssPtCloudFiles,1)>1
        [FileName,FilePath,FilterIndex] = uigetfile(MrstructPath,'There are several ''Wss_point_cloud'' files, please select the right one');
        FILENAME3 = FileName;    % 3: Load WSS
    else
        FILENAME3 = WssPtCloudFiles(1,:);        % 1: Load WSS
    end
    cd(tempDir)
end

if nargin < 3 || isempty(plotFlag)
    plotFlag = 1;
end

if nargin < 4 || isempty(calculateIE_Flag)
    calculateIE_Flag = 0;
end

if nargin < 5 || isempty(calculate_area_of_higherlowerFlag)
    calculate_area_of_higherlowerFlag = 1;
end

if nargin < 6 || isempty(peak_systolicFlag)
    peak_systolicFlag = 1;
end

if nargin < 7 || isempty(images_for_surgeryFlag)
    images_for_surgeryFlag = 0;
end

global mrstruct_mask
global mrStruct
global mrstruct % this is to correct line 294 error, i.e. diff names for mrstruct_mask (TODO: don't declare lots of global vars) 
global Wss_point_cloud
global atlas
global angles
global count3

cd([MrstructPath '\..'])

Rotation_Translation = [];

load(strcat(AtlasPath,'\',FILENAME_atlas))
mask1 = atlas.mask;
gray_colormap = colormap(gray);

if plotFlag == 1
    
    atlas_matrix = zeros(size(atlas.mask));
    L = (atlas.mask~=0);
    atlas_matrix(L) = atlas.mean_vel;
    figure('Name', 'Velocity and WSS atlas')
    subplot(2,2,1); L_figure = (squeeze(max(atlas_matrix,[],3))~=0);
    imagesc(flipdim(squeeze(max(atlas_matrix,[],3)),2),'Alphadata',double(flipdim(L_figure,2)));
    axis tight;axis equal;axis ij;axis off;caxis([0 1.5]);
    title('mean velocity atlas (m/s)')

    atlas_matrix(L) = atlas.std_vel;
    subplot(2,2,2); L_figure = (squeeze(max(atlas_matrix,[],3))~=0);
    imagesc(flipdim(squeeze(max(atlas_matrix,[],3)),2),'Alphadata',double(flipdim(L_figure,2)));
    axis tight;axis equal;axis ij;axis off;caxis([0 1.5]);
    title('std velocity atlas (m/s)')
    h_cb = colorbar;
    pos_cb = get(h_cb,'Position');
    set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])
    
    subplot(2,2,3);patch('Faces',atlas.faces,'Vertices',atlas.vertices,'EdgeColor','none', 'FaceVertexCData',atlas.mean_wss,'FaceColor','interp','FaceAlpha',1);
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
    title('mean WSS atlas (Pa)')
    subplot(2,2,4);patch('Faces',atlas.faces,'Vertices',atlas.vertices,'EdgeColor','none', 'FaceVertexCData',atlas.std_wss,'FaceColor','interp','FaceAlpha',1);
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
    title('std WSS atlas (Pa)')
    h_cb = colorbar;
    pos_cb = get(h_cb,'Position');
    set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])
    
end

if calculateIE_Flag == 1;
    
    if ~exist(strcat(AtlasPath,'\interpolation_error_ROI\mask1.mat'),'file')
        
        F1=figure('Name','Atlas shape: Paused after finishing a region so press space when finished!');
        patch('Faces',atlas.faces,'Vertices',atlas.vertices,'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.25);
        view([-180 -90]);axis ij;axis equal;axis off
        
        for i = 1:3
            %Polygon and mask for AAo
            polyAAo = impoly;
            wait(polyAAo);
            region = getPosition(polyAAo);
            
            %          disp('saving, pausing')
            mkdir(AtlasPath,'\interpolation_error_ROI')
            save(strcat([AtlasPath '\interpolation_error_ROI\mask' num2str(i)]),'region');
            pause
        end
        
        close(F1)
    end
    
    load(strcat(AtlasPath,'\interpolation_error_ROI\mask1'))
    atlas_mask_AAo_inner_vel = inpolygon(atlas.x_coor_vel, atlas.y_coor_vel, region(:,1), region(:,2));
    mean_vel_asc_inner = mean(atlas.mean_vel(atlas_mask_AAo_inner_vel));
    atlas_mask_AAo_inner_wss = inpolygon(atlas.x_coor_wss, atlas.y_coor_wss, region(:,1), region(:,2));
    mean_wss_asc_inner = mean(atlas.mean_wss(atlas_mask_AAo_inner_wss));
    load(strcat(AtlasPath,'\interpolation_error_ROI\mask2'))
    atlas_mask_AAo_outer_vel = inpolygon(atlas.x_coor_vel, atlas.y_coor_vel, region(:,1), region(:,2));
    mean_vel_asc_outer = mean(atlas.mean_vel(atlas_mask_AAo_outer_vel));
    atlas_mask_AAo_outer_wss = inpolygon(atlas.x_coor_wss, atlas.y_coor_wss, region(:,1), region(:,2));
    mean_wss_asc_outer = mean(atlas.mean_wss(atlas_mask_AAo_outer_wss));
    load(strcat(AtlasPath,'\interpolation_error_ROI\mask3'))
    atlas_mask_arch_inner_vel = inpolygon(atlas.x_coor_vel, atlas.y_coor_vel, region(:,1), region(:,2));
    mean_vel_arch_inner = mean(atlas.mean_vel(atlas_mask_arch_inner_vel));
    atlas_mask_arch_inner_wss = inpolygon(atlas.x_coor_wss, atlas.y_coor_wss, region(:,1), region(:,2));
    mean_wss_arch_inner = mean(atlas.mean_wss(atlas_mask_arch_inner_wss));

    mean_vel_atlas_total_before_interpolation = mean(atlas.mean_vel)
    mean_wss_atlas_total_before_interpolation = mean(atlas.mean_wss)
    
    mean_vel_before_interpolation(1,1) = mean_vel_asc_inner;
    mean_vel_before_interpolation(2,1) = mean_vel_asc_outer;
    mean_vel_before_interpolation(3,1) = mean_vel_arch_inner
    
    mean_wss_before_interpolation(1,1) = mean_wss_asc_inner;
    mean_wss_before_interpolation(2,1) = mean_wss_asc_outer;
    mean_wss_before_interpolation(3,1) = mean_wss_arch_inner
    
    if plotFlag == 1
        figure('Name','Velocity before interpolation: inner AAo')
        scatter3(atlas.x_coor_vel(atlas_mask_AAo_inner_vel),atlas.y_coor_vel(atlas_mask_AAo_inner_vel),atlas.z_coor_vel(atlas_mask_AAo_inner_vel),20,atlas.mean_vel(atlas_mask_AAo_inner_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);caxis([0 1.5])
        figure('Name','Velocity before interpolation: outer AAo')
        scatter3(atlas.x_coor_vel(atlas_mask_AAo_outer_vel),atlas.y_coor_vel(atlas_mask_AAo_outer_vel),atlas.z_coor_vel(atlas_mask_AAo_outer_vel),20,atlas.mean_vel(atlas_mask_AAo_outer_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);caxis([0 1.5])
        figure('Name','Velocity before interpolation: inner arch')
        scatter3(atlas.x_coor_vel(atlas_mask_arch_inner_vel),atlas.y_coor_vel(atlas_mask_arch_inner_vel),atlas.z_coor_vel(atlas_mask_arch_inner_vel),20,atlas.mean_vel(atlas_mask_arch_inner_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);caxis([0 1.5])

        figure('Name','WSS before interpolation: inner AAo')
        scatter3(atlas.x_coor_wss(atlas_mask_AAo_inner_wss),atlas.y_coor_wss(atlas_mask_AAo_inner_wss),atlas.z_coor_wss(atlas_mask_AAo_inner_wss),20,atlas.mean_wss(atlas_mask_AAo_inner_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);caxis([0 1.5])
        figure('Name','WSS before interpolation: outer AAo')
        scatter3(atlas.x_coor_wss(atlas_mask_AAo_outer_wss),atlas.y_coor_wss(atlas_mask_AAo_outer_wss),atlas.z_coor_wss(atlas_mask_AAo_outer_wss),20,atlas.mean_wss(atlas_mask_AAo_outer_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);caxis([0 1.5])
        figure('Name','WSS before interpolation: inner arch')
        scatter3(atlas.x_coor_wss(atlas_mask_arch_inner_wss),atlas.y_coor_wss(atlas_mask_arch_inner_wss),atlas.z_coor_wss(atlas_mask_arch_inner_wss),20,atlas.mean_wss(atlas_mask_arch_inner_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);caxis([0 1.5])

    end
end

load(strcat(MrstructPath,'\',FILENAME1)) % ugly handling of mrStruct variants (TODO: clean up)
if ~isempty(mrstruct_mask)
    mask2     = mrstruct_mask.dataAy;
    mask2_vox = mrstruct_mask.vox;
elseif ~isempty(mrStruct)
    mask2     = mrStruct.dataAy;
    mask2_vox = mrStruct.vox;
else
    mask2     = mrstruct.dataAy; % lower case 'S'
    mask2_vox = mrstruct.vox;
end
clear mrstruct_mask mrStruct mrstruct

L2 = (mask2 ~= 0);
% create velocity coordinates
[x,y,z] = meshgrid((1:size(mask2,2)).* mask2_vox(2), ...
    (1:size(mask2,1)).*mask2_vox(1),(1:size(mask2,3)).* mask2_vox(3));

data2.x_coor_vel = x(L2);data2.y_coor_vel = y(L2);data2.z_coor_vel = z(L2);
contours = zeros(size(L2));
contours(L2==0) = -1;
contours(L2==1) = 1;
[F,V] = isosurface(contours,0); % make a surface from the detected contours
V = V .* (ones(size(V,1),1) * mask2_vox(1:3));
[data2.F,data2.V] = SmoothLaplacian(F,V,15); %laplacian smoothing for surface (Kevin Moerman)
clear F, clear V

load(strcat(MrstructPath,'\',FILENAME2))
velocity = mrStruct.dataAy; clear mrstruct

for t = 1:size(velocity,5)
    vx = squeeze(velocity(:,:,:,1,t));
    vy = squeeze(velocity(:,:,:,2,t));
    vz = squeeze(velocity(:,:,:,3,t));
    vmagn = sqrt(vx.^2 + vy.^2 + vz.^2);
    mean_velo(t) = mean(vmagn(L2));
end

if plotFlag == 1
    h_meanVel=figure('Name','Patient mean velocity waveform');
    plot(1:size(velocity,5),mean_velo,'-ko','LineWidth',4,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',14);
    ylabel('Mean velocity (m/s)');
    xlabel('Time frame #');
end

[I,time] = find(mean_velo==max(mean_velo));

%%% Wall shear stress coordinates
data2.x_coor_wss = data2.V(:,1);
data2.y_coor_wss = data2.V(:,2);
data2.z_coor_wss = data2.V(:,3);

load(strcat(MrstructPath,'\',FILENAME3))
WSS = Wss_point_cloud; clear Wss_point_cloud

%%% What follows is a horrible piece of code, so if you're reading this and feel like cleaning it up, please do,
%%% I'll buy you a beer next time we meet. PvO
if peak_systolicFlag == 1
    data2.x_value_vel = velocity(:,:,:,1,time);%medfilt3(velocity(:,:,:,1,time));
    data2.y_value_vel = velocity(:,:,:,2,time);%medfilt3(velocity(:,:,:,2,time));
    data2.z_value_vel = velocity(:,:,:,3,time);%medfilt3(velocity(:,:,:,3,time));
    data2.x_value_vel = data2.x_value_vel(L2);
    data2.y_value_vel = data2.y_value_vel(L2);
    data2.z_value_vel = data2.z_value_vel(L2);
    if plotFlag == 1
        figure(h_meanVel)
        hold on, plot(time,mean_velo(time),'-ko','LineWidth',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    end
    if size(WSS,2) > 5
        data2.x_value_wss = WSS{time}(:,1);
        data2.y_value_wss = WSS{time}(:,2);
        data2.z_value_wss = WSS{time}(:,3);
    elseif size(WSS,2) == 5
        data2.x_value_wss = WSS{3}(:,1);
        data2.y_value_wss = WSS{3}(:,2);
        data2.z_value_wss = WSS{3}(:,3);
    elseif size(WSS,2) == 4
        data2.x_value_wss = WSS{2}(:,1);
        data2.y_value_wss = WSS{2}(:,2);
        data2.z_value_wss = WSS{2}(:,3);
    elseif size(WSS,2) == 3
        data2.x_value_wss = WSS{1}(:,1);
        data2.y_value_wss = WSS{1}(:,2);
        data2.z_value_wss = WSS{1}(:,3);
    elseif size(WSS,2) == 1    % Emilie: calculated at only one time (peak systole)
        data2.x_value_wss = WSS{1}(:,1);
        data2.y_value_wss = WSS{1}(:,2);
        data2.z_value_wss = WSS{1}(:,3);
    end
elseif peak_systolicFlag == 0
    % Emilie: if WSS was previously calculated only at peak systole
    if size(WSS,2) == 1
        warndlg('WSS was previously calculated only at peak systole! Please check the ''Compute at peak systole'' box.');
        return;
    % Velocity averaged over 5 systolic time frames
    elseif time == 2    % mistriggering: second time frame is peak systole, averaging over 5 timesteps is not possible
        if plotFlag == 1
            figure(h_meanVel)
            hold on, plot(time-1:time+2,mean_velo(time-1:time+2),'-ko','LineWidth',4,...
                'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
        end
        disp('TIME FRAMES AVERAGED OVER 4 TIME FRAMES!')
        data2.x_value_vel_t1 = velocity(:,:,:,1,time-1);data2.y_value_vel_t1 = velocity(:,:,:,2,time-1);data2.z_value_vel_t1 = velocity(:,:,:,3,time-1);
        data2.x_value_vel_t2 = velocity(:,:,:,1,time);  data2.y_value_vel_t2 = velocity(:,:,:,2,time);  data2.z_value_vel_t2 = velocity(:,:,:,3,time);
        data2.x_value_vel_t3 = velocity(:,:,:,1,time+1);data2.y_value_vel_t3 = velocity(:,:,:,2,time+1);data2.z_value_vel_t3 = velocity(:,:,:,3,time+1);
        data2.x_value_vel_t4 = velocity(:,:,:,1,time+2);data2.y_value_vel_t4 = velocity(:,:,:,2,time+2);data2.z_value_vel_t4 = velocity(:,:,:,3,time+2);
        data2.x_value_vel = (data2.x_value_vel_t1(L2) + data2.x_value_vel_t2(L2) + data2.x_value_vel_t3(L2) + data2.x_value_vel_t4(L2))./4;
        data2.y_value_vel = (data2.y_value_vel_t1(L2) + data2.y_value_vel_t2(L2) + data2.y_value_vel_t3(L2) + data2.y_value_vel_t4(L2))./4;
        data2.z_value_vel = (data2.z_value_vel_t1(L2) + data2.z_value_vel_t2(L2) + data2.z_value_vel_t3(L2) + data2.z_value_vel_t4(L2))./4;
        if size(WSS,2) > 5
            data2.x_value_wss_t1 = WSS{time-1}(:,1);data2.y_value_wss_t1 = WSS{time-1}(:,2);data2.z_value_wss_t1 = WSS{time-1}(:,3);
            data2.x_value_wss_t2 = WSS{time}(:,1);  data2.y_value_wss_t2 = WSS{time}(:,2);  data2.z_value_wss_t2 = WSS{time}(:,3);
            data2.x_value_wss_t3 = WSS{time+1}(:,1);data2.y_value_wss_t3 = WSS{time+1}(:,2);data2.z_value_wss_t3 = WSS{time+1}(:,3);
            data2.x_value_wss_t4 = WSS{time+2}(:,1);data2.y_value_wss_t4 = WSS{time+2}(:,2);data2.z_value_wss_t4 = WSS{time+2}(:,3);
            data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4)./4;
            data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4)./4;
            data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4)./4;
        elseif size(WSS,2) == 4
            data2.x_value_wss_t1 = WSS{time-1}(:,1);data2.y_value_wss_t1 = WSS{time-1}(:,2);data2.z_value_wss_t1 = WSS{time-1}(:,3);
            data2.x_value_wss_t2 = WSS{time}(:,1);  data2.y_value_wss_t2 = WSS{time}(:,2);  data2.z_value_wss_t2 = WSS{time}(:,3);
            data2.x_value_wss_t3 = WSS{time+1}(:,1);data2.y_value_wss_t3 = WSS{time+1}(:,2);data2.z_value_wss_t3 = WSS{time+1}(:,3);
            data2.x_value_wss_t4 = WSS{time+2}(:,1);data2.y_value_wss_t4 = WSS{time+2}(:,2);data2.z_value_wss_t4 = WSS{time+2}(:,3);
            data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4)./4;
            data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4)./4;
            data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4)./4;
        end
    elseif time == 1 % mistriggering: first time frame is peak systole, averaging over 5 timesteps is not possible
        if plotFlag == 1
            figure(h_meanVel)
            hold on, plot(time:time+2,mean_velo(time:time+2),'-ko','LineWidth',4,...
                'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
        end
        disp('TIME FRAMES AVERAGED OVER 3 TIME FRAMES!')
        data2.x_value_vel_t1 = velocity(:,:,:,1,time);  data2.y_value_vel_t1 = velocity(:,:,:,2,time);  data2.z_value_vel_t1 = velocity(:,:,:,3,time);
        data2.x_value_vel_t2 = velocity(:,:,:,1,time+1);data2.y_value_vel_t2 = velocity(:,:,:,2,time+1);data2.z_value_vel_t2 = velocity(:,:,:,3,time+1);
        data2.x_value_vel_t3 = velocity(:,:,:,1,time+2);data2.y_value_vel_t3 = velocity(:,:,:,2,time+2);data2.z_value_vel_t3 = velocity(:,:,:,3,time+2);
        data2.x_value_vel = (data2.x_value_vel_t1(L2) + data2.x_value_vel_t2(L2) + data2.x_value_vel_t3(L2))./3;
        data2.y_value_vel = (data2.y_value_vel_t1(L2) + data2.y_value_vel_t2(L2) + data2.y_value_vel_t3(L2))./3;
        data2.z_value_vel = (data2.z_value_vel_t1(L2) + data2.z_value_vel_t2(L2) + data2.z_value_vel_t3(L2))./3;
        if size(WSS,2) > 5
            data2.x_value_wss_t1 = WSS{time}(:,1);  data2.y_value_wss_t1 = WSS{time}(:,2);  data2.z_value_wss_t1 = WSS{time}(:,3);
            data2.x_value_wss_t2 = WSS{time+1}(:,1);data2.y_value_wss_t2 = WSS{time+1}(:,2);data2.z_value_wss_t2 = WSS{time+1}(:,3);
            data2.x_value_wss_t3 = WSS{time+2}(:,1);data2.y_value_wss_t3 = WSS{time+2}(:,2);data2.z_value_wss_t3 = WSS{time+2}(:,3);
            data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3)./3;
            data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3)./3;
            data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3)./3;
        elseif size(WSS,2) == 3
            data2.x_value_wss_t1 = WSS{time}(:,1);  data2.y_value_wss_t1 = WSS{time}(:,2);  data2.z_value_wss_t1 = WSS{time}(:,3);
            data2.x_value_wss_t2 = WSS{time+1}(:,1);data2.y_value_wss_t2 = WSS{time+1}(:,2);data2.z_value_wss_t2 = WSS{time+1}(:,3);
            data2.x_value_wss_t3 = WSS{time+2}(:,1);data2.y_value_wss_t3 = WSS{time+2}(:,2);data2.z_value_wss_t3 = WSS{time+2}(:,3);
            data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3)./3;
            data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3)./3;
            data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3)./3;
        end
    else % normal triggering timestep > 2 is peak systole
        if plotFlag == 1
            figure(h_meanVel)
            hold on, plot(time-2:time+2,mean_velo(time-2:time+2),'-ko','LineWidth',4,...
                'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
        end
        data2.x_value_vel_t1 = velocity(:,:,:,1,time-2);data2.y_value_vel_t1 = velocity(:,:,:,2,time-2);data2.z_value_vel_t1 = velocity(:,:,:,3,time-2);
        data2.x_value_vel_t2 = velocity(:,:,:,1,time-1);data2.y_value_vel_t2 = velocity(:,:,:,2,time-1);data2.z_value_vel_t2 = velocity(:,:,:,3,time-1);
        data2.x_value_vel_t3 = velocity(:,:,:,1,time);  data2.y_value_vel_t3 = velocity(:,:,:,2,time);  data2.z_value_vel_t3 = velocity(:,:,:,3,time);
        data2.x_value_vel_t4 = velocity(:,:,:,1,time+1);data2.y_value_vel_t4 = velocity(:,:,:,2,time+1);data2.z_value_vel_t4 = velocity(:,:,:,3,time+1);
        data2.x_value_vel_t5 = velocity(:,:,:,1,time+2);data2.y_value_vel_t5 = velocity(:,:,:,2,time+2);data2.z_value_vel_t5 = velocity(:,:,:,3,time+2);
        data2.x_value_vel = (data2.x_value_vel_t1(L2) + data2.x_value_vel_t2(L2) + data2.x_value_vel_t3(L2) + data2.x_value_vel_t4(L2) + data2.x_value_vel_t5(L2))./5;
        data2.y_value_vel = (data2.y_value_vel_t1(L2) + data2.y_value_vel_t2(L2) + data2.y_value_vel_t3(L2) + data2.y_value_vel_t4(L2) + data2.y_value_vel_t5(L2))./5;
        data2.z_value_vel = (data2.z_value_vel_t1(L2) + data2.z_value_vel_t2(L2) + data2.z_value_vel_t3(L2) + data2.z_value_vel_t4(L2) + data2.z_value_vel_t5(L2))./5;
        if size(WSS,2) > 5
            data2.x_value_wss_t1 = WSS{time-2}(:,1);data2.y_value_wss_t1 = WSS{time-2}(:,2);data2.z_value_wss_t1 = WSS{time-2}(:,3);
            data2.x_value_wss_t2 = WSS{time-1}(:,1);data2.y_value_wss_t2 = WSS{time-1}(:,2);data2.z_value_wss_t2 = WSS{time-1}(:,3);
            data2.x_value_wss_t3 = WSS{time}(:,1);  data2.y_value_wss_t3 = WSS{time}(:,2);  data2.z_value_wss_t3 = WSS{time}(:,3);
            data2.x_value_wss_t4 = WSS{time+1}(:,1);data2.y_value_wss_t4 = WSS{time+1}(:,2);data2.z_value_wss_t4 = WSS{time+1}(:,3);
            data2.x_value_wss_t5 = WSS{time+2}(:,1);data2.y_value_wss_t5 = WSS{time+2}(:,2);data2.z_value_wss_t5 = WSS{time+2}(:,3);
            data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4 + data2.x_value_wss_t5)./5;
            data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4 + data2.y_value_wss_t5)./5;
            data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4 + data2.z_value_wss_t5)./5;
        elseif size(WSS,2) == 5
            time = 3;
            data2.x_value_wss_t1 = WSS{time-2}(:,1);data2.y_value_wss_t1 = WSS{time-2}(:,2);data2.z_value_wss_t1 = WSS{time-2}(:,3);
            data2.x_value_wss_t2 = WSS{time-1}(:,1);data2.y_value_wss_t2 = WSS{time-1}(:,2);data2.z_value_wss_t2 = WSS{time-1}(:,3);
            data2.x_value_wss_t3 = WSS{time}(:,1);  data2.y_value_wss_t3 = WSS{time}(:,2);  data2.z_value_wss_t3 = WSS{time}(:,3);
            data2.x_value_wss_t4 = WSS{time+1}(:,1);data2.y_value_wss_t4 = WSS{time+1}(:,2);data2.z_value_wss_t4 = WSS{time+1}(:,3);
            data2.x_value_wss_t5 = WSS{time+2}(:,1);data2.y_value_wss_t5 = WSS{time+2}(:,2);data2.z_value_wss_t5 = WSS{time+2}(:,3);
            data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4 + data2.x_value_wss_t5)./5;
            data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4 + data2.y_value_wss_t5)./5;
            data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4 + data2.z_value_wss_t5)./5;
        end
    end
end

% Velocity and WSS magnitude (scalar)
data2.vel_m = sqrt(data2.x_value_vel.^2 + data2.y_value_vel.^2 + data2.z_value_vel.^2);
data2.wss_m = sqrt(data2.x_value_wss.^2 + data2.y_value_wss.^2 + data2.z_value_wss.^2);

if plotFlag == 1
    
    figure('Name', 'Patient velocity and WSS 3D vectors')
    subplot(1,2,1)
    patch('Faces',data2.F,'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss], ...
        'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.25);
    hold on
    c = [];
    a = [2 20];
    [F,V,C]=quiver3Dpatch(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,data2.y_value_vel,data2.x_value_vel,data2.z_value_vel,c,a);
    patch('Faces',F,'Vertices',V,'CData',C,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.75);
    caxis([0 1.5]);axis equal;view([-180 90]);axis off
    title('Velocity (m/s)')
    h_cb = colorbar;
    pos_cb = get(h_cb,'Position');
    set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])
    
    subplot(1,2,2)
    a = [2 15];
    c = [ ];
    patch('Faces',data2.F,'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss], ...
        'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
    hold on
    [F2,V2,C2]=quiver3Dpatch(data2.x_coor_wss,data2.y_coor_wss,data2.z_coor_wss,data2.x_value_wss ...
        ,data2.y_value_wss,data2.z_value_wss,c,a);
    patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
    caxis([0 1.5])
    axis equal;axis off; axis ij
    view([-180 -90])  
    title('WSS (Pa)')
    h_cb = colorbar;
    pos_cb = get(h_cb,'Position');
    set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])
    
    figure('Name', 'Registration between atlas and patient')
    subplot(2,3,1)
    plot3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,'r.')
    hold on
    plot3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,'b.')
    axis equal; axis ij; axis off; view([-180 -90])
    subplot(2,3,2)
    plot3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,'r.')
    hold on
    plot3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,'b.')
    axis equal; axis ij; axis off; view([-180 -90]); camorbit(-90,0,'data',[0 1 0])
    legend('to remain the same','to be transformed')
    title('before registration')
    subplot(2,3,3)
    plot3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,'r.')
    hold on
    plot3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,'b.')
    axis equal; axis ij; axis off; view([-180 0])
    
%     subplot(1,2,1)
%     L_figure = (squeeze(max(atlas_matrix,[],3))~=0);
%     imagesc(squeeze(max(atlas_matrix,[],3)),'Alphadata',double(L_figure));
%     axis tight; axis equal; axis ij; axis off;caxis([0 1.2]);
%     title('data2 velocity')
%     h_cb = colorbar;
%     pos_cb = get(h_cb,'Position');
%     set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])
%     
%     subplot(1,2,2)
%     patch('Faces',data2.F,'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss], ...
%         'EdgeColor','none','FaceVertexCData',data2.wss_m,'FaceColor','interp','FaceAlpha',1);
%     caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
%     title('data2 WSS')
%     h_cb = colorbar;
%     pos_cb = get(h_cb,'Position');
%     set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])
    
end

%%% Velocities
PSF = fspecial('gaussian',5,1);
mask1_to_register = mask1;
mask1_to_register = imfilter(mask1_to_register,PSF,'conv');

mask2 = imfilter(mask2,PSF,'conv');

disp(' ')
disp('affine registration (dof = 12)! So only scalars used (Affine registration doesnt work for vectors)...Busy registering...')
disp(' ')
tic
disp('...This can take up to 5 minutes...')

% directory with flirt.exe and cygwin1.dll (use cygwin convention)
fsldir = [path_flirt '/']

% save as nifti
cnii=make_nii(mask1_to_register,[atlas.vox(1) atlas.vox(2) atlas.vox(3)]);
save_nii(cnii,'mask1.nii');
mnii=make_nii(mask2,[mask2_vox(1) mask2_vox(2) mask2_vox(3)]);
save_nii(mnii,'mask2.nii');

% run flirt (needs cygwin installed)
infile = 'mask2.nii';
reffile = 'mask1.nii';
outfile = 'rmask2.nii';
omat = 'Rotation_Translation.mat';

flirtcmd= [...
    'flirt -searchrx -0 0 -searchry -0 0 -searchrz -0 0 -dof 12' ...
    '-in ' infile ...
    ' -ref ' reffile ' -out ' outfile ' -omat ' omat];
flirtcmd=[fsldir flirtcmd];
flirtcmd = ['export FSLOUTPUTTYPE=NIFTI;' flirtcmd];

% create file with run command
f=fopen('runflirt.sh','wt');
fprintf(f,'%s',flirtcmd);
fclose(f);

% and go! takes 4-5 mins.
%system('c:\cygwin64\bin\bash runflirt.sh');
system([path_cygwin '\bash runflirt.sh']);

% load transformation mask
load Rotation_Translation -ascii

flirtmat = 'Rotation_Translation.mat';
src = infile;
trg = reffile;
[worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
rotmat = worldmat(1:3,1:3);

%%% Velocity coordinates
yxz_coor_vel = [data2.y_coor_vel data2.x_coor_vel data2.z_coor_vel];
yxz_coor_vel(:,4) = 1;
yxz_coor_vel_new = inv(worldmat)*yxz_coor_vel'; clear yxz_coor_vel
data2.x_coor_vel_new = yxz_coor_vel_new(2,:)';
data2.y_coor_vel_new = yxz_coor_vel_new(1,:)';
data2.z_coor_vel_new = yxz_coor_vel_new(3,:)'; clear yxz_coor_vel_new

%%% WSS coordinates
yxz_coor_wss = [data2.y_coor_wss data2.x_coor_wss data2.z_coor_wss];
yxz_coor_wss(:,4) = 1;
yxz_coor_wss_new = inv(worldmat)*yxz_coor_wss'; clear yxz_coor_wss
data2.x_coor_wss_new = yxz_coor_wss_new(2,:)';
data2.y_coor_wss_new = yxz_coor_wss_new(1,:)';
data2.z_coor_wss_new = yxz_coor_wss_new(3,:)'; clear yxz_coor_wss_new
toc

if plotFlag == 1
    subplot(2,3,4)
    plot3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,'r.')
    hold on
    plot3(data2.x_coor_vel_new,data2.y_coor_vel_new,data2.z_coor_vel_new,'b.')
    axis equal; axis ij; axis off; view([-180 -90])
    subplot(2,3,5)
    plot3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,'r.')
    hold on
    plot3(data2.x_coor_vel_new,data2.y_coor_vel_new,data2.z_coor_vel_new,'b.')
    axis equal; axis ij; axis off; view([-180 -90]); camorbit(-90,0,'data',[0 1 0])
    title('registered')
    subplot(2,3,6)
    plot3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,'r.')
    hold on
    plot3(data2.x_coor_vel_new,data2.y_coor_vel_new,data2.z_coor_vel_new,'b.')
    axis equal; axis ij; axis off; view([-180 0])
end

choice = questdlg('Is registration to the atlas geometry correct?', ...
    'Heatmap creation... registration step', ...
    'Yes','No','Yes');
% Handle response
switch choice
    case 'Yes'

    case 'No'
        return;
end

offset = 100;
x_vel_round = round(data2.x_coor_vel_new./atlas.vox(1)) + offset;
y_vel_round = round(data2.y_coor_vel_new./atlas.vox(2)) + offset;
z_vel_round = round(data2.z_coor_vel_new./atlas.vox(3)) + offset;

indices_mask2 = [x_vel_round y_vel_round z_vel_round];
siz=max(indices_mask2,[],1);
IND = sub2ind(siz,x_vel_round,y_vel_round,z_vel_round);
[b, IND_double_removed, n] = unique(IND);
clear b, clear n
indices_mask2 = [x_vel_round(IND_double_removed) y_vel_round(IND_double_removed) z_vel_round(IND_double_removed)];
clear IND_double_removed, clear IND

mask_new = zeros([max(y_vel_round) max(x_vel_round) max(z_vel_round)]);

for i = 1:size(indices_mask2,1)
    mask_new(indices_mask2(i,2),indices_mask2(i,1),indices_mask2(i,3)) = 1;
end

% Due to the rounding of coordinates there are holes in the aorta, we fill them by smooth3
% and then erode the aorta to give it the right size again
se = strel('disk',1);
mask_new = imerode(smooth3(mask_new),se);

%%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
sizes = [size(mask1,1)+offset size(mask1,2)+offset size(mask1,3)+offset];
mask1b = zeros(sizes);
mask1b((offset+1):size(mask1b,1),(offset+1):size(mask1b,2),(offset+1):size(mask1b,3)) = mask1;
mask1 = mask1b;clear mask2b

% Make sure both masks have the same dimensions
if size(mask1,1) > size(mask_new,1)
    mask_new(size(mask_new,1):size(mask1,1),:,:) = 0;
elseif size(mask1,1) < size(mask_new,1)
    mask_new(size(mask1,1)+1:size(mask_new,1),:,:) = [];
end
if size(mask1,2) > size(mask_new,2)
    mask_new(:,size(mask_new,2):size(mask1,2),:) = 0;
elseif size(mask1,2) < size(mask_new,2)
    mask_new(:,size(mask1,2)+1:size(mask_new,2),:) = [];
end
if size(mask1,3) > size(mask_new,3)
    mask_new(:,:,size(mask_new,3):size(mask1,3)) = 0;
elseif size(mask1,3) < size(mask_new,3)
    mask_new(:,:,size(mask1,3)+1:size(mask_new,3)) = [];
end

L_mask2 = double(mask_new ~= 0);

difference = abs(mask1-L_mask2);
[I1,J] = find(mask1~=0);
[I2,J] = find(L_mask2~=0);
mean_I = (size(I1,1) + size(I2,1))/2;
[I_diff,J] = find(difference~=0);
diff_voxels = size(I_diff,1);
diff_percentage = ((diff_voxels / mean_I) * 100)/2;
disp(['RE: Difference between aorta and atlas = ' num2str(diff_percentage)])
disp(' ')

% Interpolate Velocity
interpolation_function = TriScatteredInterp([atlas.x_coor_vel atlas.y_coor_vel atlas.z_coor_vel],atlas.mean_vel,'nearest');
atlas_mean2 = interpolation_function([data2.x_coor_vel_new data2.y_coor_vel_new data2.z_coor_vel_new]);
atlas_mean2(isnan(atlas_mean2)) = 0;
atlas_mean_vel = atlas_mean2;

interpolation_function = TriScatteredInterp([atlas.x_coor_vel atlas.y_coor_vel atlas.z_coor_vel],atlas.std_vel,'nearest');
atlas_SD2 = interpolation_function([data2.x_coor_vel_new data2.y_coor_vel_new data2.z_coor_vel_new]);
atlas_SD2(isnan(atlas_SD2)) = 0;
atlas_std_vel = atlas_SD2;

% Interpolate WSS
interpolation_function = TriScatteredInterp([atlas.x_coor_wss atlas.y_coor_wss atlas.z_coor_wss],atlas.mean_wss,'nearest');
atlas_mean2 = interpolation_function([data2.x_coor_wss_new data2.y_coor_wss_new data2.z_coor_wss_new]);
atlas_mean2(isnan(atlas_mean2)) = 0;
atlas_mean_wss = atlas_mean2;

interpolation_function = TriScatteredInterp([atlas.x_coor_wss atlas.y_coor_wss atlas.z_coor_wss],atlas.std_wss,'nearest');
atlas_SD2 = interpolation_function([data2.x_coor_wss_new data2.y_coor_wss_new data2.z_coor_wss_new]);
atlas_SD2(isnan(atlas_SD2)) = 0;
atlas_std_wss = atlas_SD2;

if calculateIE_Flag == 1;
    
    if ~exist(strcat(PATHNAME,'\atlas_interpolation_error_ROI_after_transformation\mask1.mat'),'file')
        
        F2=figure('Name','Atlas shape: Paused after finishing a region so press space when finished!');
        plot3(data2.x_coor_wss_new,data2.y_coor_wss_new,data2.z_coor_wss_new,'r.');
        view([-180 -90]);axis ij;axis equal;axis off
        
        mkdir(PATHNAME,'\atlas_interpolation_error_ROI_after_transformation')
        
        for i = 1:3
            %Polygon and mask for AAo
            polyAAo = impoly;
            wait(polyAAo);
            region = getPosition(polyAAo);
            
            %          disp('saving, pausing')
            save(strcat([PATHNAME '\atlas_interpolation_error_ROI_after_transformation\mask' num2str(i)]),'region');
            pause
        end
        
        close(F2)
    end
    load(strcat(PATHNAME,'\atlas_interpolation_error_ROI_after_transformation\mask1'))
    transformed_atlas_mask_AAo_inner_vel = inpolygon(data2.x_coor_vel_new,data2.y_coor_vel_new, region(:,1), region(:,2));
    transformed_mean_vel_asc_inner = mean(atlas_mean_vel(transformed_atlas_mask_AAo_inner_vel));
    transformed_atlas_mask_AAo_inner_wss = inpolygon(data2.x_coor_wss_new,data2.y_coor_wss_new, region(:,1), region(:,2));
    transformed_mean_wss_asc_inner = mean(atlas_mean_wss(transformed_atlas_mask_AAo_inner_wss));
    load(strcat(PATHNAME,'\atlas_interpolation_error_ROI_after_transformation\mask2'))
    transformed_atlas_mask_AAo_outer_vel = inpolygon(data2.x_coor_vel_new,data2.y_coor_vel_new, region(:,1), region(:,2));
    transformed_mean_vel_asc_outer = mean(atlas_mean_vel(transformed_atlas_mask_AAo_outer_vel));
    transformed_atlas_mask_AAo_outer_wss = inpolygon(data2.x_coor_wss_new,data2.y_coor_wss_new, region(:,1), region(:,2));
    transformed_mean_wss_asc_outer = mean(atlas_mean_wss(transformed_atlas_mask_AAo_outer_wss));
    load(strcat(PATHNAME,'\atlas_interpolation_error_ROI_after_transformation\mask3'))
    transformed_atlas_mask_arch_inner_vel = inpolygon(data2.x_coor_vel_new,data2.y_coor_vel_new, region(:,1), region(:,2));
    transformed_mean_vel_arch_inner = mean(atlas_mean_vel(transformed_atlas_mask_arch_inner_vel));
    transformed_atlas_mask_arch_inner_wss = inpolygon(data2.x_coor_wss_new,data2.y_coor_wss_new, region(:,1), region(:,2));
    transformed_mean_wss_arch_inner = mean(atlas_mean_wss(transformed_atlas_mask_arch_inner_wss));

    mean_vel_atlas_total_after_interpolation = mean(atlas_mean_vel)
    mean_wss_atlas_total_after_interpolation = mean(atlas_mean_wss)
    
    mean_vel_after_interpolation(1,1) = transformed_mean_vel_asc_inner;
    mean_vel_after_interpolation(2,1) = transformed_mean_vel_asc_outer;
    mean_vel_after_interpolation(3,1) = transformed_mean_vel_arch_inner

    mean_wss_after_interpolation(1,1) = transformed_mean_wss_asc_inner;
    mean_wss_after_interpolation(2,1) = transformed_mean_wss_asc_outer;
    mean_wss_after_interpolation(3,1) = transformed_mean_wss_arch_inner
    
    if plotFlag == 1
        figure('Name','Velocity after interpolation: inner AAo')
        scatter3(data2.x_coor_vel_new(transformed_atlas_mask_AAo_inner_vel),data2.y_coor_vel_new(transformed_atlas_mask_AAo_inner_vel),data2.z_coor_vel_new(transformed_atlas_mask_AAo_inner_vel),20,atlas_mean_vel(transformed_atlas_mask_AAo_inner_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','Velocity after interpolation: outer AAo')
        scatter3(data2.x_coor_vel_new(transformed_atlas_mask_AAo_outer_vel),data2.y_coor_vel_new(transformed_atlas_mask_AAo_outer_vel),data2.z_coor_vel_new(transformed_atlas_mask_AAo_outer_vel),20,atlas_mean_vel(transformed_atlas_mask_AAo_outer_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','Velocity after interpolation: inner arch')
        scatter3(data2.x_coor_vel_new(transformed_atlas_mask_arch_inner_vel),data2.y_coor_vel_new(transformed_atlas_mask_arch_inner_vel),data2.z_coor_vel_new(transformed_atlas_mask_arch_inner_vel),20,atlas_mean_vel(transformed_atlas_mask_arch_inner_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);

        figure('Name','WSS after interpolation: inner AAo')
        scatter3(data2.x_coor_wss_new(transformed_atlas_mask_AAo_inner_wss),data2.y_coor_wss_new(transformed_atlas_mask_AAo_inner_wss),atlas.z_coor_wss(transformed_atlas_mask_AAo_inner_wss),20,atlas_mean_wss(transformed_atlas_mask_AAo_inner_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS after interpolation: outer AAo')
        scatter3(data2.x_coor_wss_new(transformed_atlas_mask_AAo_outer_wss),data2.y_coor_wss_new(transformed_atlas_mask_AAo_outer_wss),data2.z_coor_wss_new(transformed_atlas_mask_AAo_outer_wss),20,atlas_mean_wss(transformed_atlas_mask_AAo_outer_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS before interpolation: inner arch')
        scatter3(data2.x_coor_wss_new(transformed_atlas_mask_arch_inner_wss),data2.y_coor_wss_new(transformed_atlas_mask_arch_inner_wss),data2.z_coor_wss_new(transformed_atlas_mask_arch_inner_wss),20,atlas_mean_wss(transformed_atlas_mask_arch_inner_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);

    end
    
    IE_inner_AAo_vel = abs(mean_vel_before_interpolation(1,1)-mean_vel_after_interpolation(1,1)) / ((mean_vel_before_interpolation(1,1)+mean_vel_after_interpolation(1,1))./2)*100;
    IE_outer_AAo_vel = abs(mean_vel_before_interpolation(2,1)-mean_vel_after_interpolation(2,1)) / ((mean_vel_before_interpolation(2,1)+mean_vel_after_interpolation(2,1))./2)*100;
    IE_inner_asc_vel = abs(mean_vel_before_interpolation(3,1)-mean_vel_after_interpolation(3,1)) / ((mean_vel_before_interpolation(3,1)+mean_vel_after_interpolation(3,1))./2)*100;
    IE_total_vel = abs(mean_vel_atlas_total_before_interpolation-mean_vel_atlas_total_after_interpolation) / ((mean_vel_atlas_total_before_interpolation+mean_vel_atlas_total_after_interpolation)./2)*100;
    IE_inner_AAo_wss = abs(mean_wss_before_interpolation(1,1)-mean_wss_after_interpolation(1,1)) / ((mean_wss_before_interpolation(1,1)+mean_wss_after_interpolation(1,1))./2)*100;
    IE_outer_AAo_wss = abs(mean_wss_before_interpolation(2,1)-mean_wss_after_interpolation(2,1)) / ((mean_wss_before_interpolation(2,1)+mean_wss_after_interpolation(2,1))./2)*100;
    IE_inner_asc_wss = abs(mean_wss_before_interpolation(3,1)-mean_wss_after_interpolation(3,1)) / ((mean_wss_before_interpolation(3,1)+mean_wss_after_interpolation(3,1))./2)*100;
    IE_total_wss = abs(mean_wss_atlas_total_before_interpolation-mean_wss_atlas_total_after_interpolation) / ((mean_wss_atlas_total_before_interpolation+mean_wss_atlas_total_after_interpolation)./2)*100;
    disp(['IE velocity inner AAo = ' num2str(IE_inner_AAo_vel) ' %'])
    disp(['IE velocity outer AAo = ' num2str(IE_outer_AAo_vel) ' %'])
    disp(['IE velocity inner asc = ' num2str(IE_inner_asc_vel) ' %'])
    disp(['IE velocity total = ' num2str(IE_total_vel) ' %'])
    disp(' ')
    disp(['IE wall shear stress inner AAo = ' num2str(IE_inner_AAo_wss) ' %'])
    disp(['IE wall shear stress outer AAo = ' num2str(IE_outer_AAo_wss) ' %'])
    disp(['IE wall shear stress inner asc = ' num2str(IE_inner_asc_wss) ' %'])
    disp(['IE wss total = ' num2str(IE_total_wss) ' %'])
    disp(' ')
end

if plotFlag == 1
        
    h=figure('Name','Transformed velocity and WSS atlas to patient geometry');
    subplot(1,2,1)
    atlas_matrix = zeros(size(L2)); 
    atlas_matrix(L2) = atlas_mean_vel;
    L_figure = (squeeze(max(atlas_matrix,[],3))~=0);
    imagesc(flipdim(squeeze(max(atlas_matrix,[],3)),2),'Alphadata',double(flipdim(L_figure,2)));
    axis tight; axis equal; axis ij; axis off;caxis([0 1.5]);
    title('Velocity (m/s)');
    h_cb = colorbar;
    pos_cb = get(h_cb,'Position');
    set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])
    subplot(1,2,2)
    patch('Faces',data2.F,'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',atlas_mean_wss,'FaceColor','interp','FaceAlpha',1);
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
    title('WSS (Pa)');
    h_cb = colorbar;
    pos_cb = get(h_cb,'Position');
    set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])    
    
%     scatter3(data2.x_coor_vel_new,data2.y_coor_vel_new,data2.z_coor_vel_new,atlas_mean_vel.*50,atlas_mean_vel,'filled')
%     axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90]);
%     title('Mean transformed velocity atlas');
%     subplot(1,3,2)
%     scatter3(data2.x_coor_vel_new,data2.y_coor_vel_new,data2.z_coor_vel_new,atlas_std_vel.*50,atlas_std_vel,'filled')
%     axis equal;axis off; axis ij;
%     caxis([0 1.5])
%     view([180 -90])
%     title('std transformed velocity atlas');
%     
%     set(h,'Position', [680   558   672   420])
%     
%     patch('Faces',data2.F,'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',atlas_std_wss,'FaceColor','interp','FaceAlpha',1);
%     axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
%     title('std transformed WSS atlas');
%     subplot(1,3,3)
%     atlas_matrix = zeros(size(L2));
%     atlas_matrix(L2) = atlas_std_vel;
%     L_figure = (squeeze(max(L2,[],3))~=0);
%     imagesc(squeeze(max(atlas_matrix,[],3)),'Alphadata',double(L_figure));
%     axis tight; axis equal; axis ij; axis off;caxis([0 1.5]);
%     title('Mean transformed WSS atlas 2');
% 
%     set(h,'Position', [680   558   672   420])
    
end

% Determine thresholds
mean_plus_2SD_atlas_vel = atlas_mean_vel + 1.96.*atlas_std_vel;
mean_plus_1SD_atlas_vel = atlas_mean_vel + 0.98.*atlas_std_vel;
mean_min_2SD_atlas_vel = atlas_mean_vel - 1.96.*atlas_std_vel;
mean_plus_2SD_atlas_wss = atlas_mean_wss + 1.96.*atlas_std_wss;
mean_min_2SD_atlas_wss = atlas_mean_wss - 1.96.*atlas_std_wss;

% confidence interval figures: create switch to keep display?
% if plotFlag == 1
%     % Velocity
% %     figure('Name','mean_plus_2SD_atlas velocity')
% %     scatter3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,atlas_mean_vel.*50,mean_plus_2SD_atlas_vel,'filled')
% %     colorbar;axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90]);
% 
% %     figure('Name','mean_plus_1SD_atlas velocity')
% %     scatter3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,atlas_mean_vel.*50,mean_plus_1SD_atlas_vel,'filled')
% %     colorbar;axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90]);
% 
%     figure('Name', 'Velocity confidence interval atlas')
%     subplot(1,2,1)
%     atlas_matrix = zeros(size(L2));
%     atlas_matrix(L2) = mean_min_2SD_atlas_vel;
%     L_figure = (squeeze(max(L2,[],3))~=0);
%     imagesc(squeeze(max(atlas_matrix,[],3)),'Alphadata',double(L_figure));
%     axis tight; axis equal; axis ij; axis off;caxis([0 1.5]);
%     title('mean minus 2SD atlas velocity')
%     subplot(1,2,2)
%     atlas_matrix = zeros(size(L2)); 
%     atlas_matrix(L2) = mean_plus_2SD_atlas_vel;
%     L_figure = (squeeze(max(atlas_matrix,[],3))~=0);
%     imagesc(squeeze(max(atlas_matrix,[],3)),'Alphadata',double(L_figure));
%     axis tight; axis equal; axis ij; axis off;caxis([0 1.5]);
%     title('mean plus 2SD atlas velocity')
%     h_cb = colorbar;
%     pos_cb = get(h_cb,'Position');
%     set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])
% 
%     % WSS
%     figure('Name', 'WSS confidence interval atlas')
%     subplot(1,2,1)
%     patch('Faces',data2.F,'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',mean_min_2SD_atlas_wss,'FaceColor','interp','FaceAlpha',1);
%     axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90]);
%     title('mean minus 2SD atlas WSS')
%     subplot(1,2,2)
%     patch('Faces',data2.F,'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',mean_plus_2SD_atlas_wss,'FaceColor','interp','FaceAlpha',1);
%     axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90]);
%     title('mean plus 2SD atlas WSS')
%     h_cb = colorbar;
%     pos_cb = get(h_cb,'Position');
%     set(h_cb,'Position',[pos_cb(1)+.1 pos_cb(2) pos_cb(3) pos_cb(4)])
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TRAFFIC LIGHT MAP FOR ABNORMAL VELOCITY %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_mask_red_ = zeros(size(data2.vel_m));
new_mask_yellow_ = zeros(size(data2.vel_m));
new_mask_blue_ = zeros(size(data2.vel_m));
new_mask_green_ = zeros(size(data2.vel_m));
for i=1:size(data2.vel_m,1)
    if data2.vel_m(i) > mean_plus_2SD_atlas_vel(i)
        new_mask_red_(i,1) = 1;
        new_mask_yellow_(i,1) = -1;
        new_mask_blue_(i,1) = -1;        
        new_mask_green_(i,1) = -1;
    elseif data2.vel_m(i) > mean_plus_1SD_atlas_vel(i)
        new_mask_red_(i,1) = -1;
        new_mask_yellow_(i,1) = 1;
        new_mask_blue_(i,1) = -1;         
        new_mask_green_(i,1) = -1;
    elseif data2.vel_m(i) < mean_min_2SD_atlas_vel(i)
        new_mask_red_(i,1) = -1;
        new_mask_yellow_(i,1) = -1;
        new_mask_blue_(i,1) = 1;         
        new_mask_green_(i,1) = -1;        
    else
        new_mask_red_(i,1) = -1;
        new_mask_yellow_(i,1) = -1;
        new_mask_blue_(i,1) = -1;          
        new_mask_green_(i,1) = 1;
    end
end

% translate point cloud into matrix
new_mask_red = zeros(size(L2));
new_mask_red(L2) = new_mask_red_;
new_mask_red(~L2) = -1;
new_mask_yellow = zeros(size(L2));
new_mask_yellow(L2) = new_mask_yellow_;
new_mask_yellow(~L2) = -1;
new_mask_blue = zeros(size(L2));
new_mask_blue(L2) = new_mask_blue_;
new_mask_blue(~L2) = -1;
new_mask_green = zeros(size(L2));
new_mask_green(L2) = new_mask_green_;
new_mask_green(~L2) = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% HEAT MAP FOR ABNORMAL WSS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heat_mapp = zeros(size(data2.wss_m,1),1);
for i=1:size(data2.wss_m,1)
   if data2.wss_m(i) < mean_min_2SD_atlas_wss(i)
       heat_mapp(i,1) = 0;
   elseif data2.wss_m(i) > mean_plus_2SD_atlas_wss(i)
        heat_mapp(i,1) = 1;
    else
        heat_mapp(i,1) = 2;
    end
end

color1 = ones([63 3]);

% 0 = blue
color1(1:21,1) = color1(1:21,1).*0;%color(1,1);
color1(1:21,2) = color1(1:21,2).*0;%color(1,2);
color1(1:21,3) = color1(1:21,3).*1;%color(1,3);

% 1 = red
color1(22:42,1) = color1(22:42,1).*1;%color(17,1);
color1(22:42,2) = color1(22:42,2).*0;%color(17,2);
color1(22:42,3) = color1(22:42,3).*0;%color(17,3);

% 2 = gray
color1(43:63,1) = color1(43:63,1).*0.5;%0;%color(33,1);
color1(43:63,2) = color1(43:63,2).*0.5;%1;%color(33,2);
color1(43:63,3) = color1(43:63,3).*0.5;%0;%color(33,3);

[I,J] = find(L2~=0);
total_volume = mask2_vox(1)*mask2_vox(2)*mask2_vox(3)*size(I,1);

[I,J] = find(smooth3(new_mask_red)>0);
red_volume = mask2_vox(1)*mask2_vox(2)*mask2_vox(3)*size(I,1);
percentage_red_volume = red_volume / total_volume * 100;
[I,J] = find(smooth3(new_mask_yellow)>0);
yellow_volume = mask2_vox(1)*mask2_vox(2)*mask2_vox(3)*size(I,1);
percentage_yellow_volume = yellow_volume / total_volume * 100;
[I,J] = find(smooth3(new_mask_blue)>0);
blue_volume = mask2_vox(1)*mask2_vox(2)*mask2_vox(3)*size(I,1);
percentage_blue_volume = blue_volume / total_volume * 100;
[I,J] = find(smooth3(new_mask_green)>0);
green_volume = mask2_vox(1)*mask2_vox(2)*mask2_vox(3)*size(I,1);
percentage_green_volume = green_volume / total_volume * 100;
total_percentage = percentage_red_volume + percentage_yellow_volume + percentage_blue_volume + percentage_green_volume;

disp(['Red volume percentage of total aorta = ' num2str(round(percentage_red_volume)) ' % (' num2str(round(red_volume./1000)) ' cm3)'])
disp(['Yellow volume percentage of total aorta = ' num2str(round(percentage_yellow_volume)) ' % (' num2str(round(yellow_volume./1000)) ' cm3)'])
disp(['Blue volume percentage of total aorta = ' num2str(round(percentage_blue_volume)) ' % (' num2str(round(blue_volume./1000)) ' cm3)'])
disp(['Green volume percentage of total aorta = ' num2str(round(percentage_green_volume)) ' % (' num2str(round(green_volume./1000)) ' cm3)'])
disp(['Total percentage of total aorta = ' num2str(total_percentage) ' % (' num2str(round(total_volume./1000)) ' cm3)'])
disp(' ')

[I1,J1] = find(heat_mapp == 0);
[I2,J2] = find(heat_mapp == 1);
percentage_significant_higher_than_controls = size(I2,1) / size(heat_mapp,1) * 100;
percentage_significant_lower_than_controls = size(I1,1) / size(heat_mapp,1) * 100;

disp(['Percentage higher than controls inner AAo = ' num2str(round(percentage_significant_higher_than_controls)) '%'])
disp(['Percentage lower than controls inner AAo = ' num2str(round(percentage_significant_lower_than_controls)) '%'])
disp(' ')

color2(1,:) = [0 0 1];
color2(2,:) = [1 0 0];
color2(3,:) = [0.5 0.5 0.5];
color2(4:64,:) = gray_colormap(4:64,:);

if calculate_area_of_higherlowerFlag == 1;
    
%     vertices = [data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss];
    x1 = data2.x_coor_wss;%/mask2_vox(1);
    y1 = data2.y_coor_wss;%/mask2_vox(2);
    z1 = data2.z_coor_wss;%/mask2_vox(3);
    vertices = [x1 y1 z1];
    
    f5=figure('Name', 'Calculation of higher and lower WSS areas')
    patch('Faces',data2.F,'Vertices',vertices, ...
        'EdgeColor','none','FaceVertexCData',heat_mapp,'FaceColor','interp','FaceAlpha',1);
    axis equal;axis off; axis ij
    view([-180 -90])
    
%     gray_colormap = colormap(gray);
%     color3(1,:) = [0 0 1];
%     color3(2,:) = [1 0 0];
%     color3(3,:) = [0.5 0.5 0.5];
%     color3(4:64,:) = gray_colormap(4:64,:);
    colormap(color2);
    caxis([0 64]);
%     load(strcat(MrstructPath,'\',FILENAME4))
%     magnitude = flipdim(double(mrStruct.dataAy(:,:,:,3)),3);
%     magnitude(magnitude == 0) = 3;
%     magnitude(magnitude == 1) = 3;
%     magnitude(magnitude == 2) = 3;
%     hold on,
%     s4 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3),magnitude(:,:,1),'EdgeColor','none');
    title({'Please keep in mind that regions of interest will be numbered in the same order you draw them';'once you''re done drawing an ROI, please double-click to validate and press space to move on to the next one'})
    
%     uicontrol('Style','text',...
%         'Position',[10 200 120 70],...
%         'String','Please choose using the slider the magnitude slice')
%     uicontrol('Style','text',...
%         'Position',[10 75 120 20],...
%         'String','Slice slider')
%     sl1 = uicontrol('Style', 'slider',...
%         'Min',1,'Max',size(magnitude,3),'Value',size(magnitude,3)/2,...
%         'Position', [10 50 120 20],...
%         'SliderStep',[1/(size(magnitude,3)-1) 10/(size(magnitude,3)-1)],...
%         'Callback', {@move_slice4,gca});
    
    mkdir(strcat(PATHNAME,'\heat_map_higher_lower_masks'));
    
    choice = questdlg('Do you want to draw new ROIs or load previous ones?', ...
        'ROIs creation', ...
        'Draw new','Load existing','Draw new');
    % Handle response
    switch choice
        case 'Draw new'
            
            nbROIs = inputdlg('How many ROIs do you need for calculation of higher and lower WSS areas?','Number of ROIs',1,{'3'});
            nbROIs = (round(str2num(nbROIs{1})));
            
            figure(f5)
            for i = 1:nbROIs
                %Polygon and mask
                polyAAo = impoly;
                wait(polyAAo);
                region = getPosition(polyAAo);
                disp('saving, pausing')
                save(strcat([PATHNAME '\heat_map_higher_lower_masks\mask' num2str(i)]),'region');
                text(sum(region(:,1))/size(region,1),sum(region(:,2))/size(region,1),strcat('ROI',num2str(i)),'fontweight','bold')
                clear region
                pause
            end

            h1 = waitbar(0,'ROIs drawn, calculation of higher and lower WSS areas in progress...');
            
            % save in Excel
            indices{1,1} = 'ROI #';
            indices{2,1} = 'total area (cm2)';
            indices{3,1} = 'red area (cm2)';
            indices{4,1} = 'blue area (cm2)';
            
            for i=1:nbROIs
                
                load(strcat(PATHNAME, '\heat_map_higher_lower_masks\mask', num2str(i)));
                mask = inpolygon(x1, y1, region(:,1), region(:,2));
                new_heatmap = heat_mapp.*mask;
                new_heatmap(~mask)=3;
                new_sel_blue = (new_heatmap == 0);
                new_sel_red = (new_heatmap == 1);
                new_sel_gray = (new_heatmap == 2);
                new_sel_blue = sum([new_sel_blue(data2.F(:,1)) new_sel_blue(data2.F(:,2)) new_sel_blue(data2.F(:,3))],2)>1;
                new_sel_red = sum([new_sel_red(data2.F(:,1)) new_sel_red(data2.F(:,2)) new_sel_red(data2.F(:,3))],2)>1;
                new_sel_gray = sum([new_sel_gray(data2.F(:,1)) new_sel_gray(data2.F(:,2)) new_sel_gray(data2.F(:,3))],2)>1;
                new_area_blue = triangleArea3d(vertices(data2.F(new_sel_blue,1),:),vertices(data2.F(new_sel_blue,2),:),vertices(data2.F(new_sel_blue,3),:));
                new_area_red = triangleArea3d(vertices(data2.F(new_sel_red,1),:),vertices(data2.F(new_sel_red,2),:),vertices(data2.F(new_sel_red,3),:));
                new_area_gray = triangleArea3d(vertices(data2.F(new_sel_gray,1),:),vertices(data2.F(new_sel_gray,2),:),vertices(data2.F(new_sel_gray,3),:));
                new_area_total = (sum(new_area_blue)+sum(new_area_red)+sum(new_area_gray)); % mm2
                new_area_red_total = sum(new_area_red);% mm2
                new_area_blue_total = sum(new_area_blue);% mm2
                
                disp(['total area ' num2str(i) ' = ' num2str(new_area_total/100) ' cm2'])
                disp(['red area ' num2str(i) '= ' num2str(new_area_red_total/100) ' cm2'])
                disp(['blue area ' num2str(i) ' = ' num2str(new_area_blue_total/100) ' cm2'])
                
                % save in the Excel file
                indices{1,i+1} = i;
                indices{2,i+1} = num2str(new_area_total/100);
                indices{3,i+1} = num2str(new_area_red_total/100);
                indices{4,i+1} = num2str(new_area_blue_total/100);
                
                clear new_heatmap mask new_sel_blue new_sel_red new_sel_gray new_area_blue new_area_red new_area_gray new_area_total new_area_red_total new_area_blue_total region
                
                waitbar(i / nbROIs)
            end
            
            currDir=pwd;
            cd(strcat(PATHNAME,'\heat_map_higher_lower_masks'))
            % save in an Excel sheet
            xlswrite('higher_lower_areas',indices);
            saveas(f5,'ROIs','tif')
            cd(currDir)
            close(h1)
            close(f5)
        
        case 'Load existing'
            
            choice = questdlg('Do you want to...', ...
                'Load existing ROIs', ...
                'Recompute WSS in all existing ROIs','Modify one ROI','Modify all ROIs','Modify one ROI');
            % Handle response
            switch choice
                case 'Modify one ROI'
                    [FileName,MrstructPath1,FilterIndex] = uigetfile(PATHNAME,'Select the ROI you want to load');
                    currentDir=pwd;
                    cd(MrstructPath1);
                    masks=ls('mask*');
                    figure(f5)
                    for i=1:size(masks,1)
                        load(strcat(MrstructPath1,masks(i,:)));
                        hold on, plot([region(:,1);,region(1,1)],[region(:,2);region(1,2)])
                        text(sum(region(:,1))/size(region,1),sum(region(:,2))/size(region,1),strcat('ROI',num2str(i)),'fontweight','bold')
                        clear region
                    end
                    cd(currentDir);
                    
                    load(strcat(MrstructPath1,FileName));
                    
                    polyAAo = impoly(gca,region);
                    wait(polyAAo);
                    region = getPosition(polyAAo);
                    disp('saving, pausing')
                    save(strcat(MrstructPath1,FileName),'region');
                    pause
                    
                    h1 = waitbar(0,'ROI modified, updated calculation of higher and lower WSS areas in progress...');
                    
                    mask = inpolygon(x1, y1, region(:,1), region(:,2));
                    new_heatmap = heat_mapp.*mask;
                    new_heatmap(~mask)=3;
                    new_sel_blue = (new_heatmap == 0);
                    new_sel_red = (new_heatmap == 1);
                    new_sel_gray = (new_heatmap == 2);
                    new_sel_blue = sum([new_sel_blue(data2.F(:,1)) new_sel_blue(data2.F(:,2)) new_sel_blue(data2.F(:,3))],2)>1;
                    new_sel_red = sum([new_sel_red(data2.F(:,1)) new_sel_red(data2.F(:,2)) new_sel_red(data2.F(:,3))],2)>1;
                    new_sel_gray = sum([new_sel_gray(data2.F(:,1)) new_sel_gray(data2.F(:,2)) new_sel_gray(data2.F(:,3))],2)>1;
                    new_area_blue = triangleArea3d(vertices(data2.F(new_sel_blue,1),:),vertices(data2.F(new_sel_blue,2),:),vertices(data2.F(new_sel_blue,3),:));
                    new_area_red = triangleArea3d(vertices(data2.F(new_sel_red,1),:),vertices(data2.F(new_sel_red,2),:),vertices(data2.F(new_sel_red,3),:));
                    new_area_gray = triangleArea3d(vertices(data2.F(new_sel_gray,1),:),vertices(data2.F(new_sel_gray,2),:),vertices(data2.F(new_sel_gray,3),:));
                    new_area_total = (sum(new_area_blue)+sum(new_area_red)+sum(new_area_gray)); % mm2
                    new_area_red_total = sum(new_area_red);% mm2
                    new_area_blue_total = sum(new_area_blue);% mm2
                    
                    disp(['new total area ' '= ' num2str(new_area_total/100) ' cm2'])
                    disp(['new red area ' '= ' num2str(new_area_red_total/100) ' cm2'])
                    disp(['new blue area ' '= ' num2str(new_area_blue_total/100) ' cm2'])
                    
                    % save in the Excel file
                    new_indices{1,1} = num2str(new_area_total/100);
                    new_indices{2,1} = num2str(new_area_red_total/100);
                    new_indices{3,1} = num2str(new_area_blue_total/100);
                    
                    clear new_heatmap mask new_sel_blue new_sel_red new_sel_gray new_area_blue new_area_red new_area_gray new_area_total new_area_red_total new_area_blue_total region
                    
                    currDir=pwd;
                    cd(strcat(PATHNAME,'\heat_map_higher_lower_masks'))
                    indMask = strfind(FileName, 'mask');
                    indMat = strfind(FileName, '.mat');
                    col = char(str2num(FileName(indMask+4:indMat-1))+'A');
                    xlRange = strcat([col '2:' col '4']);
                    % save in an Excel sheet
                    xlswrite('higher_lower_areas',new_indices,xlRange);
                    saveas(f5,'ROIs','tif')
                    cd(currDir)
                    waitbar(1)
                    close(h1)
                    close(f5)
                            
            case 'Recompute WSS in all existing ROIs'
                
                h1 = waitbar(0,'ROIs loading, updated calculation of higher and lower WSS areas in progress...');
                
                currentDir=pwd;
                cd(strcat(PATHNAME,'\heat_map_higher_lower_masks'))
                masks=ls('mask*');
                figure(f5)
                for i=1:size(masks,1)
                    load(strcat([PATHNAME '\heat_map_higher_lower_masks\mask' num2str(i)]));
                    hold on, plot([region(:,1);,region(1,1)],[region(:,2);region(1,2)])
                    clear region
                end
                cd(currentDir);
                
                % save in Excel
                indices{1,1} = 'ROI #';
                indices{2,1} = 'total area (cm2)';
                indices{3,1} = 'red area (cm2)';
                indices{4,1} = 'blue area (cm2)';
                
                for i=1:size(masks,1)
                    
                    load(strcat(PATHNAME, '\heat_map_higher_lower_masks\mask', num2str(i)));
                    mask = inpolygon(x1, y1, region(:,1), region(:,2));
                    new_heatmap = heat_mapp.*mask;
                    new_heatmap(~mask)=3;
                    new_sel_blue = (new_heatmap == 0);
                    new_sel_red = (new_heatmap == 1);
                    new_sel_gray = (new_heatmap == 2);
                    new_sel_blue = sum([new_sel_blue(data2.F(:,1)) new_sel_blue(data2.F(:,2)) new_sel_blue(data2.F(:,3))],2)>1;
                    new_sel_red = sum([new_sel_red(data2.F(:,1)) new_sel_red(data2.F(:,2)) new_sel_red(data2.F(:,3))],2)>1;
                    new_sel_gray = sum([new_sel_gray(data2.F(:,1)) new_sel_gray(data2.F(:,2)) new_sel_gray(data2.F(:,3))],2)>1;
                    new_area_blue = triangleArea3d(vertices(data2.F(new_sel_blue,1),:),vertices(data2.F(new_sel_blue,2),:),vertices(data2.F(new_sel_blue,3),:));
                    new_area_red = triangleArea3d(vertices(data2.F(new_sel_red,1),:),vertices(data2.F(new_sel_red,2),:),vertices(data2.F(new_sel_red,3),:));
                    new_area_gray = triangleArea3d(vertices(data2.F(new_sel_gray,1),:),vertices(data2.F(new_sel_gray,2),:),vertices(data2.F(new_sel_gray,3),:));
                    new_area_total = (sum(new_area_blue)+sum(new_area_red)+sum(new_area_gray)); % mm2
                    new_area_red_total = sum(new_area_red);% mm2
                    new_area_blue_total = sum(new_area_blue);% mm2
                    
                    disp(['total area ' num2str(i) ' = ' num2str(new_area_total/100) ' cm2'])
                    disp(['red area ' num2str(i) '= ' num2str(new_area_red_total/100) ' cm2'])
                    disp(['blue area ' num2str(i) ' = ' num2str(new_area_blue_total/100) ' cm2'])
                    
                    % save in the Excel file
                    indices{1,i+1} = i;
                    indices{2,i+1} = num2str(new_area_total/100);
                    indices{3,i+1} = num2str(new_area_red_total/100);
                    indices{4,i+1} = num2str(new_area_blue_total/100);
                    
                    clear new_heatmap mask new_sel_blue new_sel_red new_sel_gray new_area_blue new_area_red new_area_gray new_area_total new_area_red_total new_area_blue_total region
                    
                    waitbar (i/size(masks,1));
                end                
                
                currDir=pwd;
                cd(strcat(PATHNAME,'\heat_map_higher_lower_masks'))
                % save in an Excel sheet
                xlswrite('higher_lower_areas',indices);
                cd(currDir)
                close(h1)
                
            case 'Modify all ROIs'
                    
                currentDir=pwd;
                cd(strcat([MrstructPath '\..' '\heat_map_higher_lower_masks']));
                masks=ls('mask*');
                figure(f5)
                for i=1:size(masks,1)
                    load(strcat([MrstructPath '\..' '\heat_map_higher_lower_masks\mask' num2str(i)]));
                    hold on, plot([region(:,1);,region(1,1)],[region(:,2);region(1,2)])
                    clear region
                end
                for i=1:size(masks,1)
                    load(strcat([MrstructPath '\..' '\heat_map_higher_lower_masks\mask' num2str(i)]));
                    %Polygon and mask
                    polyAAo = impoly(gca,region);
                    wait(polyAAo);
                    region = getPosition(polyAAo);
                    disp('saving, pausing')
                    save(strcat([MrstructPath '\..' '\heat_map_higher_lower_masks\mask' num2str(i)]),'region');
                    text(sum(region(:,1))/size(region,1),sum(region(:,2))/size(region,1),strcat('ROI',num2str(i)),'fontweight','bold')
                    clear region
                    pause
                end
                cd(currentDir);
                h1 = waitbar(0,'ROIs modified, updated calculation of higher and lower WSS areas in progress...');
                
                % save in Excel
                indices{1,1} = 'ROI #';
                indices{2,1} = 'total area (cm2)';
                indices{3,1} = 'red area (cm2)';
                indices{4,1} = 'blue area (cm2)';
                
                for i=1:size(masks,1)
                    
                    load(strcat(PATHNAME, '\heat_map_higher_lower_masks\mask', num2str(i)));
                    mask = inpolygon(x1, y1, region(:,1), region(:,2));
                    new_heatmap = heat_mapp.*mask;
                    new_heatmap(~mask)=3;
                    new_sel_blue = (new_heatmap == 0);
                    new_sel_red = (new_heatmap == 1);
                    new_sel_gray = (new_heatmap == 2);
                    new_sel_blue = sum([new_sel_blue(data2.F(:,1)) new_sel_blue(data2.F(:,2)) new_sel_blue(data2.F(:,3))],2)>1;
                    new_sel_red = sum([new_sel_red(data2.F(:,1)) new_sel_red(data2.F(:,2)) new_sel_red(data2.F(:,3))],2)>1;
                    new_sel_gray = sum([new_sel_gray(data2.F(:,1)) new_sel_gray(data2.F(:,2)) new_sel_gray(data2.F(:,3))],2)>1;
                    new_area_blue = triangleArea3d(vertices(data2.F(new_sel_blue,1),:),vertices(data2.F(new_sel_blue,2),:),vertices(data2.F(new_sel_blue,3),:));
                    new_area_red = triangleArea3d(vertices(data2.F(new_sel_red,1),:),vertices(data2.F(new_sel_red,2),:),vertices(data2.F(new_sel_red,3),:));
                    new_area_gray = triangleArea3d(vertices(data2.F(new_sel_gray,1),:),vertices(data2.F(new_sel_gray,2),:),vertices(data2.F(new_sel_gray,3),:));
                    new_area_total = (sum(new_area_blue)+sum(new_area_red)+sum(new_area_gray)); % mm2
                    new_area_red_total = sum(new_area_red);% mm2
                    new_area_blue_total = sum(new_area_blue);% mm2
                    
                    disp(['total area ' num2str(i) ' = ' num2str(new_area_total/100) ' cm2'])
                    disp(['red area ' num2str(i) '= ' num2str(new_area_red_total/100) ' cm2'])
                    disp(['blue area ' num2str(i) ' = ' num2str(new_area_blue_total/100) ' cm2'])
                    
                    % save in the Excel file
                    indices{1,i+1} = i;
                    indices{2,i+1} = num2str(new_area_total/100);
                    indices{3,i+1} = num2str(new_area_red_total/100);
                    indices{4,i+1} = num2str(new_area_blue_total/100);
                    
                    clear new_heatmap mask new_sel_blue new_sel_red new_sel_gray new_area_blue new_area_red new_area_gray new_area_total new_area_red_total new_area_blue_total region
                    
                    waitbar(i / size(masks,1))
                end
                
                currDir=pwd;
                cd(strcat(MrstructPath,'\..','\heat_map_higher_lower_masks'))
                % save in an Excel sheet
                xlswrite('higher_lower_areas',indices);
                saveas(f5,'ROIs','tif')
                cd(currDir)
                close(h1)
                close(f5)

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  TRAFFIC LIGHT MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count1 = 0;
angles(1) = 0;
f1 = figure('Name','Traffic Light Map');
contours = zeros(size(L2));
contours(L2==0) = -1;
contours(L2==1) = 1;
%[F,V] = isosurface(smooth3(contours),0); % make a surface from the detected contours
[F,V] = isosurface(contours,0);
[F,V] = SmoothLaplacian(F,V,15);
patch('Faces',F,'Vertices',[V(:,1) V(:,2) V(:,3)],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.1);
hold on
[F1,V1] = isosurface(x./mask2_vox(1),y./mask2_vox(2),z./mask2_vox(3),smooth3(new_mask_red),0);
[F2,V2] = isosurface(x./mask2_vox(1),y./mask2_vox(2),z./mask2_vox(3),smooth3(new_mask_yellow),0);
[F3,V3] = isosurface(x./mask2_vox(1),y./mask2_vox(2),z./mask2_vox(3),smooth3(new_mask_blue),0);
[F4,V4] = isosurface(x./mask2_vox(1),y./mask2_vox(2),z./mask2_vox(3),smooth3(new_mask_green),0);
p11=patch('Faces',F1,'Vertices',V1,'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',1);
p12=patch('Faces',F2,'Vertices',V2,'EdgeColor','none','FaceColor',[1 0.9 0],'FaceAlpha',1);
set(p12,'HandleVisibility','on','Visible','off');
p13=patch('Faces',F3,'Vertices',V3,'EdgeColor','none','FaceColor',[0 0 1],'FaceAlpha',1);
%set(p13,'HandleVisibility','on','Visible','off')
p14=patch('Faces',F4,'Vertices',V4,'EdgeColor','none','FaceColor',[0 1 0],'FaceAlpha',1);
set(p14,'HandleVisibility','on','Visible','off')

aspectRatio = 1./mask2_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))

camlight headlight;camlight(180,0); lighting phong
% set up results folder
dir_orig = pwd;
dir_new = PATHNAME; %cd(dir_new); %cd('..')
%dir_new = pwd;
mkdir(PATHNAME,'results_traffic_light_map');
dir_new = strcat(dir_new,'\results_traffic_light_map');
saveas(gcf,[dir_new '\traffic_light_map.fig']);
load(strcat(MrstructPath,'\',FILENAME4))
magnitude = flipdim(double(mrStruct.dataAy(:,:,:,time)),3);clear mrStruct
magnitude(magnitude == 0) = 3;
magnitude(magnitude == 1) = 3;
magnitude(magnitude == 2) = 3;
hold on
s1 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3),magnitude(:,:,1),'EdgeColor','none');
set(s1,'HandleVisibility','off','Visible','off');
axis equal;colormap(gray)
view([180 -90])
caxis([0 64]);
aspectRatio = 1./mask2_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))
camlight(-45,0); lighting phong
axis ij; axis off
%camlight(-90,0); lighting phong
text(min(x(:)./mask2_vox(1)),max(y(:)./mask2_vox(2)-35),['Red volume: ' num2str(round(red_volume/1000)) ' cm^{3}' ])
text(min(x(:)./mask2_vox(1)),max(y(:)./mask2_vox(2)-30),['Red volume: ' num2str(round(percentage_red_volume)) ' %' ])
text(min(x(:)./mask2_vox(1)),max(y(:)./mask2_vox(2)-25),['Blue volume: ' num2str(round(blue_volume/1000)) ' cm^{3}' ])
text(min(x(:)./mask2_vox(1)),max(y(:)./mask2_vox(2)-20),['Blue volume: ' num2str(round(percentage_blue_volume)) ' %' ])
print(f1,'-djpeg','-r600',strcat(dir_new,'\traffic_light_map_front.jpg'));
axis ij; axis off
view([0 90]);camlight(45,0);%axis ij; %axis vis3d
print(f1,'-djpeg','-r600',strcat(dir_new,'\traffic_light_map_back.jpg'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PvO: Images for the Brief Report for NEJM were created with:
%%% 1. NO axis ij: I couldn't get the lighting right with axis ij on
%%%    The images therefore had to be flipped horizontally for the paper.
%%% 2. view([0 -90])
%%% 3. camlight(-90,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uicontrol('Style','text',...
    'Position',[10 375 120 20],...
    'String','Rotate')
uicontrol('Style', 'slider',...
    'Min',-90,'Max',90,'Value',0,...
    'Position', [10 350 120 20],...
    'SliderStep',[1/180 10/180],...
    'Callback', {@rotater1,gca});

uicontrol('Style','text',...
    'Position',[10 325 120 20],...
    'String','Rotate Light')
uicontrol('Style', 'slider',...
    'Min',-180,'Max',180,'Value',0,...
    'Position', [10 300 120 20],...
    'SliderStep',[1/360 120/360],...
    'Callback', {@rotate_light,gca});

uicontrol('Style','text',...
    'Position',[15 275 120 20],...
    'String','Show Red Map')
uicontrol('Style','checkbox',...
    'Value',1, 'Position', [15 275 20 20], ...
    'Callback', {@show_red_region,gca});

uicontrol('Style','text',...
    'Position',[15 250 120 20],...
    'String','Show Yellow Map')
uicontrol('Style','checkbox',...
    'Value',0, 'Position', [10 250 20 20], ...
    'Callback', {@show_yellow_region,gca});

uicontrol('Style','text',...
    'Position',[15 225 120 20],...
    'String','Show Blue Map')
uicontrol('Style','checkbox',...
    'Value',0, 'Position', [10 225 20 20], ...
    'Callback', {@show_blue_region,gca});

uicontrol('Style','text',...
    'Position',[15 200 120 20],...
    'String','Show Green Map')
uicontrol('Style','checkbox',...
    'Value',0, 'Position', [10 200 20 20], ...
    'Callback', {@show_green_region,gca});

uicontrol('Style','text',...
    'Position',[15 175 120 20],...
    'String','Show Anatomy')
uicontrol('Style','checkbox',...
    'Value',0, 'Position', [10 175 20 20], ...
    'Callback', {@show_anatomy1,gca});

uicontrol('Style','text',...
    'Position',[10 125 120 20],...
    'String','Contrast')
uicontrol('Style', 'slider',...
    'Min',0,'Max',max(magnitude(:)),'Value',1,...
    'Position', [10 100 120 20],...
    'SliderStep',[1/max(magnitude(:)) 1/max(magnitude(:))],...
    'Callback', {@change_contrast1,gca});

uicontrol('Style','text',...
    'Position',[10 75 120 20],...
    'String','Slice slider')
uicontrol('Style', 'slider',...
    'Min',1,'Max',size(magnitude,3),'Value',1,...
    'Position', [10 50 120 20],...
    'SliderStep',[1/(size(magnitude,3)-1) 10/(size(magnitude,3)-1)],...
    'Callback', {@move_slice1,gca});

    function move_slice1(hObj,event,ax)
        sliceobj = findobj(s1);
        delete(sliceobj)
        slice = round(get(hObj,'Value'));
        s1 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]).*(size(magnitude,3)-(slice-1)),magnitude(:,:,slice),'EdgeColor','none');
    end

    function show_red_region(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
            patchobj = findobj(p11);
            set(patchobj,'HandleVisibility','on','Visible','on');
        elseif show == 0
            patchobj = findobj(p11);
            set(patchobj,'HandleVisibility','off','Visible','off');
        end
    end

    function show_yellow_region(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
            patchobj = findobj(p12);
            set(patchobj,'HandleVisibility','on','Visible','on');
        elseif show == 0
            patchobj = findobj(p12);
            set(patchobj,'HandleVisibility','off','Visible','off');
        end
    end

    function show_blue_region(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
            patchobj = findobj(p13);
            set(patchobj,'HandleVisibility','on','Visible','on');
        elseif show == 0
            patchobj = findobj(p13);
            set(patchobj,'HandleVisibility','off','Visible','off');
        end
    end

    function show_green_region(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
            patchobj = findobj(p14);
            set(patchobj,'HandleVisibility','on','Visible','on');
        elseif show == 0
            patchobj = findobj(p14);
            set(patchobj,'HandleVisibility','off','Visible','off');
        end
    end

    function show_anatomy1(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
            patchobj = findobj(s1);
            set(patchobj,'HandleVisibility','on','Visible','on');
        elseif show == 0
            patchobj = findobj(s1);
            set(patchobj,'HandleVisibility','off','Visible','off');
        end
    end

    function change_contrast1(hObj,event,ax)
        contrast1 = round(get(hObj,'Value'));
        caxis([0 contrast1])
    end

    function rotater1(hObj,event,ax)
        count1 = count1 + 1;
        dtheta = get(hObj,'Value');
        dphi = 0;
        
        if count1 == 1;
            dtheta2 = dtheta*3;
        else
            dtheta2 = (dtheta - angles(count1-1))*3;
        end
        
        camorbit(dtheta2,dphi,'data',[0 1 0])
        angles(count1) = dtheta;
    end

    function rotate_light(hObj,event,ax)
        clear c
        dthetas = get(hObj,'Value');
        dphi = 0;
        
        c=camlight(dthetas,dphi);
        lighting phong
    end

set(f1,'toolbar','figure');

traffic_light.vox = mask2_vox;
traffic_light.mask = mask2;
traffic_light.red = new_mask_red;
traffic_light.yellow = new_mask_yellow;
traffic_light.blue = new_mask_blue;
traffic_light.green = new_mask_green;
traffic_light.vertices = [V(:,1) V(:,2) V(:,3)];
traffic_light.faces = F;

% save results in results folder
save(strcat(dir_new,'\results_trafficlight_map'),'traffic_light');
%savefig(f2,strcat(dir_new,'\heat_map'))
%cd(dir_orig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Heat map (WSS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count3 = 0;
angles(1) = 0;
f2 = figure('Name','Heat map');
x = data2.x_coor_wss;%./mask2_vox(1);
y = data2.y_coor_wss;%./mask2_vox(2);
z = data2.z_coor_wss;%./mask2_vox(3);
p2=patch('Faces',data2.F,'Vertices',[x y z],'EdgeColor','none', 'FaceVertexCData',heat_mapp,'FaceColor','interp','FaceAlpha',1);
sel_blue = (heat_mapp == 0);
sel_red = (heat_mapp == 1);
sel_gray = (heat_mapp == 2);
vertices = [x y z];
sel_blue = sum([sel_blue(data2.F(:,1)) sel_blue(data2.F(:,2)) sel_blue(data2.F(:,3))],2)>1;
sel_red = sum([sel_red(data2.F(:,1)) sel_red(data2.F(:,2)) sel_red(data2.F(:,3))],2)>1;
sel_gray = sum([sel_gray(data2.F(:,1)) sel_gray(data2.F(:,2)) sel_gray(data2.F(:,3))],2)>1;
area_blue = triangleArea3d(vertices(data2.F(sel_blue,1),:),vertices(data2.F(sel_blue,2),:),vertices(data2.F(sel_blue,3),:));
area_red = triangleArea3d(vertices(data2.F(sel_red,1),:),vertices(data2.F(sel_red,2),:),vertices(data2.F(sel_red,3),:));
area_gray = triangleArea3d(vertices(data2.F(sel_gray,1),:),vertices(data2.F(sel_gray,2),:),vertices(data2.F(sel_gray,3),:));
area_total = triangleArea3d(vertices(data2.F(:,1),:),vertices(data2.F(:,2),:),vertices(data2.F(:,3),:));
%hold on
% figure('Name','Area map')
% p5=patch('Faces',data2.F,'Vertices',[x y z],'EdgeColor','none','CData',area_total.*sel_red,'FaceColor','flat','FaceAlpha',1);
area_total = (sum(area_blue)+sum(area_red)+sum(area_gray)); % mm2
area_red_total = sum(area_red);% mm2
area_blue_total = sum(area_blue);% mm2

disp(['Total area = ' num2str(area_total/100) ' cm2'])
disp(['Red area = ' num2str(area_red_total/100) ' cm2'])
disp(['Blue area = ' num2str(area_blue_total/100) ' cm2'])

area_red_total_perc = area_red_total/area_total*100;
area_blue_total_perc = area_blue_total/area_total*100;
%return
colormap(color2);
caxis([0 64]);
axis equal; axis ij; axis off;
aspectRatio = 1./mask2_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))
view([-180 -90]);
% set up results folder
dir_orig = pwd;
dir_new = PATHNAME; %cd(dir_new); %cd('..')
%dir_new = pwd;
mkdir(PATHNAME,'results_heatmap');
dir_new = strcat(dir_new,'\results_heatmap');
saveas(gcf,[dir_new '\heat_map.fig'])
load(strcat(MrstructPath,'\',FILENAME4))
magnitude = flipdim(double(mrStruct.dataAy(:,:,:,time)),3);
magnitude(magnitude == 0) = 3;
magnitude(magnitude == 1) = 3;
magnitude(magnitude == 2) = 3;
hold on
s2 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3),magnitude(:,:,1),'EdgeColor','none');
set(s2,'HandleVisibility','off','Visible','off');
axis equal;
%view([-180 -90])
aspectRatio = 1./mask2_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))
text(min(x(:))-40,max(y(:))-10,['Red area: ' num2str(round(percentage_significant_higher_than_controls)) '%' ])
text(min(x(:))-40,max(y(:)),['Blue area: ' num2str(round(percentage_significant_lower_than_controls)) '%' ])
%print(f2,'-dtif','-r600',strcat(dir_new,'\heat_map_front.tif'));
print(f2,'-djpeg','-r600',strcat(dir_new,'\heat_map_front.jpg'));
axis equal; axis ij; axis off;%axis vis3d
view([0 90]);
print(f2,'-djpeg','-r600',strcat(dir_new,'\heat_map_back.jpg'));
aspectRatio = 1./mask2_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))

uicontrol('Style','text',...
    'Position',[10 375 120 20],...
    'String','Rotate')
uicontrol('Style', 'slider',...
    'Min',-90,'Max',90,'Value',0,...
    'Position', [10 350 120 20],...
    'SliderStep',[1/180 10/180],...
    'Callback', {@rotater2,gca});

uicontrol('Style','text',...
    'Position',[10 125 120 20],...
    'String','Contrast')
uicontrol('Style', 'slider',...
    'Min',0,'Max',max(magnitude(:)),'Value',1,...
    'Position', [10 100 120 20],...
    'SliderStep',[1/max(magnitude(:)) 1/max(magnitude(:))],...
    'Callback', {@change_contrast2,gca});

uicontrol('Style','text',...
    'Position',[10 75 120 20],...
    'String','Slice slider')
uicontrol('Style', 'slider',...
    'Min',1,'Max',size(magnitude,3),'Value',1,...
    'Position', [10 50 120 20],...
    'SliderStep',[1/(size(magnitude,3)-1) 10/(size(magnitude,3)-1)],...
    'Callback', {@move_slice2,gca});

uicontrol('Style','text',...
    'Position',[15 225 120 20],...
    'String','Show Anatomy')
uicontrol('Style','checkbox',...
    'Value',0, 'Position', [10 225 20 20], ...
    'Callback', {@show_anatomy2,gca});

% uicontrol('Style','text',...
%     'Position',[15 300 120 20],...
%     'String','Show Heat Map')
% uicontrol('Style','checkbox',...
%     'Value',1, 'Position', [15 300 20 20], ...
%     'Callback', {@show_heat_mapp,gca});

    function move_slice2(hObj,event,ax)
        sliceobj = findobj(s2);
        delete(sliceobj)
        slice = round(get(hObj,'Value'));
        s2 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]).*(size(magnitude,3)-(slice-1)),magnitude(:,:,slice),'EdgeColor','none');
    end

%     function show_heat_mapp(hObj,event,ax)
%         show = round(get(hObj,'Value'));
%         if show == 1
%             patchobj = findobj(p2);
%             set(patchobj,'HandleVisibility','on','Visible','on');
%         elseif show == 0
%             patchobj = findobj(p2);
%             set(patchobj,'HandleVisibility','off','Visible','off');
%         end
%     end

    function show_anatomy2(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
%             delete(p2)
%             x = data2.x_coor_wss./mask2_vox(1);
%             y = data2.y_coor_wss./mask2_vox(2);
%             z = data2.z_coor_wss./mask2_vox(3);
%             p2=patch('Faces',data2.F,'Vertices',[x y z],'EdgeColor','none', 'FaceVertexCData',heat_mapp,'FaceColor','interp','FaceAlpha',1);
%             
            patchobj = findobj(s2);
            set(patchobj,'HandleVisibility','on','Visible','on');
%            caxis([0 64])            
        elseif show == 0
            
%             delete(p2)
%             x = data2.x_coor_wss;
%             y = data2.y_coor_wss;
%             z = data2.z_coor_wss;
%             p2=patch('Faces',data2.F,'Vertices',[x y z],'EdgeColor','none', 'FaceVertexCData',heat_mapp,'FaceColor','interp','FaceAlpha',1);
%                 
            patchobj = findobj(s2);
            set(patchobj,'HandleVisibility','off','Visible','off');
 %           caxis([0 64])
        end
    end

    function change_contrast2(hObj,event,ax)
        contrast = round(get(hObj,'Value'))
        caxis([0 contrast])
    end

    function rotater2(hObj,event,ax)
        count3 = count3 + 1;
        dtheta3 = get(hObj,'Value');
        dphi = 0;
        
        if count3 == 1;
            dtheta4 = dtheta3*3;
        else
            dtheta4 = (dtheta3 - angles(count3-1))*3;
        end
        
        camorbit(dtheta4,dphi,'data',[0 1 0])
        angles(count3) = dtheta3;
    end

set(f2,'toolbar','figure');
axis vis3d

heat_map.heat_map = heat_mapp;
heat_map.vertices = [x y z];
heat_map.faces = data2.F;
heat_map.color = color1;

% save results in results folder
save(strcat(dir_new,'\heat_map'),'heat_map');
%savefig(f2,strcat(dir_new,'\heat_map'))
%cd(dir_orig)

% % Emilie: segmented red areas (ISMRM abstract)
% figure, patch('Faces',data2.F,'Vertices',[x y z], ...
%     'EdgeColor','none', 'FaceVertexCData',heat_mapp,'FaceColor','interp','FaceAlpha',1);
% axis equal;axis off; axis ij
% view([-180 -90])
% colormap(color2);
% caxis([0 64]);
% axis equal; axis ij; axis off;
% if exist(strcat(PATHNAME,'\heat_map_higher_lower_masks\mask1.mat'),'file')~=2
%     
%     mkdir(PATHNAME,'heat_map_higher_lower_masks')
%     for i = 1:3
%         %Polygon and mask for AAo
%         polyAAo = impoly;
%         wait(polyAAo);
%         region = getPosition(polyAAo);
%         %mask = inpolygon(x, y, AAo(:,1), AAo(:,2));
%         
%         disp('saving, pausing')
%         save(strcat([PATHNAME '\heat_map_higher_lower_masks\mask' num2str(i)]),'region');
%         pause
%     end
% 
% end
% load(strcat(PATHNAME,'\heat_map_higher_lower_masks\mask1'));
% hold on, plot([region(:,1);,region(1,1)],[region(:,2);region(1,2)])
% mask_AAo = inpolygon(x,y, region(:,1), region(:,2));
% new_heatmap = heat_mapp.*mask_AAo;
% new_heatmap(~mask_AAo)=3;
% new_sel_blue = (new_heatmap == 0);
% new_sel_red = (new_heatmap == 1);
% new_sel_gray = (new_heatmap == 2);
% new_sel_blue = sum([new_sel_blue(data2.F(:,1)) new_sel_blue(data2.F(:,2)) new_sel_blue(data2.F(:,3))],2)>1;
% new_sel_red = sum([new_sel_red(data2.F(:,1)) new_sel_red(data2.F(:,2)) new_sel_red(data2.F(:,3))],2)>1;
% new_sel_gray = sum([new_sel_gray(data2.F(:,1)) new_sel_gray(data2.F(:,2)) new_sel_gray(data2.F(:,3))],2)>1;
% new_area_blue = triangleArea3d(vertices(data2.F(new_sel_blue,1),:),vertices(data2.F(new_sel_blue,2),:),vertices(data2.F(new_sel_blue,3),:));
% new_area_red = triangleArea3d(vertices(data2.F(new_sel_red,1),:),vertices(data2.F(new_sel_red,2),:),vertices(data2.F(new_sel_red,3),:));
% new_area_gray = triangleArea3d(vertices(data2.F(new_sel_gray,1),:),vertices(data2.F(new_sel_gray,2),:),vertices(data2.F(new_sel_gray,3),:));
% new_area_total = (sum(new_area_blue)+sum(new_area_red)+sum(new_area_gray)); % mm2
% new_area_red_total = sum(new_area_red);% mm2
% new_area_blue_total = sum(new_area_blue);% mm2
% 
% disp(['total area = ' num2str(new_area_total/100) ' cm2'])
% disp(['red area = ' num2str(new_area_red_total/100) ' cm2'])
% % disp(['New blue area = ' num2str(new_area_blue_total/100) ' cm2'])
% area_all(1)=new_area_red_total/100;
% area_all(2)=new_area_total/100;
% 
% load(strcat(PATHNAME,'\heat_map_higher_lower_masks\mask2'));
% hold on, plot([region(:,1);,region(1,1)],[region(:,2);region(1,2)])
% mask_arch = inpolygon(x,y, region(:,1), region(:,2));
% clear new_heatmap
% new_heatmap = heat_mapp.*mask_arch;
% new_heatmap(~mask_arch)=3;
% new_sel_blue = (new_heatmap == 0);
% new_sel_red = (new_heatmap == 1);
% new_sel_gray = (new_heatmap == 2);
% new_sel_blue = sum([new_sel_blue(data2.F(:,1)) new_sel_blue(data2.F(:,2)) new_sel_blue(data2.F(:,3))],2)>1;
% new_sel_red = sum([new_sel_red(data2.F(:,1)) new_sel_red(data2.F(:,2)) new_sel_red(data2.F(:,3))],2)>1;
% new_sel_gray = sum([new_sel_gray(data2.F(:,1)) new_sel_gray(data2.F(:,2)) new_sel_gray(data2.F(:,3))],2)>1;
% new_area_blue = triangleArea3d(vertices(data2.F(new_sel_blue,1),:),vertices(data2.F(new_sel_blue,2),:),vertices(data2.F(new_sel_blue,3),:));
% new_area_red = triangleArea3d(vertices(data2.F(new_sel_red,1),:),vertices(data2.F(new_sel_red,2),:),vertices(data2.F(new_sel_red,3),:));
% new_area_gray = triangleArea3d(vertices(data2.F(new_sel_gray,1),:),vertices(data2.F(new_sel_gray,2),:),vertices(data2.F(new_sel_gray,3),:));
% new_area_total = (sum(new_area_blue)+sum(new_area_red)+sum(new_area_gray)); % mm2
% new_area_red_total = sum(new_area_red);% mm2
% new_area_blue_total = sum(new_area_blue);% mm2
% 
% disp(['total area arch = ' num2str(new_area_total/100) ' cm2'])
% disp(['red area arch = ' num2str(new_area_red_total/100) ' cm2'])
% area_all(3)=new_area_red_total/100;
% area_all(4)=new_area_total/100;
% 
% load(strcat(PATHNAME,'\heat_map_higher_lower_masks\mask3'));
% hold on, plot([region(:,1);,region(1,1)],[region(:,2);region(1,2)])
% mask_DAo = inpolygon(x,y, region(:,1), region(:,2));
% clear new_heatmap
% new_heatmap = heat_mapp.*mask_DAo;
% new_heatmap(~mask_DAo)=3;
% new_sel_blue = (new_heatmap == 0);
% new_sel_red = (new_heatmap == 1);
% new_sel_gray = (new_heatmap == 2);
% new_sel_blue = sum([new_sel_blue(data2.F(:,1)) new_sel_blue(data2.F(:,2)) new_sel_blue(data2.F(:,3))],2)>1;
% new_sel_red = sum([new_sel_red(data2.F(:,1)) new_sel_red(data2.F(:,2)) new_sel_red(data2.F(:,3))],2)>1;
% new_sel_gray = sum([new_sel_gray(data2.F(:,1)) new_sel_gray(data2.F(:,2)) new_sel_gray(data2.F(:,3))],2)>1;
% new_area_blue = triangleArea3d(vertices(data2.F(new_sel_blue,1),:),vertices(data2.F(new_sel_blue,2),:),vertices(data2.F(new_sel_blue,3),:));
% new_area_red = triangleArea3d(vertices(data2.F(new_sel_red,1),:),vertices(data2.F(new_sel_red,2),:),vertices(data2.F(new_sel_red,3),:));
% new_area_gray = triangleArea3d(vertices(data2.F(new_sel_gray,1),:),vertices(data2.F(new_sel_gray,2),:),vertices(data2.F(new_sel_gray,3),:));
% new_area_total = (sum(new_area_blue)+sum(new_area_red)+sum(new_area_gray)); % mm2
% new_area_red_total = sum(new_area_red);% mm2
% new_area_blue_total = sum(new_area_blue);% mm2
% 
% disp(['total area DAo = ' num2str(new_area_total/100) ' cm2'])
% disp(['red area DAo = ' num2str(new_area_red_total/100) ' cm2'])
% area_all(5)=new_area_red_total/100;
% area_all(6)=new_area_total/100;

if images_for_surgeryFlag
    f3 = figure('Name','Heat map');
    x = data2.x_coor_wss/mask2_vox(1);
    y = data2.y_coor_wss/mask2_vox(2);
    z = data2.z_coor_wss/mask2_vox(3);
    p3=patch('Faces',data2.F,'Vertices',[x y z],'EdgeColor','none', 'FaceVertexCData',heat_mapp,'FaceColor','interp','FaceAlpha',1);
    color3(1,:) = [0 0 1];
    color3(2,:) = [1 0 0];
    color3(3,:) = [0.5 0.5 0.5];
    color3(4:64,:) = gray_colormap(4:64,:);
    colormap(color3);
    caxis([0 64]);
    axis equal; axis ij; axis off;
    aspectRatio = 1./mask2_vox;
    set(gca,'dataaspectRatio',aspectRatio(1:3))
    view([-180 -90]);
    % set up results folder
    dir_orig = pwd;
    dir_new = PATHNAME; %cd(dir_new); %cd('..')
    %dir_new = pwd;
    mkdir([dir_new filesep 'images_for_surgery']);
    dir_new = strcat(dir_new,'\images_for_surgery');
    load(strcat(MrstructPath,'\',FILENAME4))
    magnitude = flipdim(double(mrStruct.dataAy(:,:,:,time)),3);
    magnitude(magnitude == 0) = 3;
    magnitude(magnitude == 1) = 3;
    magnitude(magnitude == 2) = 3;
    hold on
    s3 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3),magnitude(:,:,1),'EdgeColor','none');
    set(s3,'HandleVisibility','on','Visible','on');
    axis equal;
    %view([-180 -90])
    aspectRatio = 1./mask2_vox;
    set(gca,'dataaspectRatio',aspectRatio(1:3))
    print(f3,'-djpeg','-r600',strcat(dir_new,'\image1'));
    
    
    hax_f3 = gca;
    h = figure('Name', 'Report for surgery'); % Emilie: figure including all views for report
    subP1 = subplot(1,5,1);
    posP1 = get(subP1,'Position');
    delete(subP1);
    hax_h1 = copyobj(hax_f3,h); % copy figure 1 axes into report figure
    set(hax_h1, 'Position', posP1);
    colormap(color3);
    % Emilie: manual interaction to chose magnitude slice on which RPA can be visualized, if needed
    h2 = figure('Name','Selection of the magnitude slice');
    colormap(color3);
    caxis([0 64]);
    axis equal; axis ij; axis off;
    set(gca,'dataaspectRatio',aspectRatio(1:3))
    view([-180 -90]);
    hold on
    s4 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3)/2,magnitude(:,:,size(magnitude,3)/2),'EdgeColor','none');
    uicontrol('Style','text',...
        'Position',[10 200 120 70],...
        'String','Please choose using the slider a magnitude slice on which RPA can be visualized and then click ok button')
    uicontrol('Style','text',...
        'Position',[10 75 120 20],...
        'String','Slice slider')
    sl1 = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(magnitude,3),'Value',size(magnitude,3)/2,...
        'Position', [10 50 120 20],...
        'SliderStep',[1/(size(magnitude,3)-1) 10/(size(magnitude,3)-1)],...
        'Callback', {@move_slice3,gca});
    but_ok = uicontrol('Position',[10 25 120 20],...
        'String','ok',...
        'Callback', 'uiresume(gcbf)');
    uiwait(gcf);
    RPAslice=get(sl1, 'Value');
    close(h2)
    % End Emilie: manual interaction to chose magnitude slice on which RPA can be visualized
    figure(f3)
    
    axis equal; axis ij; axis off;axis vis3d
    delete(s3)
%     s3 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3)/2,magnitude(:,:,size(magnitude,3)/2),'EdgeColor','none');
    s3 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* RPAslice,magnitude(:,:,RPAslice),'EdgeColor','none');  % Emilie: display the slice previously selected by user
    set(s3,'HandleVisibility','on','Visible','on');
    aspectRatio = 1./mask2_vox;
    set(gca,'dataaspectRatio',aspectRatio(1:3))
    print(f3,'-djpeg','-r600',strcat(dir_new,'\image2'));
    % Emilie: copy figure 2 axes into report figure
    hax_f3 = gca;
    figure(h)
    subP2 = subplot(1,5,2);
    posP2 = get(subP2,'Position');
    delete(subP2);
    hax_h2 = copyobj(hax_f3,h);
    set(hax_h2, 'Position', posP2);
    figure(f3)
    camorbit(-90,0,'data',[0 1 0])
aspectRatio = 1./mask2_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))    
    print(f3,'-djpeg','-r600',strcat(dir_new,'\image3'));
    % Emilie: copy figure 3 axes into report figure
    hax_f3 = gca;
    figure(h)
    subP3 = subplot(1,5,3);
    posP3 = get(subP3,'Position');
    delete(subP3);
    hax_h3 = copyobj(hax_f3,h);
    set(hax_h3, 'Position', posP3);

    figure(f3)
    view([0 90])
aspectRatio = 1./mask2_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))    
    print(f3,'-djpeg','-r600',strcat(dir_new,'\image4'));
    % Emilie: copy figure 4 axes into report figure
    hax_f3 = gca;
    figure(h)
    subP4 = subplot(1,5,4);
    posP4=get(subP4,'Position');
    delete(subP4);
    hax_h4=copyobj(hax_f3,h);
    set(hax_h4, 'Position', posP4);
    %
    figure(f3)
    delete(s3)
    s3 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* 1,magnitude(:,:,size(magnitude,3)),'EdgeColor','none');
   aspectRatio = 1./mask2_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))
    print(f3,'-djpeg','-r600',strcat(dir_new,'\image5'));
    % Emilie: copy figure 5 axes into report figure
    hax_f3 = gca;
    figure(h)
    subP5 = subplot(1,5,5);
    posP5 = get(subP5,'Position');
    delete(subP5);
    hax_h5 = copyobj(hax_f3,h);
    set(hax_h5, 'Position', posP5);
    % Display of report figure
    posP1(1)=-0.075;
    posP1(3)=.38;
    set(hax_h1, 'Position', posP1);
    posP2(1)=posP2(1)-.1;
    posP2(3)=.3;
    set(hax_h2, 'Position', posP2);
    posP3(1)=posP3(1)-0.11;
    posP3(3)=.3;
    set(hax_h3, 'Position', posP3);
    posP4(1)=posP4(1)-.11;
    posP4(3)=.3;
    set(hax_h4, 'Position', posP4);
    posP5(1)=posP5(1)-.045;
    posP5(3)=.3;
    set(hax_h5, 'Position', posP5);
    % Save report figure
    print(h,'-djpeg','-r600',strcat(dir_new,'\report'));
    delete(f3)
end

function move_slice3(hObj,event,ax) % Emilie: for manual interaction to chose magnitude slice on which RPA can be visualized
	sliceobj = findobj(s4);
	delete(sliceobj)
	slice = round(get(hObj,'Value'));
	s4 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]).*(size(magnitude,3)-(slice-1)),magnitude(:,:,slice),'EdgeColor','none');
end

function move_slice4(hObj,event,ax) % Emilie: for manual interaction to chose magnitude slice on which RPA can be visualized
    sliceobj = findobj(s4);
    delete(sliceobj)
    slice = round(get(hObj,'Value'));
    s4 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]).*(size(magnitude,3)-(slice-1)),magnitude(:,:,slice),'EdgeColor','none');
end

end