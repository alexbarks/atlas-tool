function [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(offset,plotFlag,calculateRE_Flag,calculateIE_Flag,calculate_traffic_light_volumeFlag,calculate_area_of_higherlowerFlag,peak_systolicFlag)

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
% 1)offset          : To be able to calculate the registration error (RE, see paper
%                       "A Methodology to Detect Abnormal Relative Wall Shear Stress on the Full Surface of the Thoracic Aorta Using 4D Flow MRI
%                       van Ooij P, Potters WV, Nederveen AJ, Allen BD, Collins J, Carr J, Malaisrie SC, Markl M, Barker AJ
%                       Magn Reson Med. 2014 Feb 25, doi: 10.1002/mrm.25224. [Epub ahead of print])",
%                       the registered aorta needs to be transformed back to a matrix. This can only be done when all x, y and z coordinates after registration are > 0.
%                     Otherwise Matlab will return an error and all will be for nothing. We therefore need to put in an offset to prevent
%                     this from happening.
% 2)plotFlag        : If plotFlag switched, Matlab will output any possible image to the screen to check if everything happens correctly.
% 3)calculateRE_Flag: When switched on the registration error will be calculated. However, this is a different RE than the one in the paper
%                     mentioned above as the RE in this function is calculated from AFFINE registration whereas the RE reported in the paper
%                     is calculated from RIGID registration in the function 'make_geometry_point_cloud.m'
% 4)calculateIE_Flag: When switched on the interpolation error (RE, see paper mentioned above) will be calculated. Note that ROIs are needed
%                     which can be drawn manually when switched on.
% 5)calculate_area_of_significanceFlag: When switched on the area of significance will be calculated (see first paper mentioned above). Note that ROIs are
%                     needed which can be drawn manually when switched on. However, if calculateIE_Flag is switched on than you need to do this only once.
% 6)peak_systolicFlag: When switched on the comparison of the individual patient with the ensemble-averaged maps will be performed for the peak
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
% [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(offset,plotFlag,calculateRE_Flag,calculateIE_Flag,calculate_area_of_significanceFlag,peak_systolicFlag)
% [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(40,1,1,1,1,0)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all % I hate functions that don't close the stuff that was running before

if nargin < 1 || isempty(offset)
    offset = 40;
end

if nargin < 2 || isempty(plotFlag)
    plotFlag = 1;
end

if nargin < 3 || isempty(calculateRE_Flag)
    calculateRE_Flag = 1;
end

if nargin < 4 || isempty(calculateIE_Flag)
    calculateIE_Flag = 0;
end

if nargin < 6 || isempty(calculate_traffic_light_volumeFlag)
    calculate_traffic_light_volumeFlag = 1;
end

if nargin < 7 || isempty(calculate_area_of_higherlowerFlag)
    calculate_area_of_higherlowerFlag = 1;
end

if nargin < 8 || isempty(peak_systolicFlag)
    peak_systolicFlag = 0;
end

global data
global atlas
%
%data = [];
Rotation_Translation = [];

%[FILENAME, PATHNAME_atlas] = uigetfile('C:\1_Chicago\Data\MIMICS','Load atlas')
PATHNAME_atlas = 'c:\_ensightCases\bav_tissue\Controls\'
FILENAME = 'atlas.mat'
load(strcat(PATHNAME_atlas,FILENAME))
mask1 = atlas.mask;

if plotFlag == 1
    figure('Name','Mean atlas velocity')
    scatter3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,atlas.mean_vel.*50,atlas.mean_vel,'filled')
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
    
    figure('Name','std atlas velocity')
    scatter3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,atlas.std_vel.*50,atlas.std_vel,'filled')
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
    
    figure('Name','Mean atlas WSS')
    patch('Faces',atlas.faces,'Vertices',atlas.vertices,'EdgeColor','none', 'FaceVertexCData',atlas.mean_wss,'FaceColor','interp','FaceAlpha',1);colorbar;
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
    
    figure('Name','std atlas WSS')
    patch('Faces',atlas.faces,'Vertices',atlas.vertices,'EdgeColor','none', 'FaceVertexCData',atlas.std_wss,'FaceColor','interp','FaceAlpha',1);colorbar;
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
end

if calculateIE_Flag == 1;
    
    if ~exist(strcat(PATHNAME_atlas,'interpolation_error_ROI\mask1.mat'),'file')
        
        F1=figure('Name','Atlas shape: Paused after finishing a region so press space when finished!');
        patch('Faces',atlas.faces,'Vertices',atlas.vertices,'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.25);
        view([-180 -90]);axis ij;axis equal;axis off
        
        for i = 1:6
            %Polygon and mask for AAo
            polyAAo = impoly;
            wait(polyAAo);
            region = getPosition(polyAAo);
            
            %          disp('saving, pausing')
            mkdir(PATHNAME_atlas,'interpolation_error_ROI')
            save(strcat([PATHNAME_atlas 'interpolation_error_ROI\mask' num2str(i)]),'region');
            pause
        end
        
        close(F1)
    end
    
    load(strcat(PATHNAME_atlas,'interpolation_error_ROI\mask1'))
    atlas_mask_AAo_inner_vel = inpolygon(atlas.x_coor_vel, atlas.y_coor_vel, region(:,1), region(:,2));
    mean_vel_asc_inner = mean(atlas.mean_vel(atlas_mask_AAo_inner_vel));
    atlas_mask_AAo_inner_wss = inpolygon(atlas.x_coor_wss, atlas.y_coor_wss, region(:,1), region(:,2));
    mean_wss_asc_inner = mean(atlas.mean_wss(atlas_mask_AAo_inner_wss));
    load(strcat(PATHNAME_atlas,'interpolation_error_ROI\mask2'))
    atlas_mask_AAo_outer_vel = inpolygon(atlas.x_coor_vel, atlas.y_coor_vel, region(:,1), region(:,2));
    mean_vel_asc_outer = mean(atlas.mean_vel(atlas_mask_AAo_outer_vel));
    atlas_mask_AAo_outer_wss = inpolygon(atlas.x_coor_wss, atlas.y_coor_wss, region(:,1), region(:,2));
    mean_wss_asc_outer = mean(atlas.mean_wss(atlas_mask_AAo_outer_wss));
    load(strcat(PATHNAME_atlas,'interpolation_error_ROI\mask3'))
    atlas_mask_arch_inner_vel = inpolygon(atlas.x_coor_vel, atlas.y_coor_vel, region(:,1), region(:,2));
    mean_vel_arch_inner = mean(atlas.mean_vel(atlas_mask_arch_inner_vel));
    atlas_mask_arch_inner_wss = inpolygon(atlas.x_coor_wss, atlas.y_coor_wss, region(:,1), region(:,2));
    mean_wss_arch_inner = mean(atlas.mean_wss(atlas_mask_arch_inner_wss));
    load(strcat(PATHNAME_atlas,'interpolation_error_ROI\mask4'))
    atlas_mask_arch_outer_vel = inpolygon(atlas.x_coor_vel, atlas.y_coor_vel, region(:,1), region(:,2));
    mean_vel_arch_outer = mean(atlas.mean_vel(atlas_mask_arch_outer_vel));
    atlas_mask_arch_outer_wss = inpolygon(atlas.x_coor_wss, atlas.y_coor_wss, region(:,1), region(:,2));
    mean_wss_arch_outer = mean(atlas.mean_wss(atlas_mask_arch_outer_wss));
    load(strcat(PATHNAME_atlas,'interpolation_error_ROI\mask5'))
    atlas_mask_DAo_inner_vel = inpolygon(atlas.x_coor_vel, atlas.y_coor_vel, region(:,1), region(:,2));
    mean_vel_DAo_inner = mean(atlas.mean_vel(atlas_mask_DAo_inner_vel));
    atlas_mask_DAo_inner_wss = inpolygon(atlas.x_coor_wss, atlas.y_coor_wss, region(:,1), region(:,2));
    mean_wss_DAo_inner = mean(atlas.mean_wss(atlas_mask_DAo_inner_wss));
    load(strcat(PATHNAME_atlas,'interpolation_error_ROI\mask6'))
    atlas_mask_DAo_outer_vel = inpolygon(atlas.x_coor_vel, atlas.y_coor_vel, region(:,1), region(:,2));
    mean_vel_DAo_outer = mean(atlas.mean_vel(atlas_mask_DAo_outer_vel));
    atlas_mask_DAo_outer_wss = inpolygon(atlas.x_coor_wss, atlas.y_coor_wss, region(:,1), region(:,2));
    mean_wss_DAo_outer = mean(atlas.mean_wss(atlas_mask_DAo_outer_wss));
    
    mean_vel_before_interpolation(1,1) = mean_vel_asc_inner;
    mean_vel_before_interpolation(2,1) = mean_vel_asc_outer;
    mean_vel_before_interpolation(3,1) = mean_vel_arch_inner;
    mean_vel_before_interpolation(4,1) = mean_vel_arch_outer;
    mean_vel_before_interpolation(5,1) = mean_vel_DAo_inner;
    mean_vel_before_interpolation(6,1) = mean_vel_DAo_outer;
    
    mean_wss_before_interpolation(1,1) = mean_wss_asc_inner;
    mean_wss_before_interpolation(2,1) = mean_wss_asc_outer;
    mean_wss_before_interpolation(3,1) = mean_wss_arch_inner;
    mean_wss_before_interpolation(4,1) = mean_wss_arch_outer;
    mean_wss_before_interpolation(5,1) = mean_wss_DAo_inner;
    mean_wss_before_interpolation(6,1) = mean_wss_DAo_outer;
    
    if plotFlag == 1
        figure('Name','Velocity before interpolation: inner AAo')
        scatter3(atlas.x_coor_vel(atlas_mask_AAo_inner_vel),atlas.y_coor_vel(atlas_mask_AAo_inner_vel),atlas.z_coor_vel(atlas_mask_AAo_inner_vel),20,atlas.mean_vel(atlas_mask_AAo_inner_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','Velocity before interpolation: outer AAo')
        scatter3(atlas.x_coor_vel(atlas_mask_AAo_outer_vel),atlas.y_coor_vel(atlas_mask_AAo_outer_vel),atlas.z_coor_vel(atlas_mask_AAo_outer_vel),20,atlas.mean_vel(atlas_mask_AAo_outer_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','Velocity before interpolation: inner arch')
        scatter3(atlas.x_coor_vel(atlas_mask_arch_inner_vel),atlas.y_coor_vel(atlas_mask_arch_inner_vel),atlas.z_coor_vel(atlas_mask_arch_inner_vel),20,atlas.mean_vel(atlas_mask_arch_inner_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','Velocity before interpolation: outer arch')
        scatter3(atlas.x_coor_vel(atlas_mask_arch_outer_vel),atlas.y_coor_vel(atlas_mask_arch_outer_vel),atlas.z_coor_vel(atlas_mask_arch_outer_vel),20,atlas.mean_vel(atlas_mask_arch_outer_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','Velocity before interpolation: inner DAo')
        scatter3(atlas.x_coor_vel(atlas_mask_DAo_inner_vel),atlas.y_coor_vel(atlas_mask_DAo_inner_vel),atlas.z_coor_vel(atlas_mask_DAo_inner_vel),20,atlas.mean_vel(atlas_mask_DAo_inner_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','Velocity before interpolation: outer DAo')
        scatter3(atlas.x_coor_vel(atlas_mask_DAo_outer_vel),atlas.y_coor_vel(atlas_mask_DAo_outer_vel),atlas.z_coor_vel(atlas_mask_DAo_outer_vel),20,atlas.mean_vel(atlas_mask_DAo_outer_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS before interpolation: inner AAo')
        scatter3(atlas.x_coor_wss(atlas_mask_AAo_inner_wss),atlas.y_coor_wss(atlas_mask_AAo_inner_wss),atlas.z_coor_wss(atlas_mask_AAo_inner_wss),20,atlas.mean_wss(atlas_mask_AAo_inner_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS before interpolation: outer AAo')
        scatter3(atlas.x_coor_wss(atlas_mask_AAo_outer_wss),atlas.y_coor_wss(atlas_mask_AAo_outer_wss),atlas.z_coor_wss(atlas_mask_AAo_outer_wss),20,atlas.mean_wss(atlas_mask_AAo_outer_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS before interpolation: inner arch')
        scatter3(atlas.x_coor_wss(atlas_mask_arch_inner_wss),atlas.y_coor_wss(atlas_mask_arch_inner_wss),atlas.z_coor_wss(atlas_mask_arch_inner_wss),20,atlas.mean_wss(atlas_mask_arch_inner_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS before interpolation: outer arch')
        scatter3(atlas.x_coor_wss(atlas_mask_arch_outer_wss),atlas.y_coor_wss(atlas_mask_arch_outer_wss),atlas.z_coor_wss(atlas_mask_arch_outer_wss),20,atlas.mean_wss(atlas_mask_arch_outer_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS before interpolation: inner DAo')
        scatter3(atlas.x_coor_wss(atlas_mask_DAo_inner_wss),atlas.y_coor_wss(atlas_mask_DAo_inner_wss),atlas.z_coor_wss(atlas_mask_DAo_inner_wss),20,atlas.mean_wss(atlas_mask_DAo_inner_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS before interpolation: outer DAo')
        scatter3(atlas.x_coor_wss(atlas_mask_DAo_outer_wss),atlas.y_coor_wss(atlas_mask_DAo_outer_wss),atlas.z_coor_wss(atlas_mask_DAo_outer_wss),20,atlas.mean_wss(atlas_mask_DAo_outer_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
    end
end

[FILENAME, PATHNAME] = uigetfile('c:\_ensightCases\bav_tissue\','Load aorta data')
%PATHNAME = 'C:\1_Chicago\Data\MIMICS\traffic_light\2_RN\PT256-MW\'
%FILENAME = 'data_done'
load(strcat(PATHNAME,FILENAME))
data2 = data; clear data

if ~isfield(data2,'PC_unaliased')
    data2.PC_unaliased = data2.PC_zeros;
end
mask2 = (data2.PC_unaliased(:,:,:,1,1) ~= 0);

L2a = (mask2 ~= 0); % Using two logicals for the second mask since the size changes after the translation which is needed for coordinates but not for velocity

%%% translate the atlasmetry away from the origin to prevent coordinates < 0 after registration, otherwise the atlasmetry can not be transformed back to a matrix
sizes = [size(mask2,1)+offset size(mask2,2)+offset size(mask2,3)+offset];
mask2b = zeros(sizes);
mask2b((offset+1):size(mask2b,1),(offset+1):size(mask2b,2),(offset+1):size(mask2b,3)) = mask2;
mask2 = mask2b;clear mask2b

L2b = (mask2 ~= 0);
% create velocity coordinates
[x,y,z] = meshgrid((1:size(mask2,2)).* data2.vox(2), ...
    (1:size(mask2,1)).*data2.vox(1),(1:size(mask2,3)).* data2.vox(3));
data2.x_coor_vel = x(L2b);data2.y_coor_vel = y(L2b);data2.z_coor_vel = z(L2b);
contours = zeros(size(L2b));
contours(L2b==0) = -1;
contours(L2b==1) = 1;
[data2.F,V] = isosurface(contours,0); % make a surface from the detected contours
data2.V = V .* (ones(size(V,1),1) * data2.vox(1:3));

for t = 1:size(data2.PC_unaliased,5)
    vx = squeeze(data2.PC_unaliased(:,:,:,1,t));
    vy = squeeze(data2.PC_unaliased(:,:,:,2,t));
    vz = squeeze(data2.PC_unaliased(:,:,:,3,t));
    vmagn = sqrt(vx.^2 + vy.^2 + vz.^2);
    L =(vmagn ~= 0);
    mean_velo(t) = mean(vmagn(L));
end

if plotFlag == 1
    figure('Name','Mean velocity')
    plot(1:size(data2.PC_unaliased,5),mean_velo,'-ro','LineWidth',5,...
        'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',16);
end

[I,time] = find(mean_velo==max(mean_velo));

%%% Wall shear stress coordinates
data2.x_coor_wss = data2.V_matrix{1}(:,1) + (offset*data2.vox(1));
data2.y_coor_wss = data2.V_matrix{1}(:,2) + (offset*data2.vox(2));
data2.z_coor_wss = data2.V_matrix{1}(:,3) + (offset*data2.vox(3));

if peak_systolicFlag == 1
    % Peak systolic velocity
    data2.x_value_vel = data2.PC_unaliased(:,:,:,1,time);
    data2.y_value_vel = data2.PC_unaliased(:,:,:,2,time);
    data2.z_value_vel = data2.PC_unaliased(:,:,:,3,time);
    data2.x_value_vel = data2.x_value_vel(L2a);
    data2.y_value_vel = data2.y_value_vel(L2a);
    data2.z_value_vel = data2.z_value_vel(L2a);
    % Peak systolic WSS
    data2.x_value_wss = data2.WSS_matrix{time}(:,1);
    data2.y_value_wss = data2.WSS_matrix{time}(:,2);
    data2.z_value_wss = data2.WSS_matrix{time}(:,3);
    data2.x_value_wss = data2.x_value_wss;
    data2.y_value_wss = data2.y_value_wss;
    data2.z_value_wss = data2.z_value_wss;
elseif peak_systolicFlag == 0
    % Velocity averaged over 5 systolic time frames
    data2.x_value_vel_t1 = data2.PC_unaliased(:,:,:,1,time-2);data2.y_value_vel_t1 = data2.PC_unaliased(:,:,:,2,time-2);data2.z_value_vel_t1 = data2.PC_unaliased(:,:,:,3,time-2);
    data2.x_value_vel_t2 = data2.PC_unaliased(:,:,:,1,time-1);data2.y_value_vel_t2 = data2.PC_unaliased(:,:,:,2,time-1);data2.z_value_vel_t2 = data2.PC_unaliased(:,:,:,3,time-1);
    data2.x_value_vel_t3 = data2.PC_unaliased(:,:,:,1,time);  data2.y_value_vel_t3 = data2.PC_unaliased(:,:,:,2,time);  data2.z_value_vel_t3 = data2.PC_unaliased(:,:,:,3,time);
    data2.x_value_vel_t4 = data2.PC_unaliased(:,:,:,1,time+1);data2.y_value_vel_t4 = data2.PC_unaliased(:,:,:,2,time+1);data2.z_value_vel_t4 = data2.PC_unaliased(:,:,:,3,time+1);
    data2.x_value_vel_t5 = data2.PC_unaliased(:,:,:,1,time+2);data2.y_value_vel_t5 = data2.PC_unaliased(:,:,:,2,time+2);data2.z_value_vel_t5 = data2.PC_unaliased(:,:,:,3,time+2);
    data2.x_value_vel = (data2.x_value_vel_t1(L2a) + data2.x_value_vel_t2(L2a) + data2.x_value_vel_t3(L2a) + data2.x_value_vel_t4(L2a) + data2.x_value_vel_t5(L2a))./5;
    data2.y_value_vel = (data2.y_value_vel_t1(L2a) + data2.y_value_vel_t2(L2a) + data2.y_value_vel_t3(L2a) + data2.y_value_vel_t4(L2a) + data2.y_value_vel_t5(L2a))./5;
    data2.z_value_vel = (data2.z_value_vel_t1(L2a) + data2.z_value_vel_t2(L2a) + data2.z_value_vel_t3(L2a) + data2.z_value_vel_t4(L2a) + data2.z_value_vel_t5(L2a))./5;
    % WSS averaged over 5 systolic time frames
    data2.x_value_wss_t1 = data2.WSS_matrix{time-2}(:,1);data2.y_value_wss_t1 = data2.WSS_matrix{time-2}(:,2);data2.z_value_wss_t1 = data2.WSS_matrix{time-2}(:,3);
    data2.x_value_wss_t2 = data2.WSS_matrix{time-1}(:,1);data2.y_value_wss_t2 = data2.WSS_matrix{time-1}(:,2);data2.z_value_wss_t2 = data2.WSS_matrix{time-1}(:,3);
    data2.x_value_wss_t3 = data2.WSS_matrix{time}(:,1);  data2.y_value_wss_t3 = data2.WSS_matrix{time}(:,2);  data2.z_value_wss_t3 = data2.WSS_matrix{time}(:,3);
    data2.x_value_wss_t4 = data2.WSS_matrix{time+1}(:,1);data2.y_value_wss_t4 = data2.WSS_matrix{time+1}(:,2);data2.z_value_wss_t4 = data2.WSS_matrix{time+1}(:,3);
    data2.x_value_wss_t5 = data2.WSS_matrix{time+2}(:,1);data2.y_value_wss_t5 = data2.WSS_matrix{time+2}(:,2);data2.z_value_wss_t5 = data2.WSS_matrix{time+2}(:,3);
    data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4 + data2.x_value_wss_t5)./5;
    data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4 + data2.y_value_wss_t5)./5;
    data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4 + data2.z_value_wss_t5)./5;
end

% Velocity and WSS magnitude (scalar)
data2.vel_m = sqrt(data2.x_value_vel.^2 + data2.y_value_vel.^2 + data2.z_value_vel.^2);
data2.wss_m = sqrt(data2.x_value_wss.^2 + data2.y_value_wss.^2 + data2.z_value_wss.^2);

if plotFlag == 1
    figure('Name','data2 velocity')
    patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss], ...
        'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.1);
    hold on
    scatter3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,data2.vel_m.*50,data2.vel_m,'filled')
    colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    
    figure('Name','data2 WSS')
    patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss], ...
        'EdgeColor','none','FaceVertexCData',data2.wss_m,'FaceColor','interp','FaceAlpha',1);
    colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    
    figure('Name','To be registered')
    plot3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,'r.')
    hold on
    plot3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,'b.')
    legend('to remain the same','to be transformed')
    axis equal; axis ij;% axis off; view([-180 -90])
end

%%% Velocities
PSF = fspecial('gaussian',5,1);
mask1_to_register = mask1;
mask1_to_register = imfilter(mask1_to_register,PSF,'conv');

mask2 = imfilter(mask2,PSF,'conv');

% figure('Name','mask1')
% V = vol3d('cdata',mask1_to_register,'texture','3D','texturemap',mask1);
% axis tight; axis equal
%
% figure('Name','mask2')
% V = vol3d('cdata',mask2,'texture','3D','texturemap',mask2);
% axis tight; axis equal

disp(' ')
disp('affine registration (dof = 12)! So only scalars used (Affine registration doesnt work for vectors)...Busy registering...')
disp(' ')
tic
disp('...This can take up to 5 minutes...')

% directory with flirt.exe and cygwin1.dll (use cygwin convention)
fsldir = '/cygdrive/d/research/matlabcode/matlab_registration/flirt/';

% save as nifti
cnii=make_nii(mask1_to_register,[atlas.vox(1) atlas.vox(2) atlas.vox(3)]);
save_nii(cnii,'mask1.nii');
mnii=make_nii(mask2,[data2.vox(1) data2.vox(2) data2.vox(3)]);
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
system('c:\cygwin64\bin\bash runflirt.sh');

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
    figure('Name','Registered')
    plot3(atlas.x_coor_vel,atlas.y_coor_vel,atlas.z_coor_vel,'r.')
    hold on
    plot3(data2.x_coor_vel_new,data2.y_coor_vel_new,data2.z_coor_vel_new,'b.')
    legend('remains the same','transformed')
    axis equal; axis ij; axis off; view([-180 -90])
end

if calculateRE_Flag == 1
    x_vel_round = round(data2.x_coor_vel_new./atlas.vox(1));
    y_vel_round = round(data2.y_coor_vel_new./atlas.vox(2));
    z_vel_round = round(data2.z_coor_vel_new./atlas.vox(3));
    
    indices_mask2 = [x_vel_round y_vel_round z_vel_round];
    siz=max(indices_mask2,[],1);
    IND = sub2ind(siz,x_vel_round,y_vel_round,z_vel_round);
    [b, IND_double_removed, n] = unique(IND);
    clear b, clear n
    indices_mask2 = [x_vel_round(IND_double_removed) y_vel_round(IND_double_removed) z_vel_round(IND_double_removed)];
    clear IND_double_removed, clear IND
    
    mask2_new = zeros([max(y_vel_round) max(x_vel_round) max(z_vel_round)]);
    
    for i = 1:size(indices_mask2,1)
        mask2_new(indices_mask2(i,2),indices_mask2(i,1),indices_mask2(i,3)) = 1;
    end
    
    % Due to the rounding of coordinates there are holes in the aorta, we fill them by smooth3
    % and then erode the aorta to give it the right size again
    se = strel('disk',1);
    mask2_new = imerode(smooth3(mask2_new),se);
    
    % Make sure both masks have the same dimensions
    if size(mask1,1) > size(mask2_new,1)
        mask2_new(size(mask2_new,1):size(mask1,1),:,:) = 0;
    elseif size(mask1,1) < size(mask2_new,1)
        mask1(size(mask1,1):size(mask2_new,1),:,:) = 0;
    end
    if size(mask1,2) > size(mask2_new,2)
        mask2_new(:,size(mask2_new,2):size(mask1,2),:) = 0;
    elseif size(mask1,2) < size(mask2_new,2)
        mask1(:,size(mask1,2):size(mask2_new,2),:) = 0;
    end
    if size(mask1,3) > size(mask2_new,3)
        mask2_new(:,:,size(mask2_new,3):size(mask1,3)) = 0;
    elseif size(mask1,3) < size(mask2_new,3)
        mask1(:,:,size(mask1,3):size(mask2_new,3)) = 0;
    end
    
    L_mask2 = double(mask2_new ~= 0);
    
    difference = abs(mask1-L_mask2);
    [I1,J] = find(mask1~=0);
    [I2,J] = find(L_mask2~=0);
    mean_I = (size(I1,1) + size(I2,1))/2;
    [I_diff,J] = find(difference~=0);
    diff_voxels = size(I_diff,1);
    diff_percentage = ((diff_voxels / mean_I) * 100)/2;
    disp(['RE: Difference between aorta and atlas = ' num2str(diff_percentage)])
    disp(' ')
end

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
    
    if ~exist(strcat(PATHNAME,'atlas_interpolation_error_ROI_after_transformation\mask1.mat'),'file')
        
        F2=figure('Name','Atlas shape: Paused after finishing a region so press space when finished!');
        plot3(data2.x_coor_wss_new,data2.y_coor_wss_new,data2.z_coor_wss_new,'r.');
        %  patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss_new data2.y_coor_wss_new data2.z_coor_wss_new],'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.25);
        view([-180 -90]);axis ij;axis equal;axis off
        
        mkdir(PATHNAME,'atlas_interpolation_error_ROI_after_transformation')
        
        for i = 1:6
            %Polygon and mask for AAo
            polyAAo = impoly;
            wait(polyAAo);
            region = getPosition(polyAAo);
            
            %          disp('saving, pausing')
            save(strcat([PATHNAME 'atlas_interpolation_error_ROI_after_transformation\mask' num2str(i)]),'region');
            pause
        end
        
        close(F2)
    end
    load(strcat(PATHNAME,'atlas_interpolation_error_ROI_after_transformation\mask1'))
    transformed_atlas_mask_AAo_inner_vel = inpolygon(data2.x_coor_vel_new,data2.y_coor_vel_new, region(:,1), region(:,2));
    transformed_mean_vel_asc_inner = mean(atlas_mean_vel(transformed_atlas_mask_AAo_inner_vel));
    transformed_atlas_mask_AAo_inner_wss = inpolygon(data2.x_coor_wss_new,data2.y_coor_wss_new, region(:,1), region(:,2));
    transformed_mean_wss_asc_inner = mean(atlas_mean_wss(transformed_atlas_mask_AAo_inner_wss));
    load(strcat(PATHNAME,'atlas_interpolation_error_ROI_after_transformation\mask2'))
    transformed_atlas_mask_AAo_outer_vel = inpolygon(data2.x_coor_vel_new,data2.y_coor_vel_new, region(:,1), region(:,2));
    transformed_mean_vel_asc_outer = mean(atlas_mean_vel(transformed_atlas_mask_AAo_outer_vel));
    transformed_atlas_mask_AAo_outer_wss = inpolygon(data2.x_coor_wss_new,data2.y_coor_wss_new, region(:,1), region(:,2));
    transformed_mean_wss_asc_outer = mean(atlas_mean_wss(transformed_atlas_mask_AAo_outer_wss));
    load(strcat(PATHNAME,'atlas_interpolation_error_ROI_after_transformation\mask3'))
    transformed_atlas_mask_arch_inner_vel = inpolygon(data2.x_coor_vel_new,data2.y_coor_vel_new, region(:,1), region(:,2));
    transformed_mean_vel_arch_inner = mean(atlas_mean_vel(transformed_atlas_mask_arch_inner_vel));
    transformed_atlas_mask_arch_inner_wss = inpolygon(data2.x_coor_wss_new,data2.y_coor_wss_new, region(:,1), region(:,2));
    transformed_mean_wss_arch_inner = mean(atlas_mean_wss(transformed_atlas_mask_arch_inner_wss));
    load(strcat(PATHNAME,'atlas_interpolation_error_ROI_after_transformation\mask4'))
    transformed_atlas_mask_arch_outer_vel = inpolygon(data2.x_coor_vel_new,data2.y_coor_vel_new, region(:,1), region(:,2));
    transformed_mean_vel_arch_outer = mean(atlas_mean_vel(transformed_atlas_mask_arch_outer_vel));
    transformed_atlas_mask_arch_outer_wss = inpolygon(data2.x_coor_wss_new,data2.y_coor_wss_new, region(:,1), region(:,2));
    transformed_mean_wss_arch_outer = mean(atlas_mean_wss(transformed_atlas_mask_arch_outer_wss));
    load(strcat(PATHNAME,'atlas_interpolation_error_ROI_after_transformation\mask5'))
    transformed_atlas_mask_DAo_inner_vel = inpolygon(data2.x_coor_vel_new,data2.y_coor_vel_new, region(:,1), region(:,2));
    transformed_mean_vel_DAo_inner = mean(atlas_mean_vel(transformed_atlas_mask_DAo_inner_vel));
    transformed_atlas_mask_DAo_inner_wss = inpolygon(data2.x_coor_wss_new,data2.y_coor_wss_new, region(:,1), region(:,2));
    transformed_mean_wss_DAo_inner = mean(atlas_mean_wss(transformed_atlas_mask_DAo_inner_wss));
    load(strcat(PATHNAME,'atlas_interpolation_error_ROI_after_transformation\mask6'))
    transformed_atlas_mask_DAo_outer_vel = inpolygon(data2.x_coor_vel_new,data2.y_coor_vel_new, region(:,1), region(:,2));
    transformed_mean_vel_DAo_outer = mean(atlas_mean_vel(transformed_atlas_mask_DAo_outer_vel));
    transformed_atlas_mask_DAo_outer_wss = inpolygon(data2.x_coor_wss_new,data2.y_coor_wss_new, region(:,1), region(:,2));
    transformed_mean_wss_DAo_outer = mean(atlas_mean_wss(transformed_atlas_mask_DAo_outer_wss));
    
    mean_vel_after_interpolation(1,1) = transformed_mean_vel_asc_inner;
    mean_vel_after_interpolation(2,1) = transformed_mean_vel_asc_outer;
    mean_vel_after_interpolation(3,1) = transformed_mean_vel_arch_inner;
    mean_vel_after_interpolation(4,1) = transformed_mean_vel_arch_outer;
    mean_vel_after_interpolation(5,1) = transformed_mean_vel_DAo_inner;
    mean_vel_after_interpolation(6,1) = transformed_mean_vel_DAo_outer;
    mean_wss_after_interpolation(1,1) = transformed_mean_wss_asc_inner;
    mean_wss_after_interpolation(2,1) = transformed_mean_wss_asc_outer;
    mean_wss_after_interpolation(3,1) = transformed_mean_wss_arch_inner;
    mean_wss_after_interpolation(4,1) = transformed_mean_wss_arch_outer;
    mean_wss_after_interpolation(5,1) = transformed_mean_wss_DAo_inner;
    mean_wss_after_interpolation(6,1) = transformed_mean_wss_DAo_outer;
    
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
        figure('Name','Velocity after interpolation: outer arch')
        scatter3(data2.x_coor_vel_new(transformed_atlas_mask_arch_outer_vel),data2.y_coor_vel_new(transformed_atlas_mask_arch_outer_vel),data2.z_coor_vel_new(transformed_atlas_mask_arch_outer_vel),20,atlas_mean_vel(transformed_atlas_mask_arch_outer_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','Velocity after interpolation: inner DAo')
        scatter3(data2.x_coor_vel_new(transformed_atlas_mask_DAo_inner_vel),data2.y_coor_vel_new(transformed_atlas_mask_DAo_inner_vel),data2.z_coor_vel_new(transformed_atlas_mask_DAo_inner_vel),20,atlas_mean_vel(transformed_atlas_mask_DAo_inner_vel),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','Velocity after interpolation: outer DAo')
        scatter3(data2.x_coor_vel_new(transformed_atlas_mask_DAo_outer_vel),data2.y_coor_vel_new(transformed_atlas_mask_DAo_outer_vel),data2.z_coor_vel_new(transformed_atlas_mask_DAo_outer_vel),20,atlas_mean_vel(transformed_atlas_mask_DAo_outer_vel),'filled');axis equal;caxis([0 1.5])
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
        figure('Name','WSS after interpolation: outer arch')
        scatter3(data2.x_coor_wss_new(transformed_atlas_mask_arch_outer_wss),data2.y_coor_wss_new(transformed_atlas_mask_arch_outer_wss),data2.z_coor_wss_new(transformed_atlas_mask_arch_outer_wss),20,atlas_mean_wss(transformed_atlas_mask_arch_outer_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS after interpolation: inner DAo')
        scatter3(data2.x_coor_wss_new(transformed_atlas_mask_DAo_inner_wss),data2.y_coor_wss_new(transformed_atlas_mask_DAo_inner_wss),data2.z_coor_wss_new(transformed_atlas_mask_DAo_inner_wss),20,atlas_mean_wss(transformed_atlas_mask_DAo_inner_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure('Name','WSS after interpolation: outer DAo')
        scatter3(data2.x_coor_wss_new(transformed_atlas_mask_DAo_outer_wss),data2.y_coor_wss_new(transformed_atlas_mask_DAo_outer_wss),data2.z_coor_wss_new(transformed_atlas_mask_DAo_outer_wss),20,atlas_mean_wss(transformed_atlas_mask_DAo_outer_wss),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
    end
    
    IE_inner_AAo_vel = abs(mean_vel_before_interpolation(1,1)-mean_vel_after_interpolation(1,1)) / ((mean_vel_before_interpolation(1,1)+mean_vel_after_interpolation(1,1))./2)*100;
    IE_outer_AAo_vel = abs(mean_vel_before_interpolation(2,1)-mean_vel_after_interpolation(2,1)) / ((mean_vel_before_interpolation(2,1)+mean_vel_after_interpolation(2,1))./2)*100;
    IE_inner_asc_vel = abs(mean_vel_before_interpolation(3,1)-mean_vel_after_interpolation(3,1)) / ((mean_vel_before_interpolation(3,1)+mean_vel_after_interpolation(3,1))./2)*100;
    IE_outer_asc_vel = abs(mean_vel_before_interpolation(4,1)-mean_vel_after_interpolation(4,1)) / ((mean_vel_before_interpolation(4,1)+mean_vel_after_interpolation(4,1))./2)*100;
    IE_inner_DAo_vel = abs(mean_vel_before_interpolation(5,1)-mean_vel_after_interpolation(5,1)) / ((mean_vel_before_interpolation(5,1)+mean_vel_after_interpolation(5,1))./2)*100;
    IE_outer_DAo_vel = abs(mean_vel_before_interpolation(6,1)-mean_vel_after_interpolation(6,1)) / ((mean_vel_before_interpolation(6,1)+mean_vel_after_interpolation(6,1))./2)*100;
    IE_inner_AAo_wss = abs(mean_wss_before_interpolation(1,1)-mean_wss_after_interpolation(1,1)) / ((mean_wss_before_interpolation(1,1)+mean_wss_after_interpolation(1,1))./2)*100;
    IE_outer_AAo_wss = abs(mean_wss_before_interpolation(2,1)-mean_wss_after_interpolation(2,1)) / ((mean_wss_before_interpolation(2,1)+mean_wss_after_interpolation(2,1))./2)*100;
    IE_inner_asc_wss = abs(mean_wss_before_interpolation(3,1)-mean_wss_after_interpolation(3,1)) / ((mean_wss_before_interpolation(3,1)+mean_wss_after_interpolation(3,1))./2)*100;
    IE_outer_asc_wss = abs(mean_wss_before_interpolation(4,1)-mean_wss_after_interpolation(4,1)) / ((mean_wss_before_interpolation(4,1)+mean_wss_after_interpolation(4,1))./2)*100;
    IE_inner_DAo_wss = abs(mean_wss_before_interpolation(5,1)-mean_wss_after_interpolation(5,1)) / ((mean_wss_before_interpolation(5,1)+mean_wss_after_interpolation(5,1))./2)*100;
    IE_outer_DAo_wss = abs(mean_wss_before_interpolation(6,1)-mean_wss_after_interpolation(6,1)) / ((mean_wss_before_interpolation(6,1)+mean_wss_after_interpolation(6,1))./2)*100;
    
    disp(['IE velocity inner AAo = ' num2str(IE_inner_AAo_vel) ' %'])
    disp(['IE velocity outer AAo = ' num2str(IE_outer_AAo_vel) ' %'])
    disp(['IE velocity inner asc = ' num2str(IE_inner_asc_vel) ' %'])
    disp(['IE velocity outer asc = ' num2str(IE_outer_asc_vel) ' %'])
    disp(['IE velocity inner DAo = ' num2str(IE_inner_DAo_vel) ' %'])
    disp(['IE velocity outer DAo = ' num2str(IE_outer_DAo_vel) ' %'])
    disp(' ')
    disp(['IE wall shear stress inner AAo = ' num2str(IE_inner_AAo_wss) ' %'])
    disp(['IE wall shear stress outer AAo = ' num2str(IE_outer_AAo_wss) ' %'])
    disp(['IE wall shear stress inner asc = ' num2str(IE_inner_asc_wss) ' %'])
    disp(['IE wall shear stress outer asc = ' num2str(IE_outer_asc_wss) ' %'])
    disp(['IE wall shear stress inner DAo = ' num2str(IE_inner_DAo_wss) ' %'])
    disp(['IE wall shear stress outer DAo = ' num2str(IE_outer_DAo_wss) ' %'])
    disp(' ')
end

if plotFlag == 1
    figure('Name','Mean atlas velocity')
    scatter3(data2.x_coor_vel_new,data2.y_coor_vel_new,data2.z_coor_vel_new,atlas_mean_vel.*50,atlas_mean_vel,'filled')
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90]);
    
    figure('Name','std atlas velocity')
    scatter3(data2.x_coor_vel_new,data2.y_coor_vel_new,data2.z_coor_vel_new,atlas_std_vel.*50,atlas_std_vel,'filled')
    axis equal;axis off; axis ij;
    caxis([0 1.5])
    view([180 -90])
    
    figure('Name','Mean atlas WSS')
    patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss_new data2.y_coor_wss_new data2.z_coor_wss_new],'EdgeColor','none', 'FaceVertexCData',atlas_mean_wss,'FaceColor','interp','FaceAlpha',1);colorbar;
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
    
    figure('Name','std atlas WSS')
    patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss_new data2.y_coor_wss_new data2.z_coor_wss_new],'EdgeColor','none', 'FaceVertexCData',atlas_std_wss,'FaceColor','interp','FaceAlpha',1);colorbar;
    axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90])
end

% Determine thresholds
mean_plus_2SD_atlas_vel = atlas_mean_vel + 1.96.*atlas_std_vel;
mean_plus_1SD_atlas_vel = atlas_mean_vel + 0.98.*atlas_std_vel;
mean_plus_2SD_atlas_wss = atlas_mean_wss + 1.96.*atlas_std_wss;
mean_min_2SD_atlas_wss = atlas_mean_wss - 1.96.*atlas_std_wss;

if plotFlag == 1
    % Velocity
    figure('Name','mean_plus_2SD_atlas velocity')
    scatter3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,atlas_mean_vel.*50,mean_plus_2SD_atlas_vel,'filled')
    colorbar;axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90]);
    
    figure('Name','mean_plus_1SD_atlas velocity')
    scatter3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,atlas_mean_vel.*50,mean_plus_1SD_atlas_vel,'filled')
    colorbar;axis equal;axis off; axis ij;caxis([0 1.5]);view([180 -90]);
    
    % WSS
    figure('Name','mean_plus_2SD_atlas WSS')
    patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',mean_plus_2SD_atlas_wss,'FaceColor','interp','FaceAlpha',1);
    colorbar;axis equal;axis off; axis ij;caxis([0 3]);view([180 -90]);
    
    figure('Name','mean_min_2SD_atlas WSS')
    patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',mean_min_2SD_atlas_wss,'FaceColor','interp','FaceAlpha',1);
    colorbar;axis equal;axis off; axis ij;caxis([0 3]);view([180 -90]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TRAFFIC LIGHT MAP FOR ABNORMAL VELOCITY %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_mask_red_ = zeros(size(data2.vel_m));
new_mask_yellow_ = zeros(size(data2.vel_m));
new_mask_green_ = zeros(size(data2.vel_m));
for i=1:size(data2.vel_m,1)
    if data2.vel_m(i) > mean_plus_2SD_atlas_vel(i)
        new_mask_red_(i,1) = 1;
        new_mask_yellow_(i,1) = -1;
        new_mask_green_(i,1) = -1;
    elseif data2.vel_m(i) > mean_plus_1SD_atlas_vel(i)
        new_mask_red_(i,1) = -1;
        new_mask_yellow_(i,1) = 1;
        new_mask_green_(i,1) = -1;
    else
        new_mask_red_(i,1) = -1;
        new_mask_yellow_(i,1) = -1;
        new_mask_green_(i,1) = 1;
    end
end

% translate point cloud into matrix
new_mask_red = zeros(size(L2b));
new_mask_red(L2b) = new_mask_red_;
new_mask_red(~L2b) = -1;
new_mask_yellow = zeros(size(L2b));
new_mask_yellow(L2b) = new_mask_yellow_;
new_mask_yellow(~L2b) = -1;
new_mask_green = zeros(size(L2b));
new_mask_green(L2b) = new_mask_green_;
new_mask_green(~L2b) = -1;

if calculate_traffic_light_volumeFlag == 1
    
    [I,J] = find(L2b~=0);
    total_volume = data2.vox(1)*data2.vox(2)*data2.vox(3)*size(I,1);
    
    [I,J] = find(new_mask_red==1);
    red_volume = data2.vox(1)*data2.vox(2)*data2.vox(3)*size(I,1);
    percentage_red_volume = red_volume / total_volume * 100;
    [I,J] = find(new_mask_yellow==1);
    yellow_volume = data2.vox(1)*data2.vox(2)*data2.vox(3)*size(I,1);
    percentage_yellow_volume = yellow_volume / total_volume * 100;
    [I,J] = find(new_mask_green==1);
    green_volume = data2.vox(1)*data2.vox(2)*data2.vox(3)*size(I,1);
    percentage_green_volume = green_volume / total_volume * 100;
    total_percentage = percentage_red_volume + percentage_yellow_volume + percentage_green_volume;
    
    disp(['Red volume percentage of total aorta = ' num2str(percentage_red_volume) ' %'])
    disp(['Yellow volume percentage of total aorta = ' num2str(percentage_yellow_volume) ' %'])
    disp(['Green volume percentage of total aorta = ' num2str(percentage_green_volume) ' %'])
    disp(['Total percentage of total aorta = ' num2str(total_percentage) ' %'])
    disp(' ')
end

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

if calculate_area_of_higherlowerFlag == 1;
    if ~exist(strcat(PATHNAME,'heat_map_higher_lower_masks\mask1.mat'),'file')      
      
        patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',1);        
        view([-180 -90]);axis ij;axis equal;axis off   
        
        mkdir(PATHNAME,'heat_map_higher_lower_masks')
        
        for i = 1:6
            %Polygon and mask for AAo
            polyAAo = impoly;
            wait(polyAAo);
            region = getPosition(polyAAo);
            %mask = inpolygon(x, y, AAo(:,1), AAo(:,2));
            
            disp('saving, pausing')
            save(strcat([PATHNAME 'heat_map_higher_lower_masks\mask' num2str(i)]),'region');
            pause
        end
    end
    load(strcat(PATHNAME,'heat_map_higher_lower_masks\mask1'));
    atlas_mask_AAo_inner = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
    load(strcat(PATHNAME,'heat_map_higher_lower_masks\mask2'));
    atlas_mask_AAo_outer = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
    load(strcat(PATHNAME,'heat_map_higher_lower_masks\mask3'));
    atlas_mask_arch_inner = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
    load(strcat(PATHNAME,'heat_map_higher_lower_masks\mask4'));
    atlas_mask_arch_outer = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
    load(strcat(PATHNAME,'heat_map_higher_lower_masks\mask5'));
    atlas_mask_DAo_inner = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
    load(strcat(PATHNAME,'heat_map_higher_lower_masks\mask6'));
    atlas_mask_DAo_outer = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
    
    heat_asc1 = heat_mapp(atlas_mask_AAo_inner);
    [I1,J1] = find(heat_asc1 == 0);
    [I2,J2] = find(heat_asc1 == 1);
    percentage_significant_higher_than_controls = size(I2,1) / size(heat_asc1,1) * 100;
    percentage_significant_lower_than_controls = size(I1,1) / size(heat_asc1,1) * 100;
    
    disp(['Percentage higher than controls inner AAo = ' num2str(percentage_significant_higher_than_controls) '%'])
    disp(['Percentage lower than controls inner AAo = ' num2str(percentage_significant_lower_than_controls) '%'])
    
    heat_asc2 = heat_mapp(atlas_mask_AAo_outer);
    [I1,J1] = find(heat_asc2 == 0);
    [I2,J2] = find(heat_asc2 == 1);
    percentage_significant_higher_than_controls = size(I2,1) / size(heat_asc2,1) * 100;
    percentage_significant_lower_than_controls = size(I1,1) / size(heat_asc2,1) * 100;
    
    disp(['Percentage higher than controls outer AAo = ' num2str(percentage_significant_higher_than_controls) '%'])
    disp(['Percentage lower than controls outer AAo = ' num2str(percentage_significant_lower_than_controls) '%'])
    
    heat_arch1 = heat_mapp(atlas_mask_arch_inner);
    [I1,J1] = find(heat_arch1 == 0);
    [I2,J2] = find(heat_arch1 == 1);
    percentage_significant_higher_than_controls = size(I2,1) / size(heat_arch1,1) * 100;
    percentage_significant_lower_than_controls = size(I1,1) / size(heat_arch1,1) * 100;
    
    disp(['Percentage higher than controls inner arch = ' num2str(percentage_significant_higher_than_controls) '%'])
    disp(['Percentage lower than controls inner arch = ' num2str(percentage_significant_lower_than_controls) '%'])
    
    heat_arch2 = heat_mapp(atlas_mask_arch_outer);
    [I1,J1] = find(heat_arch2 == 0);
    [I2,J2] = find(heat_arch2 == 1);
    percentage_significant_higher_than_controls = size(I2,1) / size(heat_arch2,1) * 100;
    percentage_significant_lower_than_controls = size(I1,1) / size(heat_arch2,1) * 100;
    
    disp(['Percentage higher than controls outer arch = ' num2str(percentage_significant_higher_than_controls) '%'])
    disp(['Percentage lower than controls outer arch = ' num2str(percentage_significant_lower_than_controls) '%'])
    
    heat_desc1 = heat_mapp(atlas_mask_DAo_inner);
    [I1,J1] = find(heat_desc1 == 0);
    [I2,J2] = find(heat_desc1 == 1);
    percentage_significant_higher_than_controls = size(I2,1) / size(heat_desc1,1) * 100;
    percentage_significant_lower_than_controls = size(I1,1) / size(heat_desc1,1) * 100;
    
    disp(['Percentage higher than controls inner DAo = ' num2str(percentage_significant_higher_than_controls) '%'])
    disp(['Percentage lower than controls inner DAo = ' num2str(percentage_significant_lower_than_controls) '%'])
    
    heat_desc2 = heat_mapp(atlas_mask_DAo_outer);
    [I1,J1] = find(heat_desc2 == 0);
    [I2,J2] = find(heat_desc2 == 1);
    percentage_significant_higher_than_controls = size(I2,1) / size(heat_desc2,1) * 100;
    percentage_significant_lower_than_controls = size(I1,1) / size(heat_desc2,1) * 100;
    
    disp(['Percentage higher than controls outer DAo = ' num2str(percentage_significant_higher_than_controls) '%'])
    disp(['Percentage lower than controls outer DAo = ' num2str(percentage_significant_lower_than_controls) '%'])
    
    if plotFlag == 1
    figure('Name','higher/lower: inner AAo')
    scatter3(data2.x_coor_wss(atlas_mask_AAo_inner),data2.y_coor_wss(atlas_mask_AAo_inner),data2.z_coor_wss(atlas_mask_AAo_inner),20,heat_asc1,'filled');axis equal;colormap(color1);colorbar;caxis([0 2])
    xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
    ffigure('Name','higher/lower: outer AAo')
    scatter3(data2.x_coor_wss(atlas_mask_AAo_outer),data2.y_coor_wss(atlas_mask_AAo_outer),data2.z_coor_wss(atlas_mask_AAo_outer),20,heat_asc2,'filled');axis equal;colormap(color1);colorbar;caxis([0 2])
    xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
    figure('Name','higher/lower: inner arch')
    scatter3(data2.x_coor_wss(atlas_mask_arch_inner),data2.y_coor_wss(atlas_mask_arch_inner),data2.z_coor_wss(atlas_mask_arch_inner),20,heat_arch1,'filled');axis equal;colormap(color1);colorbar;caxis([0 2])
    xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
    figure('Name','higher/lower: outer arch')
    scatter3(data2.x_coor_wss(atlas_mask_arch_outer),data2.y_coor_wss(atlas_mask_arch_outer),data2.z_coor_wss(atlas_mask_arch_outer),20,heat_arch2,'filled');axis equal;colormap(color1);colorbar;caxis([0 2])
    xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
    figure('Name','higher/lower: inner DAo')
    scatter3(data2.x_coor_wss(atlas_mask_DAo_inner),data2.y_coor_wss(atlas_mask_DAo_inner),data2.z_coor_wss(atlas_mask_DAo_inner),20,heat_desc1,'filled');axis equal;colormap(color1);colorbar;caxis([0 2])
    xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
    figure('Name','higher/lower: outer DAo')
    scatter3(data2.x_coor_wss(atlas_mask_DAo_outer),data2.y_coor_wss(atlas_mask_DAo_outer),data2.z_coor_wss(atlas_mask_DAo_outer),20,heat_desc2,'filled');axis equal;colormap(color1);colorbar;caxis([0 2])
    xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  TRAFFIC LIGHT MAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count1 = 0;
angles(1) = 0;
f1 = figure('Name','Traffic Light Map');
contours = zeros(size(L2b));
contours(L2b==0) = -1;
contours(L2b==1) = 1;
[F,V] = isosurface(smooth3(contours),0); % make a surface from the detected contours
patch('Faces',F,'Vertices',[V(:,1)-offset V(:,2)-offset V(:,3)-offset],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.1);
hold on
[F1,V1] = isosurface(x./data2.vox(1)-offset,y./data2.vox(2)-offset,z./data2.vox(3)-offset,smooth3(new_mask_red),0);
[F2,V2] = isosurface(x./data2.vox(1)-offset,y./data2.vox(2)-offset,z./data2.vox(3)-offset,smooth3(new_mask_yellow),0);
[F3,V3] = isosurface(x./data2.vox(1)-offset,y./data2.vox(2)-offset,z./data2.vox(3)-offset,smooth3(new_mask_green),0);
p11=patch('Faces',F1,'Vertices',V1,'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',1);
p12=patch('Faces',F2,'Vertices',V2,'EdgeColor','none','FaceColor',[1 0.9 0],'FaceAlpha',1);
p13=patch('Faces',F3,'Vertices',V3,'EdgeColor','none','FaceColor',[0 1 0],'FaceAlpha',1);
axis equal; axis off;axis ij
magnitude = flipdim(double(data2.mag_struct.dataAy(:,:,:,time)),3);
magnitude(magnitude == 0) = 3;
magnitude(magnitude == 1) = 3;
magnitude(magnitude == 2) = 3;
hold on
s1 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3),magnitude(:,:,1),'EdgeColor','none');
set(s1,'HandleVisibility','off','Visible','off');
axis equal;colormap(gray)
view([-180 -90])
caxis([0 64]);
aspectRatio = 1./data2.vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))
camlight headlight;camlight(180,0); lighting phong

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PvO: Images for the Brief Report for NEJM were created with:
%%% 1. NO axis ij: I couldn't get the lighting right with axis ij on
%%%    The images therefore had to be flipped horizontly for the paper.
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
    'Value',1, 'Position', [10 250 20 20], ...
    'Callback', {@show_yellow_region,gca});

uicontrol('Style','text',...
    'Position',[15 225 120 20],...
    'String','Show Green Map')
uicontrol('Style','checkbox',...
    'Value',1, 'Position', [10 225 20 20], ...
    'Callback', {@show_green_region,gca});

uicontrol('Style','text',...
    'Position',[15 200 120 20],...
    'String','Show Anatomy')
uicontrol('Style','checkbox',...
    'Value',0, 'Position', [10 200 20 20], ...
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

    function show_green_region(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
            patchobj = findobj(p13);
            set(patchobj,'HandleVisibility','on','Visible','on');
        elseif show == 0
            patchobj = findobj(p13);
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

traffic_light.red = new_mask_red;
traffic_light.yellow = new_mask_yellow;
traffic_light.green = new_mask_green;
traffic_light.vertices = [V(:,1)-offset V(:,2)-offset V(:,3)-offset];
traffic_light.faces = F;

% set up results folder
dir_orig = pwd;
dir_new = PATHNAME; cd(dir_new); %cd('..')
%dir_new = pwd;
mkdir('results_traffic_light_map')
dir_new
dir_new = strcat(dir_new,'results_traffic_light_map');

% save results in results folder
save(strcat(dir_new,'\results_trafficlight_map'),'traffic_light');
%savefig(f2,strcat(dir_new,'\heat_map'))
print(f1,'-zbuffer', '-djpeg','-r600',strcat(dir_new,'\traffic_light_map.jpg'));
cd(dir_orig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Heat map (WSS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count3 = 0;
angles(1) = 0;
f2 = figure('Name','Heat map');
x = data2.x_coor_wss./data2.vox(1) - offset;
y = data2.y_coor_wss./data2.vox(2) - offset;
z = data2.z_coor_wss./data2.vox(3) - offset;
p2=patch('Faces',data2.F_matrix{1},'Vertices',[x y z],'EdgeColor','none', 'FaceVertexCData',heat_mapp,'FaceColor','interp','FaceAlpha',1);
gray_colormap = colormap(gray);
color2(1,:) = [0 0 1];
color2(2,:) = [1 0 0];
color2(3,:) = [0.5 0.5 0.5];
color2(4:64,:) = gray_colormap(4:64,:);
colormap(color2);
caxis([0 64]);
axis equal; axis ij; axis off;
magnitude = flipdim(double(data2.mag_struct.dataAy(:,:,:,time)),3);
magnitude(magnitude == 0) = 3;
magnitude(magnitude == 1) = 3;
magnitude(magnitude == 2) = 3;
hold on
s2 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3),magnitude(:,:,1),'EdgeColor','none');
set(s2,'HandleVisibility','off','Visible','off');
axis equal;
view([-180 -90])
aspectRatio = 1./data2.vox;
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
 
uicontrol('Style','text',...
    'Position',[15 300 120 20],...
    'String','Show Heat Map')
uicontrol('Style','checkbox',...
    'Value',1, 'Position', [15 300 20 20], ...
    'Callback', {@show_heat_mapp,gca});

    function move_slice2(hObj,event,ax)
        sliceobj = findobj(s2);
        delete(sliceobj)
        slice = round(get(hObj,'Value'));
        s2 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]).*(size(magnitude,3)-(slice-1)),magnitude(:,:,slice),'EdgeColor','none');
    end

    function show_heat_mapp(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
            patchobj = findobj(p2);
            set(patchobj,'HandleVisibility','on','Visible','on');
        elseif show == 0
            patchobj = findobj(p2);
            set(patchobj,'HandleVisibility','off','Visible','off');
        end
    end

    function show_anatomy2(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
            patchobj = findobj(s2);
            set(patchobj,'HandleVisibility','on','Visible','on');
        elseif show == 0
            patchobj = findobj(s2);
            set(patchobj,'HandleVisibility','off','Visible','off');
        end
    end

    function change_contrast2(hObj,event,ax)
        contrast = round(get(hObj,'Value'));
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

heat_map.heat_map = heat_mapp;
heat_map.vertices = [x y z];
heat_map.faces = data2.F_matrix{1};
heat_map.color = color1;

% set up results folder
dir_orig = pwd;
dir_new = PATHNAME; cd(dir_new); %cd('..')
%dir_new = pwd;
mkdir('results_heatmap')
dir_new
dir_new = strcat(dir_new,'results_heatmap');

% save results in results folder
save(strcat(dir_new,'\heat_map'),'heat_map')
%savefig(f2,strcat(dir_new,'\heat_map'))
print(f2,'-zbuffer', '-djpeg','-r600',strcat(dir_new,'\heat_map.jpg'));
cd(dir_orig)
end