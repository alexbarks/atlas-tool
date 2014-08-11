function [p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(offset,plotFlag,calculateRE_Flag,calculateIE_Flag,calculate_p_value_volumeFlag,calculate_area_of_significance_wss,peak_systolicFlag)

%%% [p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(offset,plotFlag,calculateRE_Flag,calculateIE_Flag,calculate_p_value_volumeFlag,calculate_area_of_significance_wss,peak_systolicFlag)
%
% This function creates the p-value map for velocity and WSS atlas (ensemble-averaged velocity/WSS map is more politically correct)
% of the batch data put into this file (see below). The p-value map shows regions of significant differences between cohorts
% The method was published in (so please cite this paper when using this method):
%
% A Methodology to Detect Abnormal Relative Wall Shear Stress on the Full Surface of the Thoracic Aorta Using 4D Flow MRI
% van Ooij P, Potters WV, Nederveen AJ, Allen BD, Collins J, Carr J, Malaisrie SC, Markl M, Barker AJ
% Magn Reson Med. 2014 Feb 25, doi: 10.1002/mrm.25224. [Epub ahead of print]
%
% And:
%
% Characterization of abnormal wall shear stress using 4D flow MRI in human bicuspid aortopathy
% van Ooij P, Potters WV, Collins J, Carr M, Carr J, Malaisrie SC, Fedak PWM, McCarthy PM, Markl M, Barker AJ
% Ann Biomed Eng. Accepted 2014 August 8. IF 2014: 3.231.
%
%
% First, each aorta in both cohorts (you can probably use this for carotids and intracranial vessels as well, but I haven't tried yet) is co-registered
% (affine registration) to the idealized aorta geometry created with the function 'make_geometry_point_cloud.m'. The velocity and WSS magnitude
% values are subsequently interpolated (nearest neighbour interpolation) to the velocity and WSS coordinates of the idealized geometry. Now you have every
% aorta in each cohort on the same geometry. You can now perform a Wilcoxon rank sum test to test for significance between both cohorts on each position
% on the geometry. Differences for cohort 2 that are significantly higher than cohort 1 are visualized in red, differences for cohort 2 that are
% significantly lower than cohort 1 are visualized in blue.
%
% 2014, Pim van Ooij, Northwestern University
%
% Input
% 1)offset          : To be able to calculate the registration error (RE, see the first paper mentioned above), the registered aorta needs to be
%                     transformed back to a matrix. This can only be done when all x, y and z coordinates after registration are > 0.
%                     Otherwise Matlab will return an error and all will be for nothing. We therefore need to put in an offset to prevent
%                     this from happening.
% 2)plotFlag        : If plotFlag is switched on, Matlab will output any possible image to the screen to check if everything happens correctly.
% 3)calculateRE_Flag: When switched on the registration error will be calculated. However, this is a different RE than the one in the paper
%                     mentioned above as the RE in this function is calculated from AFFINE registration whereas the RE reported in the paper
%                     is calculated from RIGID registration in the function 'make_geometry_point_cloud.m'
% 4)calculateIE_Flag: When switched on the interpolation error (IE, see the first paper mentioned above) will be calculated. Note that ROIs are needed
%                     which can be drawn manually when switched on.
% 5)calculate_p_value_volumeFlag: calulate the volume of significant different velocity
% 5)peak_systolicFlag: When switched on the atlas will be created for the peak systolic time frame only, default is the atlass for 5 systolic
%                     timesteps averaged.
%
% Output
% 1)p_value_map:       The p_value_map is saved to this struct
%
% Usage
% This code is for creating p-value maps from 'data_done' structs as created by Pims_postprocessing tool
% The function for creating atlases from mrStructs is under construction
%
% Examples:
% [p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(offset,plotFlag,calculateRE_Flag,calculateIE_Flag,calculate_p_value_volumeFlag,calculate_area_of_significance_wss,peak_systolicFlag)
% [p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(40,1,1,1,1,1,1)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all % I hate functions that don't close the stuff that was running before

global data
Rotation_Translation = [];

if nargin < 1 || isempty(offset)
    offset = 40;
end

if nargin < 2 || isempty(plotFlag)
    plotFlag = 1;
end

if nargin < 3 || isempty(calculateRE_Flag)
    calculateRE_Flag = 0;
end

if nargin < 4 || isempty(calculateIE_Flag)
    calculateIE_Flag = 0;
end

if nargin < 5 || isempty(calculate_p_value_volumeFlag)
    calculate_p_value_volumeFlag = 1;
end

if nargin < 6 || isempty(calculate_area_of_significance_wss)
    calculate_area_of_significance_wss = 1;
end

if nargin < 7 || isempty(peak_systolicFlag)
    peak_systolicFlag = 0;
end

PATHNAME_probability_mask = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls';
probability_mask = 'probability_mask.mat';
load(strcat(PATHNAME_probability_mask,probability_mask));

%%%% COHORT 1: datasets to load
%%%% Copy paste the folders of choice here, add more or discard if needed
PATHNAME1{1} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\1_20120420_132106\';
PATHNAME1{2} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\2_20120426_132244\';
PATHNAME1{3} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\3_20121206_115454\';
PATHNAME1{4} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\4_20120522_170003\';
PATHNAME1{5} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\5_20120502_134311\';
PATHNAME1{6} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\6_20120627_093614\';
PATHNAME1{7} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\7_20120702_092347\';
PATHNAME1{8} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\8_20120831_101148\';
PATHNAME1{9} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\9_20121109_080007\';
PATHNAME1{10}= 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\10_20130621_122315\';
FILENAME = 'data_done';

%%%% COHORT 1: datasets to load
%%%% Copy paste the folders of choice here, add more or discard if needed
PATHNAME2{1} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\1_PT221-BJ\';
PATHNAME2{2} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\2_PT23-LS\';
PATHNAME2{3} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\3_PT68-DD\';
PATHNAME2{4} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\4_PT47-MR\';
PATHNAME2{5} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\5_PT115-AG\';
PATHNAME2{6} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\6_PT151-NJ\';
PATHNAME2{7} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\7_PT193-SM\';
PATHNAME2{8} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\8_PT202-AP\';
PATHNAME2{9} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\9_PT208-KK\';
PATHNAME2{10}= 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\10_PT214-DT\';
PATHNAME2{11}= 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\11_PT6-RG_followup\';
PATHNAME2{12}= 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\12_PT242-JL\';
PATHNAME2{13}= 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\13_PT53-JJ_followup\';
FILENAME = 'data_done';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Probability mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask1 = probability_mask.matrix;
mask1_vox = probability_mask.vox;clear probability_mask

%%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
sizes = [size(mask1,1)+offset size(mask1,2)+offset size(mask1,3)+offset];% size(mask1,4) size(mask1,5)];
mask1b = zeros(sizes);
mask1b((offset+1):size(mask1b,1),(offset+1):size(mask1b,2),(offset+1):size(mask1b,3),:,:) = mask1;
mask1 = mask1b;clear mask1b

L1 = (mask1 ~= 0);
% create all coordinates
[x,y,z] = meshgrid((1:size(mask1,2)).* mask1_vox(2), ...
    (1:size(mask1,1)).* mask1_vox(1),(1:size(mask1,3)).* mask1_vox(3));
geo.x_coor_vel = x(L1);geo.y_coor_vel = y(L1);geo.z_coor_vel = z(L1);
contours = zeros(size(L1));
contours(L1==0) = -1;
contours(L1==1) = 1;
[F,V] = isosurface(contours,0); clear contours % make a surface from the detected contours
V = V .* (ones(size(V,1),1) * mask1_vox(1:3));
[geo.F,geo.V] = SmoothLaplacian(F,V,15); %laplacian smoothing for surface (Kevin Moerman)
clear F, clear V;

for n = 1:(size(PATHNAME1,2)+size(PATHNAME2,2))

    if n <= size(PATHNAME1,2)
        disp('Cohort #1')
        disp(['Aorta number ' num2str(n)])
        sswitch = 1;
        PATHNAME = PATHNAME1;
        PATHNAME1{n}
    elseif n > size(PATHNAME1,2)
        n = n-size(PATHNAME1,2);
        disp('Cohort #2')
        disp(['Aorta number ' num2str(n)])
        sswitch = 2;
        PATHNAME = PATHNAME2;
        PATHNAME2{n}
    end

    tic
    disp('...Busy loading data_done...')
    load(strcat(PATHNAME{n},FILENAME))
    data2 = data; clear data;
    disp('...Done loading data...')
    toc

    if ~isfield(data2,'PC_unaliased')
        data2.PC_unaliased = data2.PC_zeros;
    end
    mask2 = (data2.PC_unaliased(:,:,:,1,1) ~= 0);

    L2a = (mask2 ~= 0); % Using two logicals for the second mask since the size changes after the translation which is needed for coordinates but not for velocity

    %%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
    sizes = [size(mask2,1)+offset size(mask2,2)+offset size(mask2,3)+offset];% size(mask1,4) size(mask1,5)];
    mask2b = zeros(sizes);
    mask2b((offset+1):size(mask2b,1),(offset+1):size(mask2b,2),(offset+1):size(mask2b,3),:,:) = mask2;
    mask2 = mask2b;clear mask2b

    L2b = (mask2 ~= 0);
    % create velocity coordinates
    [x2,y2,z2] = meshgrid((1:size(mask2,2)).* data2.vox(2), ...
        (1:size(mask2,1)).*data2.vox(1),(1:size(mask2,3)).* data2.vox(3));
    data2.x_coor_vel = x2(L2b);data2.y_coor_vel = y2(L2b);data2.z_coor_vel = z2(L2b);
    clear x2, clear y2, clear z2
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
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',16);
    end

    [I,time] = find(mean_velo==max(mean_velo));
    clear mean_velo

    %%% Wall shear stress coordinates for both datasets
    geo.x_coor_wss = geo.V(:,1);
    geo.y_coor_wss = geo.V(:,2);
    geo.z_coor_wss = geo.V(:,3);

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
    end

    if calculateIE_Flag == 1;
        if ~exist(strcat(PATHNAME{n},'interpolation_error_ROI\mask1.mat'),'file')

            mkdir(PATHNAME{n},'interpolation_error_ROI')

            F1=figure('Name','Aorta shape: Paused after finishing a region so press space when finished!');
            patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.25);
            view([-180 -90]);axis ij;axis equal;axis off

            for i = 1:6
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo);
                region = getPosition(polyAAo);

                %          disp('saving, pausing')
                save(strcat([PATHNAME{n} 'interpolation_error_ROI\mask' num2str(i)]),'region');
                pause
            end
            close(F1)
        end
        load(strcat(PATHNAME{n},'interpolation_error_ROI\mask1'))
        data2_mask_AAo_inner_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_inner = mean(data2.vel_m(data2_mask_AAo_inner_vel));
        data2_mask_AAo_inner_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_inner = mean(data2.wss_m(data2_mask_AAo_inner_wss));
        load(strcat(PATHNAME{n},'interpolation_error_ROI\mask2'))
        data2_mask_AAo_outer_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_outer = mean(data2.vel_m(data2_mask_AAo_outer_vel));
        data2_mask_AAo_outer_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_outer = mean(data2.wss_m(data2_mask_AAo_outer_wss));
        load(strcat(PATHNAME{n},'interpolation_error_ROI\mask3'))
        data2_mask_arch_inner_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_inner = mean(data2.vel_m(data2_mask_arch_inner_vel));
        data2_mask_arch_inner_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_inner = mean(data2.wss_m(data2_mask_arch_inner_wss));
        load(strcat(PATHNAME{n},'interpolation_error_ROI\mask4'))
        data2_mask_arch_outer_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_outer = mean(data2.vel_m(data2_mask_arch_outer_vel));
        data2_mask_arch_outer_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_outer = mean(data2.wss_m(data2_mask_arch_outer_wss));
        load(strcat(PATHNAME{n},'interpolation_error_ROI\mask5'))
        data2_mask_DAo_inner_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_inner = mean(data2.vel_m(data2_mask_DAo_inner_vel));
        data2_mask_DAo_inner_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_inner = mean(data2.wss_m(data2_mask_DAo_inner_wss));
        load(strcat(PATHNAME{n},'interpolation_error_ROI\mask6'))
        data2_mask_DAo_outer_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_outer = mean(data2.vel_m(data2_mask_DAo_outer_vel));
        data2_mask_DAo_outer_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_outer = mean(data2.wss_m(data2_mask_DAo_outer_wss));

        mean_vel_before_interpolation(n,1) = mean_vel_asc_inner;
        mean_vel_before_interpolation(n,2) = mean_vel_asc_outer;
        mean_vel_before_interpolation(n,3) = mean_vel_arch_inner;
        mean_vel_before_interpolation(n,4) = mean_vel_arch_outer;
        mean_vel_before_interpolation(n,5) = mean_vel_DAo_inner;
        mean_vel_before_interpolation(n,6) = mean_vel_DAo_outer;

        mean_wss_before_interpolation(n,1) = mean_wss_asc_inner;
        mean_wss_before_interpolation(n,2) = mean_wss_asc_outer;
        mean_wss_before_interpolation(n,3) = mean_wss_arch_inner;
        mean_wss_before_interpolation(n,4) = mean_wss_arch_outer;
        mean_wss_before_interpolation(n,5) = mean_wss_DAo_inner;
        mean_wss_before_interpolation(n,6) = mean_wss_DAo_outer;

        if plotFlag == 1
            figure('Name','Velocity before interpolation: inner AAo')
            scatter3(data2.x_coor_vel(data2_mask_AAo_inner_vel),data2.y_coor_vel(data2_mask_AAo_inner_vel),data2.z_coor_vel(data2_mask_AAo_inner_vel),20,data2.vel_m(data2_mask_AAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','Velocity before interpolation: outer AAo')
            scatter3(data2.x_coor_vel(data2_mask_AAo_outer_vel),data2.y_coor_vel(data2_mask_AAo_outer_vel),data2.z_coor_vel(data2_mask_AAo_outer_vel),20,data2.vel_m(data2_mask_AAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','Velocity before interpolation: inner arch')
            scatter3(data2.x_coor_vel(data2_mask_arch_inner_vel),data2.y_coor_vel(data2_mask_arch_inner_vel),data2.z_coor_vel(data2_mask_arch_inner_vel),20,data2.vel_m(data2_mask_arch_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','Velocity before interpolation: outer arch')
            scatter3(data2.x_coor_vel(data2_mask_arch_outer_vel),data2.y_coor_vel(data2_mask_arch_outer_vel),data2.z_coor_vel(data2_mask_arch_outer_vel),20,data2.vel_m(data2_mask_arch_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','Velocity before interpolation: inner DAo')
            scatter3(data2.x_coor_vel(data2_mask_DAo_inner_vel),data2.y_coor_vel(data2_mask_DAo_inner_vel),data2.z_coor_vel(data2_mask_DAo_inner_vel),20,data2.vel_m(data2_mask_DAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','Velocity before interpolation: outer DAo')
            scatter3(data2.x_coor_vel(data2_mask_DAo_outer_vel),data2.y_coor_vel(data2_mask_DAo_outer_vel),data2.z_coor_vel(data2_mask_DAo_outer_vel),20,data2.vel_m(data2_mask_DAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','WSS before interpolation: inner AAo')
            scatter3(data2.x_coor_wss(data2_mask_AAo_inner_wss),data2.y_coor_wss(data2_mask_AAo_inner_wss),data2.z_coor_wss(data2_mask_AAo_inner_wss),20,data2.wss_m(data2_mask_AAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','WSS before interpolation: outer AAo')
            scatter3(data2.x_coor_wss(data2_mask_AAo_outer_wss),data2.y_coor_wss(data2_mask_AAo_outer_wss),data2.z_coor_wss(data2_mask_AAo_outer_wss),20,data2.wss_m(data2_mask_AAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','WSS before interpolation: inner arch')
            scatter3(data2.x_coor_wss(data2_mask_arch_inner_wss),data2.y_coor_wss(data2_mask_arch_inner_wss),data2.z_coor_wss(data2_mask_arch_inner_wss),20,data2.wss_m(data2_mask_arch_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','WSS before interpolation: outer arch')
            scatter3(data2.x_coor_wss(data2_mask_arch_outer_wss),data2.y_coor_wss(data2_mask_arch_outer_wss),data2.z_coor_wss(data2_mask_arch_outer_wss),20,data2.wss_m(data2_mask_arch_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','WSS before interpolation: inner DAo')
            scatter3(data2.x_coor_wss(data2_mask_DAo_inner_wss),data2.y_coor_wss(data2_mask_DAo_inner_wss),data2.z_coor_wss(data2_mask_DAo_inner_wss),20,data2.wss_m(data2_mask_DAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','WSS before interpolation: outer DAo')
            scatter3(data2.x_coor_wss(data2_mask_DAo_outer_wss),data2.y_coor_wss(data2_mask_DAo_outer_wss),data2.z_coor_wss(data2_mask_DAo_outer_wss),20,data2.wss_m(data2_mask_DAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
        end
    end

    if plotFlag == 1
        figure('Name',strcat('To be registered aorta ',num2str(n)))
        plot3(geo.x_coor_vel,geo.y_coor_vel,geo.z_coor_vel,'r.')
        hold on
        % plot3(geo.x_coor_wss,geo.y_coor_wss,geo.z_coor_wss,'g.')
        plot3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,'b.')
        % plot3(data2.x_coor_wss,data2.y_coor_wss,data2.z_coor_wss,'k.')
        legend('to remain the same','to be transformed')
        axis equal; axis off;view([-180 -90]); axis ij
    end

    %%% Registration
    PSF = fspecial('gaussian',1,1);
    mask1_to_register = mask1;
    mask1_to_register = imfilter(mask1_to_register,PSF,'conv');

    mask2 = imfilter(mask2,PSF,'conv');

    disp(' ')
    disp('...Busy registering...Affine registration (dof = 12)! So only scalars used (Affine registration doesnt work for vectors)')
    disp(' ')
    disp('...This can take up to 5 minutes...')

    tic
    % directory with flirt.exe and cygwin1.dll (use cygwin convention)
    fsldir = '/cygdrive/c/1_Chicago/WSSprojectWithAmsterdam/flirt/';

    % save as nifti
    cnii=make_nii(mask1_to_register,[mask1_vox(1) mask1_vox(2) mask1_vox(3)]);
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
    system('c:\cygwin\bin\bash runflirt.sh');

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
    x_coor_vel = yxz_coor_vel_new(2,:)';
    y_coor_vel = yxz_coor_vel_new(1,:)';
    z_coor_vel = yxz_coor_vel_new(3,:)';clear yxz_coor_vel_new

    %%% WSS coordinates
    yxz_coor_wss = [data2.y_coor_wss data2.x_coor_wss data2.z_coor_wss];
    yxz_coor_wss(:,4) = 1;
    yxz_coor_wss_new = inv(worldmat)*yxz_coor_wss'; clear yxz_coor_wss
    x_coor_wss = yxz_coor_wss_new(2,:)';
    y_coor_wss = yxz_coor_wss_new(1,:)';
    z_coor_wss = yxz_coor_wss_new(3,:)';clear yxz_coor_wss_new
    toc
    disp('..Done registering...')

    % NOTE HERE THAT VELOCITY VECTORS AND WSS VECTORS ARE NOT ROTATED!
    % That's because affine registration includes not only translating and rotating
    % but also SCALING. We don't want our velocity and WSS vectors to be scaled.
    % That's why we work only with scalars in this code. Registration of the velocity
    % and WSS vectors is possible, but then only RIGID registration should be used.
    % If you want this, then you should use the function 'make_atlas_point_cloud_vectors_rigid_registration',
    % which can also be found in this tool. However, I expect that the interpolation works better
    % with affine registration. PvO

    if plotFlag == 1
        figure('Name','Registered')
        plot3(geo.x_coor_vel,geo.y_coor_vel,geo.z_coor_vel,'r.')
        hold on
        %  plot3(geo.x_coor_wss,geo.y_coor_wss,geo.z_coor_wss,'g.')
        plot3(x_coor_vel,y_coor_vel,z_coor_vel,'b.')
        %  plot3(x_coor_wss,y_coor_wss,z_coor_wss,'k.')
        legend('Prob mask','Aorta')
        %legend('probability mask/atlas','individual aorta')
        axis equal; axis off;view([-180 -90]); axis ij
        pause(10)

        figure('Name','transformed velocity')
        patch('Faces',data2.F_matrix{1},'Vertices',[x_coor_wss y_coor_wss z_coor_wss], ...
            'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.1);
        hold on
        scatter3(x_coor_vel,y_coor_vel,z_coor_vel,50,data2.vel_m,'filled')
        colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])

        figure('Name','transformed WSS')
        patch('Faces',data2.F_matrix{1},'Vertices',[x_coor_wss y_coor_wss z_coor_wss], ...
            'EdgeColor','none','FaceVertexCData',data2.wss_m,'FaceColor','interp','FaceAlpha',1);
        colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    end

    if calculateRE_Flag == 1
        % Transform point cloud back to matrix by rounding the coordinates
        x_vel_round = round(x_coor_vel./mask1_vox(1));
        y_vel_round = round(y_coor_vel./mask1_vox(2));
        z_vel_round = round(z_coor_vel./mask1_vox(3));

        % remove doubles ue to rounding and put them in a matrix
        indices_mask2 = [x_vel_round y_vel_round z_vel_round];
        siz=max(indices_mask2,[],1);
        IND = sub2ind(siz,x_vel_round,y_vel_round,z_vel_round);
        [b, IND_double_removed, m] = unique(IND);
        clear b, clear m
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

        % due to smoothing values of mask2_new are between 0 and 1, so set
        % everything > 0 to 1
        L_mask2 = double(mask2_new ~= 0);

        % calculate where voxels exist for both masks
        difference = abs(mask1-L_mask2);
        [I1,J] = find(mask1~=0);
        [I2,J] = find(L_mask2~=0);
        mean_I = (size(I1,1) + size(I2,1))/2;
        [I_diff,J] = find(difference~=0);
        diff_voxels = size(I_diff,1);
        diff_percentage(n,1) = ((diff_voxels / mean_I) * 100)/2;
        disp(' ')
        disp(['RE: Difference between aorta and atlas geometry = ' num2str(diff_percentage(n,1))]);
        disp(' ')
    end

    %%% Interpolate velocities to co-registered coordinates
    interpolation_function = TriScatteredInterp([x_coor_vel y_coor_vel z_coor_vel],data2.x_value_vel,'nearest');
    x_value_vel = interpolation_function([geo.x_coor_vel geo.y_coor_vel geo.z_coor_vel]);
    x_value_vel(isnan(x_value_vel)) = 0;

    interpolation_function = TriScatteredInterp([x_coor_vel y_coor_vel z_coor_vel],data2.y_value_vel,'nearest');
    y_value_vel = interpolation_function([geo.x_coor_vel geo.y_coor_vel geo.z_coor_vel]);
    y_value_vel(isnan(y_value_vel)) = 0;

    interpolation_function = TriScatteredInterp([x_coor_vel y_coor_vel z_coor_vel],data2.z_value_vel,'nearest');
    z_value_vel = interpolation_function([geo.x_coor_vel geo.y_coor_vel geo.z_coor_vel]);
    z_value_vel(isnan(z_value_vel)) = 0;

    data2.vel_m = sqrt(x_value_vel.^2+y_value_vel.^2+z_value_vel.^2);

    %%% Interpolate WSS to co-registered coordinates
    interpolation_function = TriScatteredInterp([x_coor_wss y_coor_wss z_coor_wss],data2.x_value_wss,'nearest');
    x_value_wss = interpolation_function([geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss]);
    x_value_wss(isnan(x_value_wss)) = 0;

    interpolation_function = TriScatteredInterp([x_coor_wss y_coor_wss z_coor_wss],data2.y_value_wss,'nearest');
    y_value_wss = interpolation_function([geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss]);
    y_value_wss(isnan(y_value_wss)) = 0;

    interpolation_function = TriScatteredInterp([x_coor_wss y_coor_wss z_coor_wss],data2.z_value_wss,'nearest');
    z_value_wss = interpolation_function([geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss]);
    z_value_wss(isnan(z_value_wss)) = 0;

    data2.wss_m = sqrt(x_value_wss.^2+y_value_wss.^2+z_value_wss.^2);

    if plotFlag == 1
        figure('Name','interpolated to atlas velocity')
        scatter3(geo.x_coor_vel,geo.y_coor_vel,geo.z_coor_vel,50,data2.vel_m,'filled')
        colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])

        figure('Name','interpolated to atlas WSS')
        patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none','FaceVertexCData',data2.wss_m,'FaceColor','interp','FaceAlpha',1);
        colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    end

    if calculateIE_Flag == 1;

        if ~exist(strcat(PATHNAME_probability_mask,'interpolation_error_ROI\mask1.mat'),'file')

            F2=figure('Name','Atlas shape: Paused after finishing a region so press space when finished!');
            patch('Faces',geo.F,'Vertices',[geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss],'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.25);
            view([-180 -90]);axis ij;axis equal;axis off

            for i = 1:6
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo);
                region = getPosition(polyAAo);

                %          disp('saving, pausing')
                mkdir(PATHNAME_probability_mask,'interpolation_error_ROI')
                save(strcat([PATHNAME_probability_mask 'interpolation_error_ROI\mask' num2str(i)]),'region');
                pause
            end

            close(F2)
        end

        load(strcat(PATHNAME_probability_mask,'interpolation_error_ROI\mask1'))
        geo_mask_AAo_inner_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_inner = mean(data2.vel_m(geo_mask_AAo_inner_vel));
        geo_mask_AAo_inner_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_inner = mean(data2.wss_m(geo_mask_AAo_inner_wss));
        load(strcat(PATHNAME_probability_mask,'interpolation_error_ROI\mask2'))
        geo_mask_AAo_outer_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_outer = mean(data2.vel_m(geo_mask_AAo_outer_vel));
        geo_mask_AAo_outer_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_outer = mean(data2.wss_m(geo_mask_AAo_outer_wss));
        load(strcat(PATHNAME_probability_mask,'interpolation_error_ROI\mask3'))
        geo_mask_arch_inner_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_inner = mean(data2.vel_m(geo_mask_arch_inner_vel));
        geo_mask_arch_inner_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_inner = mean(data2.wss_m(geo_mask_arch_inner_wss));
        load(strcat(PATHNAME_probability_mask,'interpolation_error_ROI\mask4'))
        geo_mask_arch_outer_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_outer = mean(data2.vel_m(geo_mask_arch_outer_vel));
        geo_mask_arch_outer_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_outer = mean(data2.wss_m(geo_mask_arch_outer_wss));
        load(strcat(PATHNAME_probability_mask,'interpolation_error_ROI\mask5'))
        geo_mask_DAo_inner_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_inner = mean(data2.vel_m(geo_mask_DAo_inner_vel));
        geo_mask_DAo_inner_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_inner = mean(data2.wss_m(geo_mask_DAo_inner_wss));
        load(strcat(PATHNAME_probability_mask,'interpolation_error_ROI\mask6'))
        geo_mask_DAo_outer_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_outer = mean(data2.vel_m(geo_mask_DAo_outer_vel));
        geo_mask_DAo_outer_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_outer = mean(data2.wss_m(geo_mask_DAo_outer_wss));

        if plotFlag == 1
            figure('Name','Velocity before interpolation: inner AAo')
            scatter3(geo.x_coor_vel(geo_mask_AAo_inner_vel),geo.y_coor_vel(geo_mask_AAo_inner_vel),geo.z_coor_vel(geo_mask_AAo_inner_vel),20,data2.vel_m(geo_mask_AAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity before interpolation: outer AAo')
            scatter3(geo.x_coor_vel(geo_mask_AAo_outer_vel),geo.y_coor_vel(geo_mask_AAo_outer_vel),geo.z_coor_vel(geo_mask_AAo_outer_vel),20,data2.vel_m(geo_mask_AAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity before interpolation: inner arch')
            scatter3(geo.x_coor_vel(geo_mask_arch_inner_vel),geo.y_coor_vel(geo_mask_arch_inner_vel),geo.z_coor_vel(geo_mask_arch_inner_vel),20,data2.vel_m(geo_mask_arch_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity before interpolation: outer arch')
            scatter3(geo.x_coor_vel(geo_mask_arch_outer_vel),geo.y_coor_vel(geo_mask_arch_outer_vel),geo.z_coor_vel(geo_mask_arch_outer_vel),20,data2.vel_m(geo_mask_arch_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity before interpolation: inner DAo')
            scatter3(geo.x_coor_vel(geo_mask_DAo_inner_vel),geo.y_coor_vel(geo_mask_DAo_inner_vel),geo.z_coor_vel(geo_mask_DAo_inner_vel),20,data2.vel_m(geo_mask_DAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity before interpolation: outer DAo')
            scatter3(geo.x_coor_vel(geo_mask_DAo_outer_vel),geo.y_coor_vel(geo_mask_DAo_outer_vel),geo.z_coor_vel(geo_mask_DAo_outer_vel),20,data2.vel_m(geo_mask_DAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS before interpolation: inner AAo')
            scatter3(geo.x_coor_wss(geo_mask_AAo_inner_wss),geo.y_coor_wss(geo_mask_AAo_inner_wss),geo.z_coor_wss(geo_mask_AAo_inner_wss),20,data2.wss_m(geo_mask_AAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS before interpolation: outer AAo')
            scatter3(geo.x_coor_wss(geo_mask_AAo_outer_wss),geo.y_coor_wss(geo_mask_AAo_outer_wss),geo.z_coor_wss(geo_mask_AAo_outer_wss),20,data2.wss_m(geo_mask_AAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS before interpolation: inner arch')
            scatter3(geo.x_coor_wss(geo_mask_arch_inner_wss),geo.y_coor_wss(geo_mask_arch_inner_wss),geo.z_coor_wss(geo_mask_arch_inner_wss),20,data2.wss_m(geo_mask_arch_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS before interpolation: outer arch')
            scatter3(geo.x_coor_wss(geo_mask_arch_outer_wss),geo.y_coor_wss(geo_mask_arch_outer_wss),geo.z_coor_wss(geo_mask_arch_outer_wss),20,data2.wss_m(geo_mask_arch_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS before interpolation: inner DAo')
            scatter3(geo.x_coor_wss(geo_mask_DAo_inner_wss),geo.y_coor_wss(geo_mask_DAo_inner_wss),geo.z_coor_wss(geo_mask_DAo_inner_wss),20,data2.wss_m(geo_mask_DAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS before interpolation: outer DAo')
            scatter3(geo.x_coor_wss(geo_mask_DAo_outer_wss),geo.y_coor_wss(geo_mask_DAo_outer_wss),geo.z_coor_wss(geo_mask_DAo_outer_wss),20,data2.wss_m(geo_mask_DAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        end

        mean_vel_after_interpolation(n,1) = mean_vel_asc_inner;
        mean_vel_after_interpolation(n,2) = mean_vel_asc_outer;
        mean_vel_after_interpolation(n,3) = mean_vel_arch_inner;
        mean_vel_after_interpolation(n,4) = mean_vel_arch_outer;
        mean_vel_after_interpolation(n,5) = mean_vel_DAo_inner;
        mean_vel_after_interpolation(n,6) = mean_vel_DAo_outer;
        mean_wss_after_interpolation(n,1) = mean_wss_asc_inner;
        mean_wss_after_interpolation(n,2) = mean_wss_asc_outer;
        mean_wss_after_interpolation(n,3) = mean_wss_arch_inner;
        mean_wss_after_interpolation(n,4) = mean_wss_arch_outer;
        mean_wss_after_interpolation(n,5) = mean_wss_DAo_inner;
        mean_wss_after_interpolation(n,6) = mean_wss_DAo_outer;

        IE_inner_AAo_vel = abs(mean_vel_before_interpolation(n,1)-mean_vel_after_interpolation(n,1)) / ((mean_vel_before_interpolation(n,1)+mean_vel_after_interpolation(n,1))./2)*100;
        IE_outer_AAo_vel = abs(mean_vel_before_interpolation(n,2)-mean_vel_after_interpolation(n,2)) / ((mean_vel_before_interpolation(n,2)+mean_vel_after_interpolation(n,2))./2)*100;
        IE_inner_asc_vel = abs(mean_vel_before_interpolation(n,3)-mean_vel_after_interpolation(n,3)) / ((mean_vel_before_interpolation(n,3)+mean_vel_after_interpolation(n,3))./2)*100;
        IE_outer_asc_vel = abs(mean_vel_before_interpolation(n,4)-mean_vel_after_interpolation(n,4)) / ((mean_vel_before_interpolation(n,4)+mean_vel_after_interpolation(n,4))./2)*100;
        IE_inner_DAo_vel = abs(mean_vel_before_interpolation(n,5)-mean_vel_after_interpolation(n,5)) / ((mean_vel_before_interpolation(n,5)+mean_vel_after_interpolation(n,5))./2)*100;
        IE_outer_DAo_vel = abs(mean_vel_before_interpolation(n,6)-mean_vel_after_interpolation(n,6)) / ((mean_vel_before_interpolation(n,6)+mean_vel_after_interpolation(n,6))./2)*100;
        IE_inner_AAo_wss = abs(mean_wss_before_interpolation(n,1)-mean_wss_after_interpolation(n,1)) / ((mean_wss_before_interpolation(n,1)+mean_wss_after_interpolation(n,1))./2)*100;
        IE_outer_AAo_wss = abs(mean_wss_before_interpolation(n,2)-mean_wss_after_interpolation(n,2)) / ((mean_wss_before_interpolation(n,2)+mean_wss_after_interpolation(n,2))./2)*100;
        IE_inner_asc_wss = abs(mean_wss_before_interpolation(n,3)-mean_wss_after_interpolation(n,3)) / ((mean_wss_before_interpolation(n,3)+mean_wss_after_interpolation(n,3))./2)*100;
        IE_outer_asc_wss = abs(mean_wss_before_interpolation(n,4)-mean_wss_after_interpolation(n,4)) / ((mean_wss_before_interpolation(n,4)+mean_wss_after_interpolation(n,4))./2)*100;
        IE_inner_DAo_wss = abs(mean_wss_before_interpolation(n,5)-mean_wss_after_interpolation(n,5)) / ((mean_wss_before_interpolation(n,5)+mean_wss_after_interpolation(n,5))./2)*100;
        IE_outer_DAo_wss = abs(mean_wss_before_interpolation(n,6)-mean_wss_after_interpolation(n,6)) / ((mean_wss_before_interpolation(n,6)+mean_wss_after_interpolation(n,6))./2)*100;

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
    end

    pause(10)
    close all

    if sswitch == 1;
        disp('VEL1')
        VEL1(:,n) = data2.vel_m;
        WSS1(:,n) = data2.wss_m;
    elseif sswitch == 2;
        disp('VEL2')
        VEL2(:,n) = data2.vel_m;
        WSS2(:,n) = data2.wss_m;
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

save C:\1_Chicago\VEL1 VEL1
save C:\1_Chicago\VEL2 VEL2
save C:\1_Chicago\WSS1 WSS1
save C:\1_Chicago\WSS2 WSS2

load C:\1_Chicago\VEL1 VEL1
load C:\1_Chicago\VEL2 VEL2
load C:\1_Chicago\WSS1 WSS1
load C:\1_Chicago\WSS2 WSS2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VELOCITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cohort1_vel = zeros([size(PATHNAME1,2) 1]);
cohort2_vel = zeros([size(PATHNAME2,2) 1]);
new_mask_red_ = zeros([size(VEL1,1) 1]);
new_mask_blue_ = zeros([size(VEL1,1) 1]);
new_mask_gray_ = zeros([size(VEL1,1) 1]);
alpha = 0.05;
disp('...Busy creating P-value map for velocity...')
tic
for i=1:size(VEL1,1)
    for n1 = 1:size(PATHNAME1,2)
        cohort1_vel(n1,1) = (VEL1(i,n1));
    end
    for n2 = 1:size(PATHNAME2,2)
        cohort2_vel(n2,1) = (VEL2(i,n2));
    end
    
    [p,h] = ranksum(cohort1_vel,cohort2_vel,alpha);
    if p < alpha && mean(cohort2_vel) > mean(cohort1_vel)
        new_mask_red_(i,1) = 1;
        new_mask_blue_(i,1) = -1;
        new_mask_gray_(i,1) = -1;
    elseif p < alpha && mean(cohort2_vel) < mean(cohort1_vel)
        new_mask_red_(i,1) = -1;
        new_mask_blue_(i,1) = 1;
        new_mask_gray_(i,1) = -1;
    elseif p > alpha && mean(cohort2_vel) > mean(cohort1_vel)
        new_mask_red_(i,1) = -1;
        new_mask_blue_(i,1) = -1;
        new_mask_gray_(i,1) = 1;
    elseif p > alpha && mean(cohort2_vel) < mean(cohort1_vel)
        new_mask_red_(i,1) = -1;
        new_mask_blue_(i,1) = -1;
        new_mask_gray_(i,1) = 1;
    end
end
toc

% translate point cloud into matrix
new_mask_red = zeros(size(L1));
new_mask_red(L1) = new_mask_red_;
new_mask_red(~L1) = -1;
new_mask_blue = zeros(size(L1));
new_mask_blue(L1) = new_mask_blue_;
new_mask_blue(~L1) = -1;
new_mask_gray = zeros(size(L1));
new_mask_gray(L1) = new_mask_gray_;
new_mask_gray(~L1) = -1;

if calculate_p_value_volumeFlag == 1
    
    [I,J] = find(L1~=0);
    total_volume = mask1_vox(1)*mask1_vox(2)*mask1_vox(3)*size(I,1);
    
    [I,J] = find(new_mask_red==1);
    red_volume = mask1_vox(1)*mask1_vox(2)*mask1_vox(3)*size(I,1);
    percentage_red_volume = red_volume / total_volume * 100;
    [I,J] = find(new_mask_blue==1);
    blue_volume = mask1_vox(1)*mask1_vox(2)*mask1_vox(3)*size(I,1);
    percentage_blue_volume = blue_volume / total_volume * 100;
    [I,J] = find(new_mask_gray==1);
    gray_volume = mask1_vox(1)*mask1_vox(2)*mask1_vox(3)*size(I,1);
    percentage_gray_volume = gray_volume / total_volume * 100;
    total_percentage = percentage_red_volume + percentage_blue_volume + percentage_gray_volume;
    
    disp(['Red volume percentage of total aorta = ' num2str(percentage_red_volume) ' %'])
    disp(['Blue volume percentage of total aorta = ' num2str(percentage_blue_volume) ' %'])
    disp(['Gray volume percentage of total aorta = ' num2str(percentage_gray_volume) ' %'])
    disp(['Total percentage of total aorta = ' num2str(total_percentage) ' %'])
    disp(' ')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% WSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cohort1_wss = zeros([size(PATHNAME1,2) 1]);
cohort2_wss = zeros([size(PATHNAME2,2) 1]);
p_value_wss = zeros([size(WSS1,1) 1]);
alpha = 0.05;
disp('...Busy creating P-value map for WSS...')
tic
for i=1:size(WSS1,1)
    for n1 = 1:size(PATHNAME1,2)
        cohort1_wss(n1,1) = (WSS1(i,n1));
    end
    for n2 = 1:size(PATHNAME2,2)
        cohort2_wss(n2,1) = (WSS2(i,n2));
    end
    
    [p,h] = ranksum(cohort1_wss,cohort2_wss,alpha);
    if p < alpha && mean(cohort2_wss) > mean(cohort1_wss)
        p_value_wss(i,1) = 0;
    elseif p < alpha && mean(cohort2_wss) < mean(cohort1_wss)
        p_value_wss(i,1) = 1;
    elseif p > alpha && mean(cohort2_wss) > mean(cohort1_wss)
        p_value_wss(i,1) = 2;
    elseif p > alpha && mean(cohort2_wss) < mean(cohort1_wss)
        p_value_wss(i,1) = 3;
    end
end
toc

%Create colormap for the p-value map
color = colormap(flipud(jet));
color2 = ones(size(color));
color2(1:16,1) = color2(1:16,1).*1;
color2(1:16,2) = color2(1:16,2).*0;
color2(1:16,3) = color2(1:16,3).*0;
color2(17:32,1) = color2(17:32,1).*0;
color2(17:32,2) = color2(17:32,2).*0;
color2(17:32,3) = color2(17:32,3).*1;
color2(33:64,1) = color2(33:64,1).*0.5;
color2(33:64,2) = color2(33:64,2).*0.5;
color2(33:64,3) = color2(33:64,3).*0.5;

if calculate_area_of_significance_wss == 1;
    if ~exist(strcat(PATHNAME_probability_mask,'six_regions_not_transformed\mask1.mat'),'file')
        
        patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',1);
        view([-180 -90]);axis ij;axis equal;axis off
        
        mkdir(PATHNAME_probability_mask,'six_regions_not_transformed')
        
        for i = 1:6
            %Polygon and mask for AAo
            polyAAo = impoly;
            wait(polyAAo);
            region = getPosition(polyAAo);
            %mask = inpolygon(x, y, AAo(:,1), AAo(:,2));
            
            disp('saving, pausing')
            save(strcat([PATHNAME_probability_mask 'six_regions_not_transformed\mask' num2str(i)]),'region');
            pause
        end
    end
    
    x_geo = geo.V(:,1);y_geo = geo.V(:,2);z_geo = geo.V(:,3);
    
    load(strcat(PATHNAME_probability_mask,'six_regions_not_transformed\mask1'));
    atlas_mask_AAo_inner = inpolygon(x_geo,y_geo, region(:,1), region(:,2));
    load(strcat(PATHNAME_probability_mask,'six_regions_not_transformed\mask2'));
    atlas_mask_AAo_outer = inpolygon(x_geo,y_geo, region(:,1), region(:,2));
    load(strcat(PATHNAME_probability_mask,'six_regions_not_transformed\mask3'));
    atlas_mask_arch_inner = inpolygon(x_geo,y_geo, region(:,1), region(:,2));
    load(strcat(PATHNAME_probability_mask,'six_regions_not_transformed\mask4'));
    atlas_mask_arch_outer = inpolygon(x_geo,y_geo, region(:,1), region(:,2));
    load(strcat(PATHNAME_probability_mask,'six_regions_not_transformed\mask5'));
    atlas_mask_DAo_inner = inpolygon(x_geo,y_geo, region(:,1), region(:,2));
    load(strcat(PATHNAME_probability_mask,'six_regions_not_transformed\mask6'));
    atlas_mask_DAo_outer = inpolygon(x_geo,y_geo, region(:,1), region(:,2));
    
    heat_asc1 = p_value_wss(atlas_mask_AAo_inner);
    [I1,J1] = find(heat_asc1 == 1);
    [I2,J2] = find(heat_asc1 == 0);
    percentage_significant_higher = size(I2,1) / size(heat_asc1,1) * 100;
    percentage_significant_lower = size(I1,1) / size(heat_asc1,1) * 100;
    
    disp(['Percentage of red significance inner AAo = ' num2str(percentage_significant_higher) '%'])
    disp(['Percentage of blue significance inner AAo = ' num2str(percentage_significant_lower) '%'])
    
    heat_asc2 = p_value_wss(atlas_mask_AAo_outer);
    [I1,J1] = find(heat_asc2 == 1);
    [I2,J2] = find(heat_asc2 == 0);
    percentage_significant_higher = size(I2,1) / size(heat_asc2,1) * 100;
    percentage_significant_lower = size(I1,1) / size(heat_asc2,1) * 100;
    
    disp(['Percentage of red significance outer AAo = ' num2str(percentage_significant_higher) '%'])
    disp(['Percentage of blue significance outer AAo = ' num2str(percentage_significant_lower) '%'])
    
    heat_arch1 = p_value_wss(atlas_mask_arch_inner);
    [I1,J1] = find(heat_arch1 == 1);
    [I2,J2] = find(heat_arch1 == 0);
    percentage_significant_higher = size(I2,1) / size(heat_arch1,1) * 100;
    percentage_significant_lower = size(I1,1) / size(heat_arch1,1) * 100;
    
    disp(['Percentage of red significance inner arch = ' num2str(percentage_significant_higher) '%'])
    disp(['Percentage of blue significance inner arch = ' num2str(percentage_significant_lower) '%'])
    
    heat_arch2 = p_value_wss(atlas_mask_arch_outer);
    [I1,J1] = find(heat_arch2 == 1);
    [I2,J2] = find(heat_arch2 == 0);
    percentage_significant_higher = size(I2,1) / size(heat_arch2,1) * 100;
    percentage_significant_lower = size(I1,1) / size(heat_arch2,1) * 100;
    
    disp(['Percentage of red significance outer arch = ' num2str(percentage_significant_higher) '%'])
    disp(['Percentage of blue significance outer arch = ' num2str(percentage_significant_lower) '%'])
    
    heat_desc1 = p_value_wss(atlas_mask_DAo_inner);
    [I1,J1] = find(heat_desc1 == 1);
    [I2,J2] = find(heat_desc1 == 0);
    percentage_significant_higher = size(I2,1) / size(heat_desc1,1) * 100;
    percentage_significant_lower = size(I1,1) / size(heat_desc1,1) * 100;
       
    disp(['Percentage of red significance inner DAo = ' num2str(percentage_significant_higher) '%'])
    disp(['Percentage of blue significance inner DAo = ' num2str(percentage_significant_lower) '%'])
        
    heat_desc2 = p_value_wss(atlas_mask_DAo_outer);
    [I1,J1] = find(heat_desc2 == 1);
    [I2,J2] = find(heat_desc2 == 0);
    percentage_significant_higher = size(I2,1) / size(heat_desc2,1) * 100;
    percentage_significant_lower = size(I1,1) / size(heat_desc2,1) * 100;
    
    disp(['Percentage of red significance outer DAo = ' num2str(percentage_significant_higher) '%'])
    disp(['Percentage of blue significance outer DAo = ' num2str(percentage_significant_lower) '%'])
    
    heat_total = p_value_wss;
    [I1,J1] = find(heat_total == 1);
    [I2,J2] = find(heat_total == 0);
    percentage_significant_higher = size(I2,1) / size(heat_total,1) * 100;
    percentage_significant_lower = size(I1,1) / size(heat_total,1) * 100;    
          
    disp(['Percentage of red significance TOTAL aorta = ' num2str(percentage_significant_higher) '%'])
    disp(['Percentage of blue significance TOTAL aorta = ' num2str(percentage_significant_lower) '%'])        
    
    if plotFlag == 1
        figure('Name','higher/lower: inner AAo')
        scatter3(x_geo(atlas_mask_AAo_inner),y_geo(atlas_mask_AAo_inner),z_geo(atlas_mask_AAo_inner),20,heat_asc1,'filled');axis equal;colormap(color2);colorbar;caxis([0 4])
        xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
        figure('Name','higher/lower: outer AAo')
        scatter3(x_geo(atlas_mask_AAo_outer),y_geo(atlas_mask_AAo_outer),z_geo(atlas_mask_AAo_outer),20,heat_asc2,'filled');axis equal;colormap(color2);colorbar;caxis([0 4])
        xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
        figure('Name','higher/lower: inner arch')
        scatter3(x_geo(atlas_mask_arch_inner),y_geo(atlas_mask_arch_inner),z_geo(atlas_mask_arch_inner),20,heat_arch1,'filled');axis equal;colormap(color2);colorbar;caxis([0 4])
        xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
        figure('Name','higher/lower: outer arch')
        scatter3(x_geo(atlas_mask_arch_outer),y_geo(atlas_mask_arch_outer),z_geo(atlas_mask_arch_outer),20,heat_arch2,'filled');axis equal;colormap(color2);colorbar;caxis([0 4])
        xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
        figure('Name','higher/lower: inner DAo')
        scatter3(x_geo(atlas_mask_DAo_inner),y_geo(atlas_mask_DAo_inner),z_geo(atlas_mask_DAo_inner),20,heat_desc1,'filled');axis equal;colormap(color2);colorbar;caxis([0 4])
        xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
        figure('Name','higher/lower: outer DAo')
        scatter3(x_geo(atlas_mask_DAo_outer),y_geo(atlas_mask_DAo_outer),z_geo(atlas_mask_DAo_outer),20,heat_desc2,'filled');axis equal;colormap(color2);colorbar;caxis([0 4])
        xlabel('x'),ylabel('y'),zlabel('z');view([180 -90]); axis ij
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VELOCITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count1 = 0;
angles(1) = 0;
f1 = figure('Name','Traffic Light Map');
contours = zeros(size(L1));
contours(L1==0) = -1;
contours(L1==1) = 1;
[F,V] = isosurface(smooth3(contours),0); % make a surface from the detected contours
patch('Faces',F,'Vertices',[V(:,1)-offset V(:,2)-offset V(:,3)-offset],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.1);
hold on
[F1,V1] = isosurface(x./mask1_vox(1)-offset,y./mask1_vox(2)-offset,z./mask1_vox(3)-offset,smooth3(new_mask_red),0);
[F2,V2] = isosurface(x./mask1_vox(1)-offset,y./mask1_vox(2)-offset,z./mask1_vox(3)-offset,smooth3(new_mask_blue),0);
p11=patch('Faces',F1,'Vertices',V1,'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',1);
p12=patch('Faces',F2,'Vertices',V2,'EdgeColor','none','FaceColor',[0 0 1],'FaceAlpha',1);
%p13=patch('Faces',F3,'Vertices',V3,'EdgeColor','none','FaceColor',[0 1 0],'FaceAlpha',1);
axis equal; axis off;axis ij;caxis([0 4])
view([-180 -90])

% set up results folder
dir_orig = pwd;
dir_new = PATHNAME_probability_mask; cd(dir_new); %cd('..')
%dir_new = pwd;
mkdir('results_pvalue_map')
dir_new
dir_new = strcat(dir_new,'results_pvalue_map');
saveas(gcf,[dir_new '\pvalue_map_velocity.fig'])

aspectRatio = 1./mask1_vox;
set(gca,'dataaspectRatio',aspectRatio(1:3))
camlight headlight;camlight(180,180); lighting phong

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
    'String','Show Red Region')
uicontrol('Style','checkbox',...
    'Value',1, 'Position', [15 275 20 20], ...
    'Callback', {@show_red_region,gca});

uicontrol('Style','text',...
    'Position',[15 250 120 20],...
    'String','Show Blue Map')
uicontrol('Style','checkbox',...
    'Value',1, 'Position', [10 250 20 20], ...
    'Callback', {@show_blue_region,gca});

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

    function show_blue_region(hObj,event,ax)
        show = round(get(hObj,'Value'));
        if show == 1
            patchobj = findobj(p12);
            set(patchobj,'HandleVisibility','on','Visible','on');
        elseif show == 0
            patchobj = findobj(p12);
            set(patchobj,'HandleVisibility','off','Visible','off');
        end
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
        dphi = 180;
        
        c=camlight(dthetas,dphi);
        lighting phong
    end

set(f1,'toolbar','figure');

%savefig(f2,strcat(dir_new,'\heat_map'))
print(f1,'-dtiff','-r600',strcat(dir_new,'\pvalue_map_velocity.tif'));
cd(dir_orig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% P-value map WSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count3 = 0;
angles(1) = 0;
f2 = figure('Name','Heat map');
contours = zeros(size(L1));
contours(L1==0) = -1;
contours(L1==1) = 1;
[F,V] = isosurface(smooth3(contours),0); % make a surface from the detected contours
patch('Faces',F,'Vertices',[V(:,1)-offset V(:,2)-offset V(:,3)-offset],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.1);
hold on
p21=patch('Faces',geo.F,'Vertices',[geo.V(:,1)./mask1_vox(1)-offset geo.V(:,2)./mask1_vox(2)-offset geo.V(:,3)./mask1_vox(3)-offset],'EdgeColor','none', 'FaceVertexCData',p_value_wss,'FaceColor','interp','FaceAlpha',1);
colormap(color2)
caxis([0 4]);
axis equal; axis ij; axis off;
view([-180 -90])
saveas(gcf,[dir_new '\pvalue_map_wss.fig'])
aspectRatio = 1./mask1_vox;
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

p_value_map.red_vel = new_mask_red;
p_value_map.blue_vel = new_mask_blue;
p_value_map.gray_vel = new_mask_gray;
p_value_map.vertices = [V(:,1)-offset V(:,2)-offset V(:,3)-offset];
p_value_map.faces = F;
p_value_map.p_value_wss = p_value_wss;

% save results in results folder
save(strcat(dir_new,'\p_value_map'),'p_value_map');
%savefig(f2,strcat(dir_new,'\heat_map'))
print(f2,'-dtiff','-r600',strcat(dir_new,'\p_value_wss.tif'));
cd(dir_orig)

end