function [atlas] = make_atlas_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,plotFlag,calculateRE_Flag,calculateIE_Flag,peak_systolicFlag)

%%% [atlas] = make_atlas_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,plotFlag,calculateRE_Flag,calculateIE_Flag,peak_systolicFlag)
%
% This function creates the velocity and WSS atlas (cohort-averaged velocity/WSS map is more politically correct)
% of the batch data put into this file (see below)
% The method developed was published in (so please cite this paper when using this method):
%
% A Methodology to Detect Abnormal Relative Wall Shear Stress on the Full Surface of the Thoracic Aorta Using 4D Flow MRI
% van Ooij P, Potters WV, Nederveen AJ, Allen BD, Collins J, Carr J, Malaisrie SC, Markl M, Barker AJ
% Magn Reson Med. 2014 Feb 25, doi: 10.1002/mrm.25224. [Epub ahead of print]
%
% First, each aorta (you can probably use this for carotids and intracranial vessels as well, but I haven't tried yet) is co-registered
% (affine registration) to the idealized aorta geometry created with the function 'make_geometry_point_cloud.m'. The velocity and WSS magnitude
% values are subsequently interpolated (nearest neighbour interpolation) to the velocity and WSS coordinates of the idealized geometry. Finally,
% an average and standard deviation (SD) of the input aortas is calculted resulting in the mean and SD velocity atlas and mean and SD WSS atlas
% that are saed to the directory of choice for further use.
%
% 2014, Pim van Ooij, Northwestern University
%
% Input
% 1)PATHNAME          : The PATHNAMES for the subjects that the atlas will consist of
% 2)PATHNAME_probability_mask  : The PATHNAME containing the probability mask (or idealized geometry) that is used for the atlas
% 2)plotFlag        : If plotFlag switched, Matlab will output any possible image to the screen to check if everything happens correctly.
% 3)calculateRE_Flag: When switched on the registration error will be calculated. However, this is a different RE than the one in the paper
%                     mentioned above as the RE in this function is calculated from AFFINE registration whereas the RE reported in the paper
%                     is calculated from RIGID registration in the function 'make_geometry_point_cloud.m'
% 4)calculateIE_Flag: When switched on the interpolation error (RE, see paper mentioned above) will be calculated. Note that ROIs are needed
%                     which can be drawn manually when switched on.
% 5)peak_systolicFlag: When switched on the atlas will be created for the peak systolic time frame only, default is the atlass for 5 systolic
%                     timesteps averaged.
%
% Output
% 1)atlas:           The mean and SD velocity and WSS atlas are saved to this struct
%
% Usage
% This code is for creating velocity and WSS atlases from 'data_done' structs as created by Pims_postprocessing tool
% The function for creating atlases from mrStructs is under construction
%
% Examples:
% [atlas] = make_atlas_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,plotFlag,calculateRE_Flag,calculateIE_Flag,peak_systolicFlag)
% [atlas] = make_atlas_point_cloud_scalars_affine_registration_NEW('','',0,0,0,0)
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

if nargin < 1 || isempty(PATHNAME)
    % In for-loop
end

if nargin < 2 || isempty(PATHNAME_probability_mask)
    [FILENAME,PATHNAME_probability_mask] = uigetfile('.mat','Load probability mask');
    load(strcat(PATHNAME_probability_mask,FILENAME))
else
    load([PATHNAME_probability_mask '\probability_mask' ])
end

if nargin < 3 || isempty(plotFlag)
    plotFlag = 1;
end

if nargin < 4 || isempty(calculateRE_Flag)
    calculateRE_Flag = 1;
end

if nargin < 5 || isempty(calculateIE_Flag)
    calculateIE_Flag = 0;
end

if nargin < 6 || isempty(peak_systolicFlag)
    peak_systolicFlag = 0;
end

%%%% datasets to load
%%%% Copy paste the folders of choice here, add more or discard if needed
FILENAME1 = 'mask_struct_aorta';        % 1: Load mask
FILENAME2 = 'vel_struct';               % 2: Load velocity
FILENAME3 = 'Wss_point_cloud_aorta';    % 3: Load WSS

mask1 = probability_mask.matrix;
mask1_vox = probability_mask.vox;clear probability_mask

if plotFlag == 1
    figure('Name','Geometry')
    L_figure = (squeeze(max(mask1,[],3))~=0);
    imagesc(squeeze(max(mask1,[],3)),'Alphadata',double(L_figure));
    colorbar
    axis tight; axis equal; axis off
end

L1 = (mask1 ~= 0);
% create all coordinates
[x,y,z] = meshgrid((1:size(mask1,2)).* mask1_vox(2), ...
    (1:size(mask1,1)).* mask1_vox(1),(1:size(mask1,3)).* mask1_vox(3));
geo.x_coor_vel = x(L1);geo.y_coor_vel = y(L1);geo.z_coor_vel = z(L1);
clear x, clear y, clear z
contours = zeros(size(L1));
contours(L1==0) = -1;
contours(L1==1) = 1;
[F,V] = isosurface(contours,0); clear contours % make a surface from the detected contours
V = V .* (ones(size(V,1),1) * mask1_vox(1:3));
[geo.F,geo.V] = SmoothLaplacian(F,V,15); %laplacian smoothing for surface (Kevin Moerman)
clear F, clear V;

% Allocate sizes
VELx = zeros(size(geo.x_coor_vel,1),size(PATHNAME,2));
VELy = zeros(size(geo.x_coor_vel,1),size(PATHNAME,2));
VELz = zeros(size(geo.x_coor_vel,1),size(PATHNAME,2));
WSSx = zeros(size(geo.V,1),size(PATHNAME,2));
WSSy = zeros(size(geo.V,1),size(PATHNAME,2));
WSSz = zeros(size(geo.V,1),size(PATHNAME,2));
for n = 1:size(PATHNAME,2)
    disp(['Aorta number ' num2str(n)])
    
    if nargin < 1 || isempty(PATHNAME)
        [FILENAME,PATHNAME{n}] = uigetfile('.mat','Load probability mask');
        load([PATHNAME{n} filesep 'mrstruct' filesep FILENAME1])
    else
        currDir=pwd;
        cd(PATHNAME{n})
        folders = ls;
        if exist('mrstruct','dir')==7
            MrstructPath = strcat(PATHNAME{n},'\mrstruct\');
        else
            for i=3:size(folders,1)
                [a,b]=find(folders(i,1:8)=='mrstruct');
                if sum(a)==8
                    MrstructPath=strcat(PATHNAME{n},'\',folders(i,:),'\');
                    break
                end
            end
        end
        cd(currDir);
        load([MrstructPath FILENAME1])
    end
    
    mask2 = mrstruct_mask.dataAy;
    mask2_vox = mrstruct_mask.vox;clear mrstruct_mask
    L2 = (mask2 ~= 0);
    
    % create velocity coordinates
    [x,y,z] = meshgrid((1:size(mask2,2)).* mask2_vox(2), ...
        (1:size(mask2,1)).*mask2_vox(1),(1:size(mask2,3)).* mask2_vox(3));
    data2.x_coor_vel = x(L2);data2.y_coor_vel = y(L2);data2.z_coor_vel = z(L2);
    clear x, clear y, clear z
    contours = zeros(size(L2));
    contours(L2==0) = -1;
    contours(L2==1) = 1;
    [F,V] = isosurface(contours,0); % make a surface from the detected contours
    V = V .* (ones(size(V,1),1) * mask2_vox(1:3));
    [data2.F,data2.V] = SmoothLaplacian(F,V,15); %laplacian smoothing for surface (Kevin Moerman)
    clear F, clear V
    
    load([MrstructPath FILENAME2])
    velocity = double(mrStruct.dataAy); clear mrStruct; % for interpolation velocity needs to be a double
    
    vx = velocity(:,:,:,1,:);vy = velocity(:,:,:,2,:);vz = velocity(:,:,:,3,:);    
    velocity(:,:,:,1,:) = vy;velocity(:,:,:,2,:) = vx;velocity(:,:,:,3,:) = vz;
    
    %%% Wall shear stress coordinates for both datasets
    geo.x_coor_wss = geo.V(:,1);
    geo.y_coor_wss = geo.V(:,2);
    geo.z_coor_wss = geo.V(:,3);
    
    data2.x_coor_wss = data2.V(:,1);
    data2.y_coor_wss = data2.V(:,2);
    data2.z_coor_wss = data2.V(:,3);
    
    load([MrstructPath FILENAME3])
    WSS = Wss_point_cloud; clear Wss_point_cloud
    
    % Find max velocity
    for t = 1:size(velocity,5)
        vx = squeeze(velocity(:,:,:,1,t));
        vy = squeeze(velocity(:,:,:,2,t));
        vz = squeeze(velocity(:,:,:,3,t));
        vmagn = sqrt(vx.^2 + vy.^2 + vz.^2);
        mean_velo(t) = mean(vmagn(L2));
    end
    
    if plotFlag == 1
        figure('Name','Mean velocity')
        plot(1:size(velocity,5),mean_velo,'-ro','LineWidth',5,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',16);
    end
    [I,time] = find(mean_velo==max(mean_velo));
    clear mean_velo
    
    %%% What follows is a horrible piece of code, so if you're reading this and feel like cleaning it up, please do,
    %%% I'll buy you a beer next time we meet. PvO
    peak_systolicFlag
    size(WSS,2)
    if peak_systolicFlag == 1
        data2.x_value_vel = velocity(:,:,:,1,time);
        data2.y_value_vel = velocity(:,:,:,2,time);
        data2.z_value_vel = velocity(:,:,:,3,time);
        data2.x_value_vel = data2.x_value_vel(L2);
        data2.y_value_vel = data2.y_value_vel(L2);
        data2.z_value_vel = data2.z_value_vel(L2);
        if size(WSS,2) == 1
            data2.x_value_wss = WSS(:,1);
            data2.y_value_wss = WSS(:,2);
            data2.z_value_wss = WSS(:,3);            
        elseif size(WSS,2) > 5
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
        end
    elseif peak_systolicFlag == 0
        % Velocity averaged over 5 systolic time frames
        if time == 2    % mistriggering: second time frame is peak systole, averaging over 5 timesteps is not possible
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
        figure('Name','data2 velocity')
        vel_matrix = zeros(size(mask2));
        L = (mask2~=0);
        vel_matrix(L) = data2.vel_m;
        L_figure = (squeeze(max(vel_matrix,[],3))~=0);
        imagesc(squeeze(max(vel_matrix,[],3)),'Alphadata',double(L_figure));
        colorbar;axis tight; axis equal; axis ij; axis off;caxis([0 1.5]);%view([180 -90])
        
        figure('Name','data2 WSS')
        patch('Faces',data2.F,'Vertices',data2.V, ...
            'EdgeColor','none','FaceVertexCData',data2.wss_m,'FaceColor','interp','FaceAlpha',1);
        colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    end
    
    if calculateIE_Flag == 1;
        if ~exist(strcat(PATHNAME{n},'\interpolation_error_ROI\mask1.mat'),'file')
            
            mkdir(PATHNAME{n},'\interpolation_error_ROI')
            
            F1=figure('Name','Aorta shape: Paused after finishing a region so press space when finished!');
            patch('Faces',data2.F,'Vertices',data2.V,'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.25);
            view([-180 -90]);axis ij;axis equal;axis off
            
            for i = 1:6
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo);
                region = getPosition(polyAAo);
                
                %          disp('saving, pausing')
                save(strcat([PATHNAME{n} '\interpolation_error_ROI\mask' num2str(i)]),'region');
                pause
            end
            close(F1)
        end
        load(strcat(PATHNAME{n},'\interpolation_error_ROI\mask1'))
        data2_mask_AAo_inner_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_inner = mean(data2.vel_m(data2_mask_AAo_inner_vel));
        data2_mask_AAo_inner_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_inner = mean(data2.wss_m(data2_mask_AAo_inner_wss));
        load(strcat(PATHNAME{n},'\interpolation_error_ROI\mask2'))
        data2_mask_AAo_outer_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_outer = mean(data2.vel_m(data2_mask_AAo_outer_vel));
        data2_mask_AAo_outer_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_outer = mean(data2.wss_m(data2_mask_AAo_outer_wss));
        load(strcat(PATHNAME{n},'\interpolation_error_ROI\mask3'))
        data2_mask_arch_inner_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_inner = mean(data2.vel_m(data2_mask_arch_inner_vel));
        data2_mask_arch_inner_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_inner = mean(data2.wss_m(data2_mask_arch_inner_wss));
        load(strcat(PATHNAME{n},'\interpolation_error_ROI\mask4'))
        data2_mask_arch_outer_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_outer = mean(data2.vel_m(data2_mask_arch_outer_vel));
        data2_mask_arch_outer_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_outer = mean(data2.wss_m(data2_mask_arch_outer_wss));
        load(strcat(PATHNAME{n},'\interpolation_error_ROI\mask5'))
        data2_mask_DAo_inner_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_inner = mean(data2.vel_m(data2_mask_DAo_inner_vel));
        data2_mask_DAo_inner_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_inner = mean(data2.wss_m(data2_mask_DAo_inner_wss));
        load(strcat(PATHNAME{n},'\interpolation_error_ROI\mask6'))
        data2_mask_DAo_outer_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_outer = mean(data2.vel_m(data2_mask_DAo_outer_vel));
        data2_mask_DAo_outer_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_outer = mean(data2.wss_m(data2_mask_DAo_outer_wss));
        
        mean_vel_total_before_interpolation(n,1) = mean(data2.vel_m)
        mean_wss_total_before_interpolation(n,1) = mean(data2.wss_m)
        
        mean_vel_before_interpolation(n,1) = mean_vel_asc_inner;
        mean_vel_before_interpolation(n,2) = mean_vel_asc_outer;
        mean_vel_before_interpolation(n,3) = mean_vel_arch_inner;
        mean_vel_before_interpolation(n,4) = mean_vel_arch_outer;
        mean_vel_before_interpolation(n,5) = mean_vel_DAo_inner;
        mean_vel_before_interpolation(n,6) = mean_vel_DAo_outer
        
        mean_wss_before_interpolation(n,1) = mean_wss_asc_inner;
        mean_wss_before_interpolation(n,2) = mean_wss_asc_outer;
        mean_wss_before_interpolation(n,3) = mean_wss_arch_inner;
        mean_wss_before_interpolation(n,4) = mean_wss_arch_outer;
        mean_wss_before_interpolation(n,5) = mean_wss_DAo_inner;
        mean_wss_before_interpolation(n,6) = mean_wss_DAo_outer
        
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
        %      plot3(geo.x_coor_wss,geo.y_coor_wss,geo.z_coor_wss,'g.')
        plot3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,'b.')
        %      plot3(data2.x_coor_wss,data2.y_coor_wss,data2.z_coor_wss,'k.')
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
    % directory with flirt.exe and cygwin1.dll as read from atlas_tool.cfg in C:\temp
    fsldir = [path_flirt '/']
    
    % save as nifti
    cnii=make_nii(mask1_to_register,[mask1_vox(1) mask1_vox(2) mask1_vox(3)]);
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
    % directory for cygwin as read from atlas_tool.cfg in C:\temp
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
    % If you want this, then this function should be adjusted to 'make_atlas_point_cloud_vectors_rigid_registration'.
    % However, I expect that the interpolation works better
    % with affine registration. PvO
    
    if plotFlag == 1
        figure('Name','Registered')
        plot3(geo.x_coor_vel,geo.y_coor_vel,geo.z_coor_vel,'r.')
        hold on
        %     plot3(geo.x_coor_wss,geo.y_coor_wss,geo.z_coor_wss,'g.')
        plot3(x_coor_vel,y_coor_vel,z_coor_vel,'b.')
        %     plot3(x_coor_wss,y_coor_wss,z_coor_wss,'k.')
        legend('Prob mask','Aorta')
        %legend('probability mask/atlas','individual aorta')
        axis equal; axis off;view([-180 -90]); axis ij
        pause(10)
        
        figure('Name','transformed WSS')
        patch('Faces',data2.F,'Vertices',[x_coor_wss y_coor_wss z_coor_wss], ...
            'EdgeColor','none','FaceVertexCData',data2.wss_m,'FaceColor','interp','FaceAlpha',1);
        colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    end
    
    if calculateRE_Flag == 1
        offset = 100;
        
        % Transform point cloud back to matrix by rounding the coordinates
        x_vel_round = round(x_coor_vel./mask1_vox(1)) + offset;
        y_vel_round = round(y_coor_vel./mask1_vox(2)) + offset;
        z_vel_round = round(z_coor_vel./mask1_vox(3)) + offset;
        
        % remove doubles ue to rounding and put them in a matrix
        indices_mask2 = [x_vel_round y_vel_round z_vel_round];
        siz=max(indices_mask2,[],1);
        IND = sub2ind(siz,x_vel_round,y_vel_round,z_vel_round);
        [b, IND_double_removed, m] = unique(IND);
        clear b, clear m
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
        %    mask1 = mask1b;clear mask2b
        
        % Make sure both masks have the same dimensions
        if size(mask1b,1) > size(mask_new,1)
            mask_new(size(mask_new,1):size(mask1b,1),:,:) = 0;
        elseif size(mask1b,1) < size(mask_new,1)
            mask_new(size(mask1b,1)+1:size(mask_new,1),:,:) = [];
        end
        if size(mask1b,2) > size(mask_new,2)
            mask_new(:,size(mask_new,2):size(mask1b,2),:) = 0;
        elseif size(mask1b,2) < size(mask_new,2)
            mask_new(:,size(mask1b,2)+1:size(mask_new,2),:) = [];
        end
        if size(mask1b,3) > size(mask_new,3)
            mask_new(:,:,size(mask_new,3):size(mask1b,3)) = 0;
        elseif size(mask1b,3) < size(mask_new,3)
            mask_new(:,:,size(mask1b,3)+1:size(mask_new,3)) = [];
        end
        
        % due to smoothing values of mask2_new are between 0 and 1, so set
        % everything > 0 to 1
        L_mask2 = double(mask_new ~= 0);
        
        % calculate where voxels exist for both masks
        difference = abs(mask1b-L_mask2);
        [I1,J] = find(mask1b~=0);
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
        vel_matrix = zeros(size(mask1));
        L = (mask1~=0);
        vel_matrix(L) = data2.vel_m;
        L_figure = (squeeze(max(vel_matrix,[],3))~=0);
        imagesc(squeeze(max(vel_matrix,[],3)),'Alphadata',double(L_figure));
        colorbar;axis tight; axis equal; axis ij; axis off;caxis([0 1.5]);%view([180 -90])
        
        figure('Name','interpolated to atlas WSS')
        patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none','FaceVertexCData',data2.wss_m ,'FaceColor','interp','FaceAlpha',1);
        colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    end
    
    if calculateIE_Flag == 1;
        
        if ~exist(strcat(PATHNAME_probability_mask,'\interpolation_error_ROI\mask1.mat'),'file')
            
            F2=figure('Name','Atlas shape: Paused after finishing a region so press space when finished!');
            patch('Faces',geo.F,'Vertices',[geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss],'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.25);
            view([-180 -90]);axis ij;axis equal;axis off
            
            for i = 1:6
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo);
                region = getPosition(polyAAo);
                
                %          disp('saving, pausing')
                mkdir(PATHNAME_probability_mask,'\interpolation_error_ROI')
                save(strcat([PATHNAME_probability_mask '\interpolation_error_ROI\mask' num2str(i)]),'region');
                pause
            end
            
            close(F2)
        end
        
        load(strcat(PATHNAME_probability_mask,'\interpolation_error_ROI\mask1'))
        geo_mask_AAo_inner_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_inner = mean(data2.vel_m(geo_mask_AAo_inner_vel));
        geo_mask_AAo_inner_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_inner = mean(data2.wss_m(geo_mask_AAo_inner_wss));
        load(strcat(PATHNAME_probability_mask,'\interpolation_error_ROI\mask2'))
        geo_mask_AAo_outer_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_outer = mean(data2.vel_m(geo_mask_AAo_outer_vel));
        geo_mask_AAo_outer_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_outer = mean(data2.wss_m(geo_mask_AAo_outer_wss));
        load(strcat(PATHNAME_probability_mask,'\interpolation_error_ROI\mask3'))
        geo_mask_arch_inner_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_inner = mean(data2.vel_m(geo_mask_arch_inner_vel));
        geo_mask_arch_inner_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_inner = mean(data2.wss_m(geo_mask_arch_inner_wss));
        load(strcat(PATHNAME_probability_mask,'\interpolation_error_ROI\mask4'))
        geo_mask_arch_outer_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_outer = mean(data2.vel_m(geo_mask_arch_outer_vel));
        geo_mask_arch_outer_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_outer = mean(data2.wss_m(geo_mask_arch_outer_wss));
        load(strcat(PATHNAME_probability_mask,'\interpolation_error_ROI\mask5'))
        geo_mask_DAo_inner_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_inner = mean(data2.vel_m(geo_mask_DAo_inner_vel));
        geo_mask_DAo_inner_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_inner = mean(data2.wss_m(geo_mask_DAo_inner_wss));
        load(strcat(PATHNAME_probability_mask,'\interpolation_error_ROI\mask6'))
        geo_mask_DAo_outer_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_outer = mean(data2.vel_m(geo_mask_DAo_outer_vel));
        geo_mask_DAo_outer_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_outer = mean(data2.wss_m(geo_mask_DAo_outer_wss));
        
        mean_vel_total_after_interpolation(n,1) = mean(data2.vel_m)
        mean_wss_total_after_interpolation(n,1) = mean(data2.wss_m)
        
        if plotFlag == 1
            figure('Name','Velocity after interpolation: inner AAo')
            scatter3(geo.x_coor_vel(geo_mask_AAo_inner_vel),geo.y_coor_vel(geo_mask_AAo_inner_vel),geo.z_coor_vel(geo_mask_AAo_inner_vel),20,data2.vel_m(geo_mask_AAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity after interpolation: outer AAo')
            scatter3(geo.x_coor_vel(geo_mask_AAo_outer_vel),geo.y_coor_vel(geo_mask_AAo_outer_vel),geo.z_coor_vel(geo_mask_AAo_outer_vel),20,data2.vel_m(geo_mask_AAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity after interpolation: inner arch')
            scatter3(geo.x_coor_vel(geo_mask_arch_inner_vel),geo.y_coor_vel(geo_mask_arch_inner_vel),geo.z_coor_vel(geo_mask_arch_inner_vel),20,data2.vel_m(geo_mask_arch_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity after interpolation: outer arch')
            scatter3(geo.x_coor_vel(geo_mask_arch_outer_vel),geo.y_coor_vel(geo_mask_arch_outer_vel),geo.z_coor_vel(geo_mask_arch_outer_vel),20,data2.vel_m(geo_mask_arch_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity after interpolation: inner DAo')
            scatter3(geo.x_coor_vel(geo_mask_DAo_inner_vel),geo.y_coor_vel(geo_mask_DAo_inner_vel),geo.z_coor_vel(geo_mask_DAo_inner_vel),20,data2.vel_m(geo_mask_DAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Velocity after interpolation: outer DAo')
            scatter3(geo.x_coor_vel(geo_mask_DAo_outer_vel),geo.y_coor_vel(geo_mask_DAo_outer_vel),geo.z_coor_vel(geo_mask_DAo_outer_vel),20,data2.vel_m(geo_mask_DAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS after interpolation: inner AAo')
            scatter3(geo.x_coor_wss(geo_mask_AAo_inner_wss),geo.y_coor_wss(geo_mask_AAo_inner_wss),geo.z_coor_wss(geo_mask_AAo_inner_wss),20,data2.wss_m(geo_mask_AAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS after interpolation: outer AAo')
            scatter3(geo.x_coor_wss(geo_mask_AAo_outer_wss),geo.y_coor_wss(geo_mask_AAo_outer_wss),geo.z_coor_wss(geo_mask_AAo_outer_wss),20,data2.wss_m(geo_mask_AAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS after interpolation: inner arch')
            scatter3(geo.x_coor_wss(geo_mask_arch_inner_wss),geo.y_coor_wss(geo_mask_arch_inner_wss),geo.z_coor_wss(geo_mask_arch_inner_wss),20,data2.wss_m(geo_mask_arch_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS after interpolation: outer arch')
            scatter3(geo.x_coor_wss(geo_mask_arch_outer_wss),geo.y_coor_wss(geo_mask_arch_outer_wss),geo.z_coor_wss(geo_mask_arch_outer_wss),20,data2.wss_m(geo_mask_arch_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS after interpolation: inner DAo')
            scatter3(geo.x_coor_wss(geo_mask_DAo_inner_wss),geo.y_coor_wss(geo_mask_DAo_inner_wss),geo.z_coor_wss(geo_mask_DAo_inner_wss),20,data2.wss_m(geo_mask_DAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','WSS after interpolation: outer DAo')
            scatter3(geo.x_coor_wss(geo_mask_DAo_outer_wss),geo.y_coor_wss(geo_mask_DAo_outer_wss),geo.z_coor_wss(geo_mask_DAo_outer_wss),20,data2.wss_m(geo_mask_DAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        end
        
        mean_vel_after_interpolation(n,1) = mean_vel_asc_inner;
        mean_vel_after_interpolation(n,2) = mean_vel_asc_outer;
        mean_vel_after_interpolation(n,3) = mean_vel_arch_inner;
        mean_vel_after_interpolation(n,4) = mean_vel_arch_outer;
        mean_vel_after_interpolation(n,5) = mean_vel_DAo_inner;
        mean_vel_after_interpolation(n,6) = mean_vel_DAo_outer
        mean_wss_after_interpolation(n,1) = mean_wss_asc_inner;
        mean_wss_after_interpolation(n,2) = mean_wss_asc_outer;
        mean_wss_after_interpolation(n,3) = mean_wss_arch_inner;
        mean_wss_after_interpolation(n,4) = mean_wss_arch_outer;
        mean_wss_after_interpolation(n,5) = mean_wss_DAo_inner;
        mean_wss_after_interpolation(n,6) = mean_wss_DAo_outer
        
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
        IE_total_vel = abs(mean_vel_total_before_interpolation(n,1)-mean_vel_total_after_interpolation(n,1)) / ((mean_vel_total_before_interpolation(n,1)+mean_vel_total_after_interpolation(n,1))./2)*100;
        IE_total_wss = abs(mean_wss_total_before_interpolation(n,1)-mean_wss_total_after_interpolation(n,1)) / ((mean_wss_total_before_interpolation(n,1)+mean_wss_total_after_interpolation(n,1))./2)*100;
        
        disp(['IE velocity inner AAo = ' num2str(round(IE_inner_AAo_vel)) ' %'])
        disp(['IE velocity outer AAo = ' num2str(round(IE_outer_AAo_vel)) ' %'])
        disp(['IE velocity inner asc = ' num2str(round(IE_inner_asc_vel)) ' %'])
        disp(['IE velocity outer asc = ' num2str(round(IE_outer_asc_vel)) ' %'])
        disp(['IE velocity inner DAo = ' num2str(round(IE_inner_DAo_vel)) ' %'])
        disp(['IE velocity outer DAo = ' num2str(round(IE_outer_DAo_vel)) ' %'])
        disp(['IE velocity total = ' num2str(round(IE_total_vel)) ' %'])
        disp(' ')
        disp(['IE wall shear stress inner AAo = ' num2str(round(IE_inner_AAo_wss)) ' %'])
        disp(['IE wall shear stress outer AAo = ' num2str(round(IE_outer_AAo_wss)) ' %'])
        disp(['IE wall shear stress inner asc = ' num2str(round(IE_inner_asc_wss)) ' %'])
        disp(['IE wall shear stress outer asc = ' num2str(round(IE_outer_asc_wss)) ' %'])
        disp(['IE wall shear stress inner DAo = ' num2str(round(IE_inner_DAo_wss)) ' %'])
        disp(['IE wall shear stress outer DAo = ' num2str(round(IE_outer_DAo_wss)) ' %'])
        disp(['IE WSS total = ' num2str(round(IE_total_wss)) ' %'])
    end  
    close all
    
    VELx(:,n) = x_value_vel;
    VELy(:,n) = y_value_vel;
    VELz(:,n) = z_value_vel;
    WSSx(:,n) = x_value_wss;
    WSSy(:,n) = y_value_wss;
    WSSz(:,n) = z_value_wss;
    
    clear MrstructPath
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

% Put everything you need for further processing into atlas struct
atlas.x_coor_vel = geo.x_coor_vel; % atlas velocity coordinates
atlas.y_coor_vel = geo.y_coor_vel;
atlas.z_coor_vel = geo.z_coor_vel;
atlas.meanx_vel = sum(VELx,2)./size(PATHNAME,2); % velocity averaged over all subjects = atlas velocity
atlas.meany_vel = sum(VELy,2)./size(PATHNAME,2);
atlas.meanz_vel = sum(VELz,2)./size(PATHNAME,2);
atlas.x_coor_wss = geo.x_coor_wss; % atlas WSS coordinates
atlas.y_coor_wss = geo.y_coor_wss;
atlas.z_coor_wss = geo.z_coor_wss;
atlas.meanx_wss = sum(WSSx,2)./size(PATHNAME,2); % WSS aceraged over all subjects = atlas WSS
atlas.meany_wss = sum(WSSy,2)./size(PATHNAME,2);
atlas.meanz_wss = sum(WSSz,2)./size(PATHNAME,2);
atlas.vox = mask1_vox; % voxel size atlas
atlas.faces = geo.F; % faces of atlas
atlas.vertices = geo.V; % vertices of atlas (identical to atlas.coor_wss so not necessarily needed)
atlas.mask = mask1; % the atlas coordinates in matrix form
atlas.mean_vel = sqrt(atlas.meanx_vel.^2 + atlas.meany_vel.^2 + atlas.meanz_vel.^2); % atlas velocity magnitude
atlas.mean_wss = sqrt(atlas.meanx_wss.^2 + atlas.meany_wss.^2 + atlas.meanz_wss.^2); % atlas WSS magnitude

% calculate the standard deviation atlas map for velocity and WSS
for n = 1:size(PATHNAME,2)
    std_x_vel(:,n) = (VELx(:,n) - atlas.meanx_vel).^2;
    std_y_vel(:,n) = (VELy(:,n) - atlas.meany_vel).^2;
    std_z_vel(:,n) = (VELz(:,n) - atlas.meanz_vel).^2;
    std_x_wss(:,n) = (WSSx(:,n) - atlas.meanx_wss).^2;
    std_y_wss(:,n) = (WSSy(:,n) - atlas.meany_wss).^2;
    std_z_wss(:,n) = (WSSz(:,n) - atlas.meanz_wss).^2;
end
atlas.stdx_vel = sqrt(sum(std_x_vel,2)./size(PATHNAME,2));
atlas.stdy_vel = sqrt(sum(std_y_vel,2)./size(PATHNAME,2));
atlas.stdz_vel = sqrt(sum(std_z_vel,2)./size(PATHNAME,2));
atlas.std_vel = sqrt(atlas.stdx_vel.^2 + atlas.stdy_vel.^2 + atlas.stdz_vel.^2); % velocity standard deviation magnitude
atlas.stdx_wss = sqrt(sum(std_x_wss,2)./size(PATHNAME,2));
atlas.stdy_wss = sqrt(sum(std_y_wss,2)./size(PATHNAME,2));
atlas.stdz_wss = sqrt(sum(std_z_wss,2)./size(PATHNAME,2));
atlas.std_wss = sqrt(atlas.stdx_wss.^2 + atlas.stdy_wss.^2 + atlas.stdz_wss.^2); % WSS standard deviation magnitude

%if plotFlag == 1
figure('Name','Mean atlas velocity')
atlas_matrix = zeros(size(atlas.mask));
L = (atlas.mask~=0);
atlas_matrix(L) = atlas.mean_vel;
L_figure = (squeeze(max(atlas_matrix,[],3))~=0);
imagesc(squeeze(max(atlas_matrix,[],3)),'Alphadata',double(L_figure));
colorbar;axis tight; axis equal; axis ij; axis off;caxis([0 1.5]);%view([180 -90])

figure('Name','SD atlas velocity')
atlas_matrix = zeros(size(atlas.mask));
L = (atlas.mask~=0);
[I,J] = find(L==1);
atlas_matrix(L) = atlas.std_vel;
L_figure = (squeeze(max(atlas_matrix,[],3))~=0);
imagesc(squeeze(max(atlas_matrix,[],3)),'Alphadata',double(L_figure));
colorbar;axis tight; axis equal; axis ij; axis off;caxis([0 1]);%view([180 -90])

figure('Name','Mean atlas WSS')
patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none','FaceVertexCData',atlas.mean_wss,'FaceColor','interp','FaceAlpha',1);
colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])

figure('Name','SD atlas wss')
patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none', 'FaceVertexCData',atlas.std_wss,'FaceColor','interp','FaceAlpha',1);colorbar;
axis equal;axis off; axis ij;caxis([0 1.5]);view([-180 -90])
%end

if calculateRE_Flag == 1;
    disp(['Mean registration error = ' num2str(mean(diff_percentage(:,1))) ' +/- ' num2str(std(diff_percentage(:,1))) '%'])
end

if calculateIE_Flag == 1;
    error_matrix_vel = abs(mean_vel_before_interpolation-mean_vel_after_interpolation) ./ ...
        ((mean_vel_before_interpolation+mean_vel_after_interpolation)./2);
    error_matrix_wss = abs(mean_wss_before_interpolation-mean_wss_after_interpolation) ./ ...
        ((mean_wss_before_interpolation+mean_wss_after_interpolation)./2);
    error_matrix_vel_total = abs( mean_vel_total_before_interpolation- mean_vel_total_after_interpolation) ./ ...
        ((mean_vel_total_before_interpolation+ mean_vel_total_after_interpolation)./2);
    error_matrix_wss_total = abs(mean_wss_total_before_interpolation-mean_wss_total_after_interpolation) ./ ...
        ((mean_wss_total_before_interpolation+mean_wss_total_after_interpolation)./2);
    disp(['Velocity: Mean interpolation error inner AAo = ' num2str(mean(error_matrix_vel(:,1),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,1),1)*100) '%'])
    disp(['Velocity: Mean interpolation error outer AAo = ' num2str(mean(error_matrix_vel(:,2),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,2),1)*100) '%'])
    disp(['Velocity: Mean interpolation error inner arch = ' num2str(mean(error_matrix_vel(:,3),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,3),1)*100) '%'])
    disp(['Velocity: Mean interpolation error outer arch = ' num2str(mean(error_matrix_vel(:,4),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,4),1)*100) '%'])
    disp(['Velocity: Mean interpolation error inner DAo = ' num2str(mean(error_matrix_vel(:,5),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,5),1)*100) '%'])
    disp(['Velocity: Mean interpolation error outer DAo = ' num2str(mean(error_matrix_vel(:,6),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,6),1)*100) '%'])
    disp(['Velocity: averaged over the six regions = ' num2str(((mean(error_matrix_vel(:,1),1)+mean(error_matrix_vel(:,2),1)+mean(error_matrix_vel(:,3),1)+mean(error_matrix_vel(:,4),1)+mean(error_matrix_vel(:,5),1)+mean(error_matrix_vel(:,6),1))/6)*100) '%'])
    disp(['Velocity: Mean interpolation error total = ' num2str(mean(error_matrix_vel_total)*100) ' +/- ' num2str(std(error_matrix_vel_total)*100) '%'])
    disp(' ')
    disp(['WSS: Mean interpolation error inner AAo = ' num2str(mean(error_matrix_wss(:,1),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,1),1)*100) '%'])
    disp(['WSS: Mean interpolation error outer AAo = ' num2str(mean(error_matrix_wss(:,2),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,2),1)*100) '%'])
    disp(['WSS: Mean interpolation error inner arch = ' num2str(mean(error_matrix_wss(:,3),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,3),1)*100) '%'])
    disp(['WSS: Mean interpolation error outer arch = ' num2str(mean(error_matrix_wss(:,4),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,4),1)*100) '%'])
    disp(['WSS: Mean interpolation error inner DAo = ' num2str(mean(error_matrix_wss(:,5),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,5),1)*100) '%'])
    disp(['WSS: Mean interpolation error outer DAo = ' num2str(mean(error_matrix_wss(:,6),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,6),1)*100) '%'])
    disp(['Velocity: averaged over the six regions = ' num2str(((mean(error_matrix_wss(:,1),1)+mean(error_matrix_wss(:,2),1)+mean(error_matrix_wss(:,3),1)+mean(error_matrix_wss(:,4),1)+mean(error_matrix_wss(:,5),1)+mean(error_matrix_wss(:,6),1))/6)*100) '%'])
    disp(['WSS: Mean interpolation error total = ' num2str(mean(error_matrix_wss_total)*100) ' +/- ' num2str(std(error_matrix_wss_total)*100) '%'])
    
end

%end
