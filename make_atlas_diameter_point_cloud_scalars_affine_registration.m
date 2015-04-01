function [atlas] = make_atlas_diameter_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,plotFlag,calculateIE_Flag)

%%% [atlas] = make_atlas_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,plotFlag,calculateRE_Flag,calculateIE_Flag,peak_systolicFlag)
%
% This function creates the diameter atlas (cohort-averaged diameter map is more politically correct)
% of the batch data put into this file (see below)
%
% First, each aorta (you can probably use this for carotids and intracranial vessels as well, but I haven't tried yet) is co-registered
% (affine registration) to the idealized aorta geometry created with the function 'make_geometry_point_cloud.m'. The 3D diameter
% values are subsequently interpolated (nearest neighbour interpolation) to the wall coordinates of the idealized geometry. Finally,
% an average and standard deviation (SD) of the input aortas is calculated resulting in the nanmean and SD diameter atlas
% that are saved to the directory of choice for further use.
%
% 2015, Pim van Ooij, Northwestern University, Academic Medical Center, Amsterdam
%
% Input
% 1)PATHNAME          : The PATHNAMES for the subjects that the atlas will consist of
% 2)PATHNAME_probability_mask  : The PATHNAME containing the probability mask (or idealized geometry) that is used for the atlas
% 3)plotFlag        : If plotFlag switched, Matlab will output any possible image to the screen to check if everything happens correctly.
% 4)calculateIE_Flag: When switched on the interpolation error (IE, see paper mentioned above) will be calculated. Note that ROIs are needed
%                     which can be drawn manually when switched on.
%
% Output
% 1)atlas:           The nanmean and SD velocity and WSS atlas are saved to this struct
%
% Usage
% This code is for creating diameter atlases from 'Diameter_point_cloud' structs as created by mimics_to_diameter.m
%
% Examples:
% [atlas] = make_atlas_diameter_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,plotFlag,calculateIE_Flag)
% [atlas] = make_atlas_diameter_point_cloud_scalars_affine_registration('','',0,0)
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

isempty(PATHNAME_probability_mask)
if nargin < 2 || isempty(PATHNAME_probability_mask)
    [FILENAME,PATHNAME_probability_mask] = uigetfile('.mat','Load probability mask');
    load(strcat(PATHNAME_probability_mask,FILENAME))
else
    load([PATHNAME_probability_mask '\probability_mask' ])
end

if nargin < 3 || isempty(plotFlag)
    plotFlag = 1;
end

if nargin < 4 || isempty(calculateIE_Flag)
    calculateIE_Flag = 0;
end

%%%% datasets to load
%%%% Copy paste the folders of choice here, add more or discard if needed
FILENAME1 = 'mask_struct_aorta';              % 1: Load mask
FILENAME2 = 'Diameter_point_cloud_aorta';     % 2: Load diameter

mask1 = probability_mask.matrix;
mask1_vox = probability_mask.vox;clear probability_mask

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

    geo.x_coor_diam = geo.V(:,1);
    geo.y_coor_diam = geo.V(:,2);
    geo.z_coor_diam = geo.V(:,3);   

% Allocate sizes
%Diameter = zeros(size(geo.V,1),size(PATHNAME,2));
for n = 1:size(PATHNAME,2)
    disp(['Aorta number ' num2str(n)])
    PATHNAME{n}
    
    if nargin < 1 || isempty(PATHNAME)
        [FILENAME,PATHNAME{n}] = uigetfile('.mat','Load mask_struct_aorta');
        load([PATHNAME{n} filesep 'mrstruct' filesep FILENAME1])
    else
        load([PATHNAME{n} filesep 'mrstruct' filesep FILENAME1])
    end
    
%     % create velocity coordinates    
    mask2 = mrstruct_mask.dataAy;
    mask2_vox = mrstruct_mask.vox;clear mrstruct_mask
%     L2 = (mask2 ~= 0);
%     contours = zeros(size(L2));
%     contours(L2==0) = -1;
%     contours(L2==1) = 1;
%     [F,V] = isosurface(contours,0); % make a surface from the detected contours
%     V = V .* (ones(size(V,1),1) * mask2_vox(1:3));
%     [data2.F,data2.V] = SmoothLaplacian(F,V,15); %laplacian smoothing for surface (Kevin Moerman)
%     clear F, clear V 
    
    load([PATHNAME{n} filesep 'mrstruct' filesep FILENAME2])
    Diameter = Diameter_point_cloud.Diameter;   
    data2.F = Diameter_point_cloud.faces;
    data2.V = Diameter_point_cloud.vertices;clear Diameter_point_cloud

    data2.x_coor_diam = data2.V(:,1);
    data2.y_coor_diam = data2.V(:,2);
    data2.z_coor_diam = data2.V(:,3);      
    
    if plotFlag == 1        
        figure('Name','Diameter')
        patch('Faces',data2.F,'Vertices',data2.V, ...
            'EdgeColor','none','FaceVertexCData',Diameter,'FaceColor','interp','FaceAlpha',1);
        colorbar;caxis([0 6]);axis equal;axis off; axis ij;view([-180 -90])
    end
    
    if calculateIE_Flag == 1;
        if ~exist(strcat(PATHNAME{n},'\diameter_ROIs\mask1.mat'),'file')
            
            mkdir(PATHNAME{n},'\diameter_ROIs')
            
            F1=figure('Name','Aorta shape: Paused after finishing a region so press space when finished!');
            patch('Faces',data2.F,'Vertices',data2.V,'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.25);
            view([-180 -90]);axis ij;axis equal;axis off
            
            for i = 1:3
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo);
                region = getPosition(polyAAo);
                
                %          disp('saving, pausing')
                save(strcat([PATHNAME{n} '\diameter_ROIs\mask' num2str(i)]),'region');
                pause
            end
            close(F1)
        end
        load(strcat(PATHNAME{n},'\diameter_ROIs\mask1'))
        data2_mask_AAo_inner_diam = inpolygon(data2.x_coor_diam, data2.y_coor_diam, region(:,1), region(:,2));
        nanmean_diam_asc_inner = nanmean(Diameter(data2_mask_AAo_inner_diam));        
        diam_asc_inner = Diameter(data2_mask_AAo_inner_diam);
        diam_asc_inner(isnan(diam_asc_inner)) = 0;
        sort_asc_inner_diam = sort(diam_asc_inner);
        lasca_diam = round(0.05*length( Diameter(data2_mask_AAo_inner_diam)));
        max_inner_asc_diam = mean(sort_asc_inner_diam((length(sort_asc_inner_diam)-lasca_diam):length(sort_asc_inner_diam)));        
        load(strcat(PATHNAME{n},'\diameter_ROIs\mask2'))
        data2_mask_AAo_outer_diam = inpolygon(data2.x_coor_diam, data2.y_coor_diam, region(:,1), region(:,2));
        nanmean_diam_asc_outer = nanmean(Diameter(data2_mask_AAo_outer_diam));
        diam_asc_outer = Diameter(data2_mask_AAo_outer_diam);
        diam_asc_outer(isnan(diam_asc_outer)) = 0;
        sort_asc_outer_diam = sort(diam_asc_outer);
        lascb_diam = round(0.05*length( Diameter(data2_mask_AAo_outer_diam)));
        max_outer_asc_diam = mean(sort_asc_outer_diam((length(sort_asc_outer_diam)-lascb_diam):length(sort_asc_outer_diam)));         
        load(strcat(PATHNAME{n},'\diameter_ROIs\mask3'))
        data2_mask_arch_inner_diam = inpolygon(data2.x_coor_diam, data2.y_coor_diam, region(:,1), region(:,2));
        nanmean_diam_arch_inner = nanmean(Diameter(data2_mask_arch_inner_diam));
        diam_arch_inner = Diameter(data2_mask_arch_inner_diam);
        diam_arch_inner(isnan(diam_arch_inner)) = 0;
        sort_arch_inner_diam = sort(diam_arch_inner);
        larcha_diam = round(0.05*length( Diameter(data2_mask_arch_inner_diam)));
        max_inner_arch_diam = mean(sort_arch_inner_diam((length(sort_arch_inner_diam)-larcha_diam):length(sort_arch_inner_diam)));         
%         load(strcat(PATHNAME{n},'\diameter_ROIs\mask4'))
%         data2_mask_arch_outer_diam = inpolygon(data2.x_coor_diam, data2.y_coor_diam, region(:,1), region(:,2));
%         nanmean_diam_arch_outer = nanmean(Diameter(data2_mask_arch_outer_diam));
%         load(strcat(PATHNAME{n},'\diameter_ROIs\mask5'))
%         data2_mask_DAo_inner_diam = inpolygon(data2.x_coor_diam, data2.y_coor_diam, region(:,1), region(:,2));
%         nanmean_diam_DAo_inner = nanmean(Diameter(data2_mask_DAo_inner_diam));
%         load(strcat(PATHNAME{n},'\diameter_ROIs\mask6'))
%         data2_mask_DAo_outer_diam = inpolygon(data2.x_coor_diam, data2.y_coor_diam, region(:,1), region(:,2));
%         nanmean_diam_DAo_outer = nanmean(Diameter(data2_mask_DAo_outer_diam));
        
        nanmean_diam_total_before_interpolation(n,1) = nanmean(Diameter)
               
        nanmean_diam_before_interpolation(n,1) = nanmean_diam_asc_inner;
        nanmean_diam_before_interpolation(n,2) = nanmean_diam_asc_outer;
        nanmean_diam_before_interpolation(n,3) = nanmean_diam_arch_inner
        
        max_diam_before_interpolation(n,1) = max_inner_asc_diam;
        max_diam_before_interpolation(n,2) = max_outer_asc_diam;
        max_diam_before_interpolation(n,3) = max_inner_arch_diam       
%         nanmean_diam_before_interpolation(n,4) = nanmean_diam_arch_outer;
%         nanmean_diam_before_interpolation(n,5) = nanmean_diam_DAo_inner;
%         nanmean_diam_before_interpolation(n,6) = nanmean_diam_DAo_outer
        
        if plotFlag == 1
            figure('Name','Diameter before interpolation: inner AAo')
            scatter3(data2.x_coor_diam(data2_mask_AAo_inner_diam),data2.y_coor_diam(data2_mask_AAo_inner_diam),data2.z_coor_diam(data2_mask_AAo_inner_diam),20,Diameter(data2_mask_AAo_inner_diam),'filled');
            axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','Diameter before interpolation: outer AAo')
            scatter3(data2.x_coor_diam(data2_mask_AAo_outer_diam),data2.y_coor_diam(data2_mask_AAo_outer_diam),data2.z_coor_diam(data2_mask_AAo_outer_diam),20,Diameter(data2_mask_AAo_outer_diam),'filled');
            axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure('Name','Diameter before interpolation: inner arch')
     
       scatter3(data2.x_coor_diam(data2_mask_arch_inner_diam),data2.y_coor_diam(data2_mask_arch_inner_diam),data2.z_coor_diam(data2_mask_arch_inner_diam),20,Diameter(data2_mask_arch_inner_diam),'filled');
            axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
%             figure('Name','Diameter before interpolation: outer arch')
%             scatter3(data2.x_coor_diam(data2_mask_arch_outer_diam),data2.y_coor_diam(data2_mask_arch_outer_diam),data2.z_coor_diam(data2_mask_arch_outer_diam),20,Diameter(data2_mask_arch_outer_diam),'filled');
%             axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
%             figure('Name','Diameter before interpolation: inner DAo')
%             scatter3(data2.x_coor_diam(data2_mask_DAo_inner_diam),data2.y_coor_diam(data2_mask_DAo_inner_diam),data2.z_coor_diam(data2_mask_DAo_inner_diam),20,Diameter(data2_mask_DAo_inner_diam),'filled');
%             axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
%             figure('Name','Diameter before interpolation: outer DAo')
%             scatter3(data2.x_coor_diam(data2_mask_DAo_outer_diam),data2.y_coor_diam(data2_mask_DAo_outer_diam),data2.z_coor_diam(data2_mask_DAo_outer_diam),20,Diameter(data2_mask_DAo_outer_diam),'filled');
%             axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
        end    
    end
    
    if plotFlag == 1
        figure('Name',strcat('To be registered aorta ',num2str(n)))
        plot3(geo.x_coor_diam,geo.y_coor_diam,geo.z_coor_diam,'r.')
        hold on
        %      plot3(geo.x_coor_wss,geo.y_coor_wss,geo.z_coor_wss,'g.')
        plot3(data2.x_coor_diam,data2.y_coor_diam,data2.z_coor_diam,'b.')
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
        
    %%% WSS coordinates
    yxz_coor_diam = [data2.y_coor_diam data2.x_coor_diam data2.z_coor_diam];
    yxz_coor_diam(:,4) = 1;
    yxz_coor_diam_new = inv(worldmat)*yxz_coor_diam'; clear yxz_coor_diam
    x_coor_diam = yxz_coor_diam_new(2,:)';
    y_coor_diam = yxz_coor_diam_new(1,:)';
    z_coor_diam = yxz_coor_diam_new(3,:)';clear yxz_coor_diam_new
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
        plot3(geo.x_coor_diam,geo.y_coor_diam,geo.z_coor_diam,'r.')
        hold on
        plot3(x_coor_diam,y_coor_diam,z_coor_diam,'b.')
        legend('Prob mask','Aorta')
        %legend('probability mask/atlas','individual aorta')
        axis equal; axis off;view([-180 -90]); axis ij
        pause(2)
        
        figure('Name','transformed Diameter')
        patch('Faces',data2.F,'Vertices',[x_coor_diam y_coor_diam z_coor_diam], ...
            'EdgeColor','none','FaceVertexCData',Diameter,'FaceColor','interp','FaceAlpha',1);
        colorbar;caxis([0 6]);axis equal;axis off; axis ij;view([-180 -90])
    end

    %%% Interpolate diameter to co-registered coordinates
    interpolation_function = TriScatteredInterp([x_coor_diam y_coor_diam z_coor_diam],Diameter,'nearest');
    diameter = interpolation_function([geo.x_coor_diam geo.y_coor_diam geo.z_coor_diam]);
    diameter(isnan(diameter)) = 0;

    if plotFlag == 1       
        figure('Name','interpolated to atlas Diameter')
        patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none','FaceVertexCData',diameter,'FaceColor','interp','FaceAlpha',1);
        colorbar;caxis([0 6]);axis equal;axis off; axis ij;view([-180 -90])
    end
    
    if calculateIE_Flag == 1;
        
        if ~exist(strcat(PATHNAME_probability_mask,'\diameter_ROIs\mask1.mat'),'file')
            
            F2=figure('Name','Atlas shape: Paused after finishing a region so press space when finished!');
            patch('Faces',geo.F,'Vertices',[geo.x_coor_diam geo.y_coor_diam geo.z_coor_diam],'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.25);
            view([-180 -90]);axis ij;axis equal;axis off
            
            for i = 1:3
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo);
                region = getPosition(polyAAo);
                
                %          disp('saving, pausing')
                mkdir(PATHNAME_probability_mask,'\diameter_ROIs')
                save(strcat([PATHNAME_probability_mask '\diameter_ROIs\mask' num2str(i)]),'region');
                pause
            end
            
            close(F2)
        end
        
        load(strcat(PATHNAME_probability_mask,'\diameter_ROIs\mask1'))
        geo_mask_AAo_inner_diam = inpolygon(geo.x_coor_diam, geo.y_coor_diam, region(:,1), region(:,2));
        nanmean_diam_asc_inner = nanmean(diameter(geo_mask_AAo_inner_diam));
        load(strcat(PATHNAME_probability_mask,'\diameter_ROIs\mask2'))
        geo_mask_AAo_outer_diam = inpolygon(geo.x_coor_diam, geo.y_coor_diam, region(:,1), region(:,2));
        nanmean_diam_asc_outer = nanmean(diameter(geo_mask_AAo_outer_diam));
        load(strcat(PATHNAME_probability_mask,'\diameter_ROIs\mask3'))
        geo_mask_arch_inner_diam = inpolygon(geo.x_coor_diam, geo.y_coor_diam, region(:,1), region(:,2));
        nanmean_diam_arch_inner = nanmean(diameter(geo_mask_arch_inner_diam));
%         load(strcat(PATHNAME_probability_mask,'\diameter_ROIs\mask4'))
%         geo_mask_arch_outer_diam = inpolygon(geo.x_coor_diam, geo.y_coor_diam, region(:,1), region(:,2));
%         nanmean_diam_arch_outer = nanmean(diameter(geo_mask_arch_outer_diam));
%         load(strcat(PATHNAME_probability_mask,'\diameter_ROIs\mask5'))
%         geo_mask_DAo_inner_diam = inpolygon(geo.x_coor_diam, geo.y_coor_diam, region(:,1), region(:,2));
%         nanmean_diam_DAo_inner = nanmean(diameter(geo_mask_DAo_inner_diam));
%         load(strcat(PATHNAME_probability_mask,'\diameter_ROIs\mask6'))
%         geo_mask_DAo_outer_diam = inpolygon(geo.x_coor_diam, geo.y_coor_diam, region(:,1), region(:,2));
%         nanmean_diam_DAo_outer = nanmean(diameter(geo_mask_DAo_outer_diam));

        nanmean_diam_total_after_interpolation(n,1) = nanmean(diameter)
        
        if plotFlag == 1
            figure('Name','Diameter after interpolation: inner AAo')
            scatter3(geo.x_coor_diam(geo_mask_AAo_inner_diam),geo.y_coor_diam(geo_mask_AAo_inner_diam),geo.z_coor_diam(geo_mask_AAo_inner_diam),20,diameter(geo_mask_AAo_inner_diam),'filled');
            axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Diameter after interpolation: outer AAo')
            scatter3(geo.x_coor_diam(geo_mask_AAo_outer_diam),geo.y_coor_diam(geo_mask_AAo_outer_diam),geo.z_coor_diam(geo_mask_AAo_outer_diam),20,diameter(geo_mask_AAo_outer_diam),'filled');
            axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure('Name','Diameter after interpolation: inner arch')
            scatter3(geo.x_coor_diam(geo_mask_arch_inner_diam),geo.y_coor_diam(geo_mask_arch_inner_diam),geo.z_coor_diam(geo_mask_arch_inner_diam),20,diameter(geo_mask_arch_inner_diam),'filled');
            axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
%             figure('Name','Diameter after interpolation: outer arch')
%             scatter3(geo.x_coor_diam(geo_mask_arch_outer_diam),geo.y_coor_diam(geo_mask_arch_outer_diam),geo.z_coor_diam(geo_mask_arch_outer_diam),20,diameter(geo_mask_arch_outer_diam),'filled');
%             axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
%             figure('Name','Diameter after interpolation: inner DAo')
%             scatter3(geo.x_coor_diam(geo_mask_DAo_inner_diam),geo.y_coor_diam(geo_mask_DAo_inner_diam),geo.z_coor_diam(geo_mask_DAo_inner_diam),20,diameter(geo_mask_DAo_inner_diam),'filled');
%             axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
%             figure('Name','Diameter after interpolation: outer DAo')
%             scatter3(geo.x_coor_diam(geo_mask_DAo_outer_diam),geo.y_coor_diam(geo_mask_DAo_outer_diam),geo.z_coor_diam(geo_mask_DAo_outer_diam),20,diameter(geo_mask_DAo_outer_diam),'filled');
%             axis equal;caxis([0 6]);xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        end

        nanmean_diam_after_interpolation(n,1) = nanmean_diam_asc_inner;
        nanmean_diam_after_interpolation(n,2) = nanmean_diam_asc_outer;
        nanmean_diam_after_interpolation(n,3) = nanmean_diam_arch_inner
%         nanmean_diam_after_interpolation(n,4) = nanmean_diam_arch_outer;
%         nanmean_diam_after_interpolation(n,5) = nanmean_diam_DAo_inner;
%         nanmean_diam_after_interpolation(n,6) = nanmean_diam_DAo_outer

        IE_inner_AAo_diam = abs(nanmean_diam_before_interpolation(n,1)-nanmean_diam_after_interpolation(n,1)) / ((nanmean_diam_before_interpolation(n,1)+nanmean_diam_after_interpolation(n,1))./2)*100;
        IE_outer_AAo_diam = abs(nanmean_diam_before_interpolation(n,2)-nanmean_diam_after_interpolation(n,2)) / ((nanmean_diam_before_interpolation(n,2)+nanmean_diam_after_interpolation(n,2))./2)*100;
        IE_inner_asc_diam = abs(nanmean_diam_before_interpolation(n,3)-nanmean_diam_after_interpolation(n,3)) / ((nanmean_diam_before_interpolation(n,3)+nanmean_diam_after_interpolation(n,3))./2)*100;
%         IE_outer_asc_diam = abs(nanmean_diam_before_interpolation(n,4)-nanmean_diam_after_interpolation(n,4)) / ((nanmean_diam_before_interpolation(n,4)+nanmean_diam_after_interpolation(n,4))./2)*100;
%         IE_inner_DAo_diam = abs(nanmean_diam_before_interpolation(n,5)-nanmean_diam_after_interpolation(n,5)) / ((nanmean_diam_before_interpolation(n,5)+nanmean_diam_after_interpolation(n,5))./2)*100;
%         IE_outer_DAo_diam = abs(nanmean_diam_before_interpolation(n,6)-nanmean_diam_after_interpolation(n,6)) / ((nanmean_diam_before_interpolation(n,6)+nanmean_diam_after_interpolation(n,6))./2)*100;
        IE_total_diam = abs(nanmean_diam_total_before_interpolation(n,1)-nanmean_diam_total_after_interpolation(n,1)) / ((nanmean_diam_total_before_interpolation(n,1)+nanmean_diam_total_after_interpolation(n,1))./2)*100;
        
        disp(['IE diameter inner AAo = ' num2str(round(IE_inner_AAo_diam)) ' %'])
        disp(['IE diameter outer AAo = ' num2str(round(IE_outer_AAo_diam)) ' %'])
        disp(['IE diameter inner asc = ' num2str(round(IE_inner_asc_diam)) ' %'])
%         disp(['IE diameter outer asc = ' num2str(round(IE_outer_asc_diam)) ' %'])
%         disp(['IE diameter inner DAo = ' num2str(round(IE_inner_DAo_diam)) ' %'])
%         disp(['IE diameter outer DAo = ' num2str(round(IE_outer_DAo_diam)) ' %'])
         disp(['IE Diameter total = ' num2str(round(IE_total_diam)) ' %'])
    end
    
    close all
    
    DIAMETER(:,n) = diameter;
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

atlas.x_coor_diam = geo.x_coor_diam; % atlas Diameter coordinates
atlas.y_coor_diam = geo.y_coor_diam;
atlas.z_coor_diam = geo.z_coor_diam;
atlas.nanmean_diameter = sum(DIAMETER,2)./size(PATHNAME,2); % Diameter aceraged over all subjects = atlas WSS
atlas.vox = mask1_vox; % voxel size atlas
atlas.faces = geo.F; % faces of atlas
atlas.vertices = geo.V; % vertices of atlas (identical to atlas.coor_diam so not necessarily needed)
atlas.mask = mask1; % the atlas coordinates in matrix form

% calculate the standard deviation atlas map for velocity and WSS
for n = 1:size(PATHNAME,2)
    std_diameter(:,n) = (DIAMETER(:,n) - atlas.nanmean_diameter).^2;
end
atlas.std_diameter = sqrt(sum(std_diameter,2)./size(PATHNAME,2));

figure('Name','Mean atlas diam')
patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none','FaceVertexCData',atlas.nanmean_diameter,'FaceColor','interp','FaceAlpha',1);
colorbar;caxis([0 5]);axis equal;axis off; axis ij;view([-180 -90])

figure('Name','SD atlas diam')
patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none', 'FaceVertexCData',atlas.std_diameter,'FaceColor','interp','FaceAlpha',1);colorbar;
axis equal;axis off; axis ij;caxis([0 1]);view([-180 -90])
%end

if calculateIE_Flag == 1;
    error_matrix_diam = abs(nanmean_diam_before_interpolation-nanmean_diam_after_interpolation) ./ ...
        ((nanmean_diam_before_interpolation+nanmean_diam_after_interpolation)./2);
    error_matrix_diam_total = abs(nanmean_diam_total_before_interpolation-nanmean_diam_total_after_interpolation) ./ ...
        ((nanmean_diam_total_before_interpolation+nanmean_diam_total_after_interpolation)./2);
    disp(['Diameter: Mean interpolation error inner AAo = ' num2str(nanmean(error_matrix_diam(:,1),1)*100) ' +/- ' num2str(std(error_matrix_diam(:,1),1)*100) '%'])
    disp(['Diameter: Mean interpolation error outer AAo = ' num2str(nanmean(error_matrix_diam(:,2),1)*100) ' +/- ' num2str(std(error_matrix_diam(:,2),1)*100) '%'])
    disp(['Diameter: Mean interpolation error inner arch = ' num2str(nanmean(error_matrix_diam(:,3),1)*100) ' +/- ' num2str(std(error_matrix_diam(:,3),1)*100) '%'])
%     disp(['Diameter: Mean interpolation error outer arch = ' num2str(nanmean(error_matrix_diam(:,4),1)*100) ' +/- ' num2str(std(error_matrix_diam(:,4),1)*100) '%'])
%     disp(['Diameter: Mean interpolation error inner DAo = ' num2str(nanmean(error_matrix_diam(:,5),1)*100) ' +/- ' num2str(std(error_matrix_diam(:,5),1)*100) '%'])
%     disp(['Diameter: Mean interpolation error outer DAo = ' num2str(nanmean(error_matrix_diam(:,6),1)*100) ' +/- ' num2str(std(error_matrix_diam(:,6),1)*100) '%'])
    disp(['Diameter: Mean interpolation error total = ' num2str(nanmean(error_matrix_diam_total)*100) ' +/- ' num2str(std(error_matrix_diam_total)*100) '%'])  
end

%end
