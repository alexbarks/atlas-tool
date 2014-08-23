function [probability_mask] = make_geometry_point_cloud(PATHNAME,plotFlag,saveFlag)

%%% [probability_mask] = make_geometry_point_cloud(PATHNAME,plotFlag,saveFlag)
%
% This function creates the probability mask (or idealized geometry) that will be used in the function
% 'make_atlas_point_cloud_scalars_affine_registration' of the batch data put into this file.
% The method developed was published in (so please cite this paper when using this method):
%
% A Methodology to Detect Abnormal Relative Wall Shear Stress on the Full Surface of the Thoracic Aorta Using 4D Flow MRI
% van Ooij P, Potters WV, Nederveen AJ, Allen BD, Collins J, Carr J, Malaisrie SC, Markl M, Barker AJ
% Magn Reson Med. 2014 Feb 25, doi: 10.1002/mrm.25224. [Epub ahead of print]
%
% First, the first two aortas (you can probably use this for carotids and intracranial vessels as well, but I haven't tried yet) are co-registered
% (rigid registration) and added to create an initial 'overlap map'. The third aorta is co-registered to the overlap map and added and so forth for
% all aortas in the cohort. Next, the threshold (O_thresh in the above mentioned paper) is varied from 1 to the number of aortas and again, all aortas
% are co-registered to the overlap map of varying thresholds and the registration error (RE) is determined. The overlap map where the RE averaged
% over all aortas is smallest is chosen as the probability mask (or idealized geometry).
%
% 2014, Pim van Ooij, Northwestern University
%
% Input
% 1)PATHNAME          : The PATHNAMES for the subjects that the probability_mask will be made of
% 2)plotFlag          : If plotFlag switched on, Matlab will output any possible image to the screen to check if everything happens correctly.
% 3)saveFlag          : You can choose to save the probability_mask to a directory of choice, but if make_geometry_point_cloud combined with 
%                       make_atlas_point_cloud_scalars_affine_registration there is really no need
% Output
% 1)probability_mask: The mask of the probability mask (or idealized aortic geometry)
%
% Usage
% This code is for creating the probability mask (or idealized aortic geometry) from a cohort that consists of 'data_done' structs as
% created by Pims_postprocessing tool
% The function for creating prbability masks from mrStructs is under construction
%
% Examples:
% [probability_mask] = make_geometry_point_cloud_NEW(PATHNAME,plotFlag,saveFlag)
% [probability_mask] = make_geometry_point_cloud('',0,1)
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

% %%% masks to load (use for debug)
%
if ~exist(PATHNAME{1}) == 2 || isempty(PATHNAME{1});
    % In for-loop %error('No Input PATHNAME given')
end

if nargin < 2 || isempty(plotFlag);
    plotFlag = 1;
end

FILENAME = 'mask_struct_aorta';

disp(['...Busy loading data_done aorta ' num2str(1)])
tic
load(strcat(PATHNAME{1},'\mrstruct\',FILENAME));
toc
disp(['Done loading data_done aorta'  num2str(1)]);disp(' ')
%data1 = data; clear data;
mask1 = mrstruct_mask.dataAy;
mask1_vox = mrstruct_mask.vox;clear mrstruct_mask
L = (mask1 ~= 0);

for n = 2:size(PATHNAME,2)
    
    [x,y,z] = meshgrid((1:size(mask1,2)).* mask1_vox(2), ...
        (1:size(mask1,1)).* mask1_vox(1),(1:size(mask1,3)).* mask1_vox(3));
    x_coor1 = x(L);y_coor1 = y(L);z_coor1 = z(L);
    clear x, clear y, clear z
    
    disp(['...Busy loading data_done aorta ' num2str(n)])
    tic
    if nargin < 1 || isempty(PATHNAME)
        [FILENAME,PATHNAME{n}] = uigetfile('.mat','Load probability mask');
        load(strcat(PATHNAME{n},FILENAME));
    else
        load(strcat(PATHNAME{n},'\mrstruct\',FILENAME));
    end
    toc
    disp(['Done loading data_done aorta '  num2str(n)])
    mask2 = mrstruct_mask.dataAy;
    mask2_vox = mrstruct_mask.vox;clear mrstruct_mask
        
    L = (mask2 ~= 0);
    [x,y,z] = meshgrid((1:size(mask2,2)).* mask2_vox(2), ...
        (1:size(mask2,1)).* mask2_vox(1),(1:size(mask2,3)).* mask2_vox(3));
    x_coor2 = x(L);y_coor2 = y(L);z_coor2 = z(L);
    clear x, clear y, clear z
    
    if plotFlag == 1
        figure('Name',strcat('To be registered aorta ',num2str(n)))
        plot3(x_coor1,y_coor1,z_coor1,'r.')
        hold on
        plot3(x_coor2,y_coor2,z_coor2,'b.')
        legend('to remain the same','to be transformed')
        axis equal; axis off;view([0 90]); axis ij
    end
    
    %%% Velocities
    PSF = fspecial('gaussian',1,1);
    mask1_to_register = mask1;
    mask1_to_register = imfilter(mask1_to_register,PSF,'conv');
    
    mask2 = imfilter(mask2,PSF,'conv');
    
    disp(' ')
    disp('...Busy registering...Rigid registration (dof = 6)')
    disp(' ')
    disp('...This can take up to 5 minutes...')
    
    tic
    % directory with flirt.exe and cygwin1.dll as read from atlas_tool.cfg in C:\temp
    fsldir = path_flirt;   
    
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
        'flirt -searchrx -0 0 -searchry -0 0 -searchrz -0 0 -dof 6' ...
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
    system(path_cygwin);
    
    % load transformation mask
    load Rotation_Translation -ascii
    
    flirtmat = 'Rotation_Translation.mat';
    src = infile;
    trg = reffile;
    [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
    
    %%% New coordinates
    yxz_coor = [y_coor2 x_coor2 z_coor2];clear x_coor2, clear y_coor2, clear z_coor2
    yxz_coor(:,4) = 1;
    yxz_coor_new = inv(worldmat)*yxz_coor'; clear yxz_coor
    x_coor = yxz_coor_new(2,:)';
    y_coor = yxz_coor_new(1,:)';
    z_coor = yxz_coor_new(3,:)'; clear yxz_coor_new
    toc
    
    if plotFlag == 1
        figure('Name',strcat('Aorta ',num2str(n),' registered '))
        plot3(x_coor1,y_coor1,z_coor1,'r.')
        hold on
        plot3(x_coor,y_coor,z_coor,'b.')
        legend('Prob mask',['Aorta ' num2str(n)])
        %legend('probability mask/atlas','individual aorta')
        axis equal;view([0 90]); axis ij
    end
    
    offset = 100;
    if n == 2
        %%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
        sizes = [size(mask1,1)+offset  size(mask1,2)+offset size(mask1,3)+offset];
        mask1b = zeros(sizes);
        mask1b((offset+1):size(mask1b,1),(offset+1):size(mask1b,2),(offset+1):size(mask1b,3)) = mask1;
        mask1c = mask1b;clear mask1b
    else
        %%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
        sizes = [size(probability_mask.matrix,1)+offset  size(probability_mask.matrix,2)+offset size(probability_mask.matrix,3)+offset];
        mask1b = zeros(sizes);
        mask1b((offset+1):size(mask1b,1),(offset+1):size(mask1b,2),(offset+1):size(mask1b,3)) = probability_mask.matrix;
        probability_mask.matrix = mask1b;clear mask1b
    end
    
    %%% Round the coordinates to be able to create a matrix
    x_round = round(x_coor./mask1_vox(1)) + offset;
    y_round = round(y_coor./mask1_vox(2)) + offset;
    z_round = round(z_coor./mask1_vox(3)) + offset;
    
    % Create matrix
    indices_mask = [x_round y_round z_round];
    siz=max(indices_mask,[],1);
    IND = sub2ind(siz,x_round,y_round,z_round);
    [b, IND_double_removed] = unique(IND);
    clear b
    indices_mask = [x_round(IND_double_removed) y_round(IND_double_removed) z_round(IND_double_removed)];
    clear IND_double_removed, clear IND
    mask_new = zeros([max(y_round) max(x_round) max(z_round)]);
    for i = 1:size(indices_mask,1)
        mask_new(indices_mask(i,2),indices_mask(i,1),indices_mask(i,3)) = 1;
    end
    
    % Due to the rounding of coordinates there are holes in the aorta, we fill them by smooth3
    % and then erode the aorta to give it the right size again
    se = strel('disk',1);
    mask_new = imerode(smooth3(mask_new),se);
    
    if n == 2
        % Make sure both masks have the same dimensions
        if size(mask1c,1) > size(mask_new,1)
            mask_new(size(mask_new,1):size(mask1c,1),:,:) = 0;
        elseif size(mask1c,1) < size(mask_new,1)
            %   mask1c(size(mask1c,1):size(mask_new,1),:,:) = 0;
            mask_new(size(mask1c,1)+1:size(mask_new,1),:,:) = [];
        end
        if size(mask1c,2) > size(mask_new,2)
            mask_new(:,size(mask_new,2):size(mask1c,2),:) = 0;
        elseif size(mask1c,2) < size(mask_new,2)
            %mask1c(:,size(mask1c,2):size(mask_new,2),:) = 0;
            mask_new(:,size(mask1c,2)+1:size(mask_new,2),:) = [];
        end
        if size(mask1c,3) > size(mask_new,3)
            mask_new(:,:,size(mask_new,3):size(mask1c,3)) = 0;
        elseif size(mask1c,3) < size(mask_new,3)
            %mask1c(:,:,size(mask1c,3):size(mask_new,3)) = 0;
            mask_new(:,:,size(mask1c,3)+1:size(mask_new,3)) = [];
        end
        
        L_mask_new = double(mask_new ~= 0);
        
        mask1c = mask1c+L_mask_new;
        probability_mask.matrix  = mask1c;
        
        probability_mask.matrix(1:offset,:,:) = [];
        probability_mask.matrix(:,1:offset,:) = [];
        probability_mask.matrix(:,:,1:offset) = [];
    else
        % Make sure both masks have the same dimensions
        if size(probability_mask.matrix,1) > size(mask_new,1)
            mask_new(size(mask_new,1):size(probability_mask.matrix,1),:,:) = 0;
        elseif size(probability_mask.matrix,1) < size(mask_new,1)
            %            probability_mask.matrix(size(probability_mask.matrix,1):size(mask_new,1),:,:) = 0;
            mask_new(size(probability_mask.matrix,1)+1:size(mask_new,1),:,:) = [];
        end
        if size(probability_mask.matrix,2) > size(mask_new,2)
            mask_new(:,size(mask_new,2):size(probability_mask.matrix,2),:) = 0;
        elseif size(probability_mask.matrix,2) < size(mask_new,2)
            mask_new(:,size(probability_mask.matrix,2)+1:size(mask_new,2),:) = [];
            %           probability_mask.matrix(:,size(probability_mask.matrix,2):size(mask_new,2),:) = 0;
        end
        if size(probability_mask.matrix,3) > size(mask_new,3)
            mask_new(:,:,size(mask_new,3):size(probability_mask.matrix,3)) = 0;
        elseif size(probability_mask.matrix,3) < size(mask_new,3)
            %           probability_mask.matrix(:,:,size(probability_mask.matrix,3):size(mask_new,3)) = 0;
            mask_new(:,:,size(probability_mask.matrix,3)+1:size(mask_new,3)) = [];
        end
        
        L_mask_new = double(mask_new ~= 0);
        
        probability_mask.matrix  = probability_mask.matrix + L_mask_new;
        
        probability_mask.matrix(1:offset,:,:) = [];
        probability_mask.matrix(:,1:offset,:) = [];
        probability_mask.matrix(:,:,1:offset) = [];
    end
    
    if plotFlag == 1
        figure('Name',num2str(n))
        L_figure = (squeeze(max(probability_mask.matrix,[],3))~=0);
        imagesc(squeeze(max(probability_mask.matrix,[],3)),'Alphadata',double(L_figure));
        colorbar; axis tight; axis equal;% axis off
    end
    
    L = (probability_mask.matrix~=0);
    mask1 = double(L);
    
    pause(8)
    close all
end

if plotFlag == 1
    figure('Name',num2str(n))
    L_figure = (squeeze(max(probability_mask.matrix,[],3))~=0);
    imagesc(squeeze(max(probability_mask.matrix,[],3)),'Alphadata',double(L_figure));
    colorbar;axis tight; axis equal; axis off;
    pause(8)
end

probability_mask.vox = mask1_vox;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The code to calculate the idealized geometry start here! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
diff_matrix = zeros(size(PATHNAME,2));

for n = 1:size(PATHNAME,2)
    
    disp(['...Busy loading data_done aorta ' num2str(n)])
    tic
    load(strcat(PATHNAME{n},'\mrstruct\',FILENAME))
    toc
    disp(['Done loading data_done aorta ' num2str(n)])
    mask2 = mrstruct_mask.dataAy;
    
    L = (mask2 ~= 0);
    [x,y,z] = meshgrid((1:size(mask2,2)).* mask2_vox(2), ...
        (1:size(mask2,1)).* mask2_vox(1),(1:size(mask2,3)).* mask2_vox(3));
    x_coor2 = x(L);y_coor2 = y(L);z_coor2 = z(L);
    clear x, clear y, clear z
    
    for threshold = 1:size(PATHNAME,2)
        
        L = (probability_mask.matrix >= threshold );
        [x,y,z] = meshgrid((1:size(probability_mask.matrix,2)).* probability_mask.vox(2), ...
            (1:size(probability_mask.matrix,1)).* probability_mask.vox(1),(1:size(probability_mask.matrix,3)).* probability_mask.vox(3));
        x_coor1 = x(L);y_coor1 = y(L);z_coor1 = z(L);
        
        mask1 = double(L);
        
        if plotFlag == 1
            figure('Name',num2str(threshold))
            L_figure = (squeeze(max(mask1,[],3))~=0);
            imagesc(squeeze(max(mask1,[],3)),'Alphadata',double(L_figure));
            colorbar;axis tight; axis equal; axis off
        end
        
        if plotFlag == 1
            figure('Name',strcat('To be registered aorta ',num2str(n)))
            plot3(x_coor1,y_coor1,z_coor1,'r.')
            hold on
            plot3(x_coor2,y_coor2,z_coor2,'b.')
            legend('to remain the same','to be transformed')
            axis equal; axis off;view([0 90]); axis ij
        end
        
        PSF = fspecial('gaussian',1,1);
        mask1_to_register = mask1;
        mask1_to_register = imfilter(mask1_to_register,PSF,'conv');
        
        mask2 = imfilter(mask2,PSF,'conv');
        
        disp(' ')
        disp('...Busy registering...Rigid registration (dof = 6)')
        disp(' ')
        disp('...This can take up to 5 minutes...')
        
        % directory with flirt.exe and cygwin1.dll as read from atlas_tool.cfg in C:\temp
        fsldir = path_flirt; 
        
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
            'flirt -searchrx -0 0 -searchry -0 0 -searchrz -0 0 -dof 6' ...
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
        system(path_cygwin);
        
        % load transformation mask
        load Rotation_Translation -ascii
        
        flirtmat = 'Rotation_Translation.mat';
        src = infile;
        trg = reffile;
        [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
        rotmat = worldmat(1:3,1:3);
        
        %%% New Coordinates
        yxz_coor = [y_coor2 x_coor2 z_coor2];%clear x_coor2, clear y_coor2, clear z_coor2
        yxz_coor(:,4) = 1;
        yxz_coor_new = inv(worldmat)*yxz_coor'; clear yxz_coor
        x_coor = yxz_coor_new(2,:)';
        y_coor = yxz_coor_new(1,:)';
        z_coor = yxz_coor_new(3,:)'; clear yxz_coor_new
        
        if plotFlag == 1
            figure('Name',strcat('Aorta ',num2str(n),' registered '))
            plot3(x_coor1,y_coor1,z_coor1,'r.')
            hold on
            plot3(x_coor,y_coor,z_coor,'b.')
            legend('Prob mask',['Aorta ' num2str(n)])
            %legend('probability mask/atlas','individual aorta')
            axis equal;view([0 90]); axis ij
            
            pause(10)
        end
        
        %%% Round the coordinates to be able to create a matrix
        x_round = round(x_coor./mask1_vox(1)) + offset;
        y_round = round(y_coor./mask1_vox(2)) + offset;
        z_round = round(z_coor./mask1_vox(3)) + offset;
        
        % Create matrix
        indices_mask = [x_round y_round z_round];
        siz=max(indices_mask,[],1);
        IND = sub2ind(siz,x_round,y_round,z_round);
        [b, IND_double_removed] = unique(IND);
        clear b
        indices_mask = [x_round(IND_double_removed) y_round(IND_double_removed) z_round(IND_double_removed)];
        clear IND_double_removed, clear IND
        mask_new = zeros([max(y_round) max(x_round) max(z_round)]);
        for i = 1:size(indices_mask,1)
            mask_new(indices_mask(i,2),indices_mask(i,1),indices_mask(i,3)) = 1;
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
        
        L_mask_new = double(mask_new ~= 0);
        
        if plotFlag == 1
            figure('Name',num2str(threshold))
            L_figure = (squeeze(max(mask1+L_mask_new,[],3))~=0);
            imagesc(squeeze(max(mask1+L_mask_new,[],3)),'Alphadata',double(L_figure));
            colorbar;axis tight; axis equal;
            pause(10)
        end
        
        difference = abs(mask1-L_mask_new);
        [I1,J] = find(mask1~=0);
        [I2,J] = find(L_mask_new~=0);
        mean_I = (size(I1,1) + size(I2,1))/2;
        [I_diff,J] = find(difference~=0);
        diff_voxels = size(I_diff,1);
        diff_percentage = ((diff_voxels / mean_I) * 100)./2;
        
        diff_matrix(n,threshold) = diff_percentage
        disp(['Difference between aorta and atlas = ' num2str(diff_percentage)])
        
        close all
    end
end

diff = squeeze(mean(diff_matrix,1));
[I,J] = find(diff==min(diff))
L = (probability_mask.matrix >= J);
mask1 = double(L);

disp(['The minimal error = ' num2str(diff(1,J)) '+/- ' num2str(std(diff_matrix(:,J))) ' %'])

if plotFlag == 1
    figure('Name',num2str(J))
    L_figure = (squeeze(max(mask1,[],3))~=0);
    imagesc(squeeze(max(mask1,[],3)),'Alphadata',double(L_figure));
    colorbar
    axis tight; axis equal; axis off
    pause(10)
end
toc

probability_mask.matrix = mask1;

if saveFlag
    directory = uigetdir('C:\1_Chicago\Data\MIMICS\');
    disp('...saving probablitity_mask...')
    save(strcat(directory,'\probability_mask'),'probability_mask')
    disp('done')
end

%end