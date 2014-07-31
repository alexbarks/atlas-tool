function make_geometry_point_cloud

% %%% masks to load
% 
PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\1_CAMRI-JEHAM\20140307_084523_Aera_NMH\3dpc\';
PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\2_CAMRI-LAJAM\20140313_102146_Aera_NMH\3dpc_nav80\';
PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\3_CAMRI-KIMCO\20140313_115319_Aera_NMH\3dpc_nav80\';
PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\4_CAMRI-THZUM\20140314_073534_Aera_NMH\3dpc_nav80\';
PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\5_CAMRI-MAGIM2\20140314_090452_Aera_NMH\3dpc_nav80\';
PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\6_CAMRI-WIANM\20140317_095816_Aera_NMH\3dpc_nav80\';
PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\7_CAMRI_ARNIF\20140404_082635_Aera_NMH\';
PATHNAME{8}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\8_CAMRI-FREPAM\20140404_094040_Aera_NMH\3dpc_nav80\';
PATHNAME{9}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\9_CAMRI_JUVIF\20140408_083611_Aera_NMH\3dpc_nav80\';
PATHNAME{10}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\10_CAMRI_JAKIM\20140414_085925_Aera_NMH\';
PATHNAME{11}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\11_CAMRI_VIVIM\20140416_081355_Aera_NMH\';
PATHNAME{12}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\12_CAMRI_RIJO\20140429_080804_Aera_NMH\';
PATHNAME{13}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\13_CAMRI_MASNF\20140502_075511_Aera_NMH\';
PATHNAME{14}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\14_CAMRI_ANMAM\20140506_081402_Aera_NMH\';
FILENAME = 'data_done';

load(strcat(PATHNAME{1},FILENAME))
data1 = data; clear data;
mask1 = squeeze(data1.PC_zeros(:,:,:,1,1)~=0);

%%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
sizes = [size(mask1,1)+40 size(mask1,2)+40 size(mask1,3)+40 size(mask1,4) size(mask1,5)];
mask1b = zeros(sizes);
mask1b(41:size(mask1b,1),41:size(mask1b,2),41:size(mask1b,3),:,:) = mask1;
mask1 = mask1b;clear mask1b
L = (mask1 ~= 0);

for n = 2:size(PATHNAME,2)
    n
    
    [x,y,z] = meshgrid((1:size(mask1,2)).* data1.vox(2), ...
        (1:size(mask1,1)).* data1.vox(1),(1:size(mask1,3)).* data1.vox(3));
    x_coor1 = x(L);y_coor1 = y(L);z_coor1 = z(L);
    clear x, clear y, clear z
    %     figure('Name','mask1')
    %     V = vol3d('cdata',L,'texture','3D','texturemap',L);
    %     xlabel('x');ylabel('y');zlabel('z')
    %     axis tight; axis equal;axis ij
    
    load(strcat(PATHNAME{n},FILENAME))
    data2 = data; clear data;
    mask2 = squeeze(data2.PC_zeros(:,:,:,1,1)~=0);
    
    %%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
    sizes = [size(mask2,1)+40 size(mask2,2)+40 size(mask2,3)+40 size(mask2,4) size(mask2,5)];
    mask2b = zeros(sizes);
    mask2b(41:size(mask2b,1),41:size(mask2b,2),41:size(mask2b,3),:,:) = mask2;
    mask2 = mask2b;clear mask2b
    
    L = (mask2 ~= 0);
    [x,y,z] = meshgrid((1:size(mask2,2)).* data2.vox(2), ...
        (1:size(mask2,1)).* data2.vox(1),(1:size(mask2,3)).* data2.vox(3));
    x_coor2 = x(L);y_coor2 = y(L);z_coor2 = z(L);
    clear x, clear y, clear z
    %     figure('Name','mask2')
    %     V = vol3d('cdata',L,'texture','3D','texturemap',L);
    %     xlabel('x');ylabel('y');zlabel('z')
    %     axis tight; axis equal;axis ij
    
    figure('Name',strcat('To be registered',num2str(n)))
    plot3(x_coor1,y_coor1,z_coor1,'r.')
    hold on
    plot3(x_coor2,y_coor2,z_coor2,'b.')
    legend('to remain the same','to be transformed')
    axis equal; axis off;view([0 90]); axis ij
    
    %%% Velocities
    PSF = fspecial('gaussian',1,1);
    mask1_to_register = mask1;
    mask1_to_register = imfilter(mask1_to_register,PSF,'conv');
    
    mask2 = imfilter(mask2,PSF,'conv');
    
    disp(' ')
    disp('...Busy registering...')
    disp(' ')
    disp('...This can take up to 5 minutes...')
    
    % directory with flirt.exe and cygwin1.dll (use cygwin convention)
    fsldir = '/cygdrive/c/1_Chicago/WSSprojectWithAmsterdam/flirt/';
    
    % save as nifti
    cnii=make_nii(mask1_to_register,[data1.vox(1) data1.vox(2) data1.vox(3)]);
    save_nii(cnii,'mask1.nii');
    mnii=make_nii(mask2,[data2.vox(1) data2.vox(2) data2.vox(3)]);
    save_nii(mnii,'mask2.nii');
    
    % run flirt (needs cygwin installed)
    infile = 'mask2.nii';
    reffile = 'mask1.nii';
    outfile = 'rmask2.nii';
    omat = 'Rotation_Translation.mat';
    
    disp('rigid registration (dof = 6)')
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
    system('c:\cygwin\bin\bash runflirt.sh');
    
    % load transformation mask
    load Rotation_Translation -ascii
    
    flirtmat = 'Rotation_Translation.mat';
    src = infile;
    trg = reffile;
    [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
    rotmat = worldmat(1:3,1:3);
    
    %%% Velocity
    yxz_coor = [y_coor2 x_coor2 z_coor2];clear x_coor2, clear y_coor2, clear z_coor2
    yxz_coor(:,4) = 1;
    yxz_coor_new = inv(worldmat)*yxz_coor'; clear yxz_coor
    
    %%% Velocity
    x_coor = yxz_coor_new(2,:)';
    y_coor = yxz_coor_new(1,:)';
    z_coor = yxz_coor_new(3,:)';
    clear yxz_coor_new
    
    figure('Name','Registered')
    plot3(x_coor1,y_coor1,z_coor1,'r.')
    hold on
    plot3(x_coor,y_coor,z_coor,'b.')
    legend('Prob mask','Aorta 3')
    %legend('probability mask/atlas','individual aorta')
    axis equal;view([0 90]); axis ij
    pause
    
    %%% Round the coordinates to be able to create a matrix
    x_round = round(x_coor./data1.vox(1));
    y_round = round(y_coor./data1.vox(2));
    z_round = round(z_coor./data1.vox(3));
    
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
    
    % Make sure both masks have the same dimensions
    if size(mask1,1) > size(mask_new,1)
        mask_new(size(mask_new,1):size(mask1,1),:,:) = 0;
    elseif size(mask1,1) < size(mask_new,1)
        mask1(size(mask1,1):size(mask_new,1),:,:) = 0;
    end
    if size(mask1,2) > size(mask_new,2)
        mask_new(:,size(mask_new,2):size(mask1,2),:) = 0;
    elseif size(mask1,2) < size(mask_new,2)
        mask1(:,size(mask1,2):size(mask_new,2),:) = 0;
    end
    if size(mask1,3) > size(mask_new,3)
        mask_new(:,:,size(mask_new,3):size(mask1,3)) = 0;
    elseif size(mask1,3) < size(mask_new,3)
        mask1(:,:,size(mask1,3):size(mask_new,3)) = 0;
    end
    
    L_mask_new = double(mask_new ~= 0);
    
    if n == 2
        mask1 = mask1+L_mask_new;
        probability_mask.matrix  = mask1;
    else
        % Make sure both masks have the same dimensions
        if size(probability_mask.matrix,1) > size(mask_new,1)
            mask_new(size(mask_new,1):size(probability_mask.matrix,1),:,:) = 0;
        elseif size(probability_mask.matrix,1) < size(mask_new,1)
            probability_mask.matrix(size(probability_mask.matrix,1):size(mask_new,1),:,:) = 0;
        end
        if size(probability_mask.matrix,2) > size(mask_new,2)
            mask_new(:,size(mask_new,2):size(probability_mask.matrix,2),:) = 0;
        elseif size(probability_mask.matrix,2) < size(mask_new,2)
            probability_mask.matrix(:,size(probability_mask.matrix,2):size(mask_new,2),:) = 0;
        end
        if size(probability_mask.matrix,3) > size(mask_new,3)
            mask_new(:,:,size(mask_new,3):size(probability_mask.matrix,3)) = 0;
        elseif size(probability_mask.matrix,3) < size(mask_new,3)
            probability_mask.matrix(:,:,size(probability_mask.matrix,3):size(mask_new,3)) = 0;
        end
        
        probability_mask.matrix  = probability_mask.matrix + L_mask_new;
    end
    
    figure('Name',num2str(n))
    L_figure = (squeeze(max(probability_mask.matrix,[],3))~=0);
    imagesc(squeeze(max(probability_mask.matrix,[],3)),'Alphadata',double(L_figure));
    colorbar; axis tight; axis equal; axis off
    
    L = (probability_mask.matrix~=0);
    mask1 = double(L);
    
    pause(1)
    clc, close all
end

figure('Name',num2str(n))
L_figure = (squeeze(max(probability_mask.matrix,[],3))~=0);
imagesc(squeeze(max(probability_mask.matrix,[],3)),'Alphadata',double(L_figure));
colorbar
axis tight; axis equal; axis off

probability_mask.vox = data1.vox;
% save C:\1_Chicago\probability_mask
% load C:\1_Chicago\probability_mask

% %Quantify percentages of overlap for the different thresholds
% for threshold = 1:size(PATHNAME,2)
%     L1 = (probability_mask.matrix >= 4 );
%     L2 = (probability_mask.matrix >= threshold );
%     
%     [I1,J1] = find(L1==1);
%     [I2,J2] = find(L2==1);
%     
%     size(I1,1)
%     size(I2,1)
%     
%     percentage_overlap(threshold,1) = (1-(size(I1,1)-size(I2,1))./(size(I1,1)))*100
%     
% end

% figure('Name',num2str(n))
% L_figure = (squeeze(max(probability_mask.matrix,[],3))>=4);
% imagesc(squeeze(max(probability_mask.matrix,[],3)),'Alphadata',double(L_figure));
% colorbar
% axis tight; axis equal; axis off
% 
% pause(1)

tic
diff_matrix = zeros(size(PATHNAME,2));
for threshold = 1:size(PATHNAME,2)
    
    L = (probability_mask.matrix >= threshold );
    
    mask1 = double(L);
    
    %     figure('Name',num2str(threshold))
    %     L_figure = (squeeze(max(mask1,[],3))~=0);
    %     imagesc(squeeze(max(mask1,[],3)),'Alphadata',double(L_figure));
    %     colorbar
    %     axis tight; axis equal; axis off
    %
    %     [x,y,z] = meshgrid((1:size(mask1,2)) .* probability_mask.vox(2), ...
    %         (1:size(mask1,1)).* probability_mask.vox(1),(1:size(mask1,3)).* probability_mask.vox(3));
    %     x_coor1 = x(L);y_coor1 = y(L);z_coor1 = z(L);
    %     figure('Name','L')
    %     V = vol3d('cdata',L,'texture','3D','texturemap',L);
    %     axis tight; axis equal, view([0 -90])
    
    for n = 1:size(PATHNAME,2)
        
        
        load(strcat(PATHNAME{n},FILENAME))
        data2 = data; clear data;
        mask2 = squeeze(data2.PC_zeros(:,:,:,1,1)~=0);
        
%         strMrstructMask = [PATHNAME{n},maskName];
%         %end
%         % Read MRstructs
%         [MaskStruct] = mrstruct_read(strMrstructMask); % load magnitude mrstruct (saved with velomap_tool)
%         
%         mask2 = MaskStruct.dataAy;
        
%         eval(['mask2 = (data' num2str(n) '.PC_zeros(:,:,:,1,7)~=0)'])
        
        %%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
        sizes = [size(mask2,1)+40 size(mask2,2)+40 size(mask2,3)+40 size(mask2,4) size(mask2,5)];
        mask2b = zeros(sizes);
        mask2b(41:size(mask2b,1),41:size(mask2b,2),41:size(mask2b,3),:,:) = mask2;
        mask2 = mask2b;clear mask2b
        
        L = (mask2 ~= 0);
        [x,y,z] = meshgrid((1:size(mask2,2)).* data2.vox(2), ...
            (1:size(mask2,1)).* data2.vox(1),(1:size(mask2,3)).* data2.vox(3));
        x_coor2 = x(L);y_coor2 = y(L);z_coor2 = z(L);
        clear x, clear y, clear z
        %     figure('Name','mask2')
        %     V = vol3d('cdata',L,'texture','3D','texturemap',L);
        %     xlabel('x');ylabel('y');zlabel('z')
        %     axis tight; axis equal;axis ij
        
        %         figure('Name',strcat('To be registered',num2str(n)))
        %         plot3(x_coor1,y_coor1,z_coor1,'r.')
        %         hold on
        %         plot3(x_coor2,y_coor2,z_coor2,'b.')
        %         legend('to remain the same','to be transformed')
        %         axis equal; axis off;view([0 90]); axis ij
        
        %%% Velocities
        PSF = fspecial('gaussian',1,1);
        mask1_to_register = mask1;
        mask1_to_register = imfilter(mask1_to_register,PSF,'conv');
        
        mask2 = imfilter(mask2,PSF,'conv');
        
        disp(' ')
        disp('...Busy registering (RIGID)...')
        disp(' ')
        disp('...This can take up to 5 minutes...')
        
        % directory with flirt.exe and cygwin1.dll (use cygwin convention)
        fsldir = '/cygdrive/c/1_Chicago/WSSprojectWithAmsterdam/flirt/';
        
        % save as nifti
        cnii=make_nii(mask1_to_register,[data1.vox(1) data1.vox(2) data1.vox(3)]);
        save_nii(cnii,'mask1.nii');
        mnii=make_nii(mask2,[data2.vox(1) data2.vox(2) data2.vox(3)]);
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
        system('c:\cygwin\bin\bash runflirt.sh');
        
        % load transformation mask
        load Rotation_Translation -ascii
        
        flirtmat = 'Rotation_Translation.mat';
        src = infile;
        trg = reffile;
        [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
        rotmat = worldmat(1:3,1:3);
        
        %%% Velocity
        yxz_coor = [y_coor2 x_coor2 z_coor2];clear x_coor2, clear y_coor2, clear z_coor2
        yxz_coor(:,4) = 1;
        yxz_coor_new = inv(worldmat)*yxz_coor'; clear yxz_coor
        
        %%% Velocity
        x_coor = yxz_coor_new(2,:)';
        y_coor = yxz_coor_new(1,:)';
        z_coor = yxz_coor_new(3,:)';
        clear yxz_coor_new
        
        %         figure('Name','Registered')
        %         plot3(x_coor1,y_coor1,z_coor1,'r.')
        %         hold on
        %         plot3(x_coor,y_coor,z_coor,'b.')
        %         legend('Prob mask','Aorta 3')
        %         %legend('probability mask/atlas','individual aorta')
        %         axis equal; axis off;view([0 90]); axis ij
        
        %%% Round the coordinates to be able to create a matrix
        x_round = round(x_coor./data1.vox(1));
        y_round = round(y_coor./data1.vox(2));
        z_round = round(z_coor./data1.vox(3));
        
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
        
        % Make sure both masks have the same dimensions
        if size(mask1,1) > size(mask_new,1)
            mask_new(size(mask_new,1):size(mask1,1),:,:) = 0;
        elseif size(mask1,1) < size(mask_new,1)
            mask1(size(mask1,1):size(mask_new,1),:,:) = 0;
        end
        if size(mask1,2) > size(mask_new,2)
            mask_new(:,size(mask_new,2):size(mask1,2),:) = 0;
        elseif size(mask1,2) < size(mask_new,2)
            mask1(:,size(mask1,2):size(mask_new,2),:) = 0;
        end
        if size(mask1,3) > size(mask_new,3)
            mask_new(:,:,size(mask_new,3):size(mask1,3)) = 0;
        elseif size(mask1,3) < size(mask_new,3)
            mask1(:,:,size(mask1,3):size(mask_new,3)) = 0;
        end
        
        L_mask_new = double(mask_new ~= 0);
        
        %         figure('Name','difference')
        difference = abs(mask1-L_mask_new);
        %         V = vol3d('cdata',difference,'texture','3D','texturemap',(difference));colorbar
        %         axis tight; axis equal
        [I1,J] = find(mask1~=0);
        [I2,J] = find(L_mask_new~=0);
        mean_I = (size(I1,1) + size(I2,1))/2;
        [I_diff,J] = find(difference~=0);
        diff_voxels = size(I_diff,1);
        diff_percentage = ((diff_voxels / mean_I) * 100)./2;
        disp(['Difference between aorta and atlas = ' num2str(diff_percentage)])
        
        diff_matrix(n,threshold) = diff_percentage
        
        %     pause
    end
    
end

diff = squeeze(mean(diff_matrix,1))
[I,J] = find(diff==min(diff))
L = (probability_mask.matrix >= J);
mask1 = double(L);

disp(['The minimal error = ' num2str(diff(1,J)) '+/- ' num2str(std(diff_matrix(:,J))) ' %'])

figure('Name',num2str(J))
L_figure = (squeeze(max(mask1,[],3))~=0);
imagesc(squeeze(max(mask1,[],3)),'Alphadata',double(L_figure));
colorbar
axis tight; axis equal; axis off
toc

% Translate back to original position for further use
mask1(1:40,:,:) = [];
mask1(:,1:40,:) = [];
mask1(:,:,1:40) = [];

probability_mask.matrix = mask1;

directory = uigetdir('C:\1_Chicago\Data\MIMICS\');
disp('...saving probablitity_mask...')
save(strcat(directory,'\probability_mask'),'probability_mask')
disp('done')
end