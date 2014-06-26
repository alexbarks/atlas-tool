function heat_map_figure

global data1
global data2

probability_mask = [];
data = [];
Rotation_Translation = [];
%clc, clear, close all

%[FILENAME, PATHNAME] = uigetfile('C:\1_Chicago\Data\MIMICS','Load probability_mask')
PATHNAME = 'c:\_ensightCases\bav_tissue\atlas_data\';
FILENAME = 'probability_mask_rigid'; %this is the atlas prob mask
load(strcat(PATHNAME,FILENAME))
data1 = probability_mask;
%L1 = (data1.matrix~=0);
L1 = (data1.matrix == 4 | data1.matrix == 5 | data1.matrix == 6 | data1.matrix == 7 | data1.matrix == 8 | data1.matrix == 9 | data1.matrix == 10);
[x,y,z] = meshgrid((1:size(data1.matrix,2)) .* data1.vox(2), ...
    (1:size(data1.matrix,1)).* data1.vox(1),(1:size(data1.matrix,3)).* data1.vox(3));
data1.x_coor_vel = x(L1);data1.y_coor_vel = y(L1);data1.z_coor_vel = z(L1);
% figure('Name','L1')
% V = vol3d('cdata',L1,'texture','3D','texturemap',L1);
% axis tight; axis equal
mask1 = double(L1);
contours = zeros(size(L1));
contours(L1==0) = -1;
contours(L1==1) = 1;
[data1.F,V] = isosurface(contours,0); % make a surface from the detected contours
data1.V = V .* (ones(size(V,1),1) * data1.vox(1:3));
%[data1.F,data1.V] = SmoothLaplacian(data1.F,data1.V,15); %laplacian smoothing for surface (Kevin Moerman)

load c:\_ensightCases\bav_tissue\atlas_data\atlas_SD % mean and std of the atlas

% figure('Name','Mean')
% patch('Faces',data1.F,'Vertices',data1.V,'EdgeColor','none', 'FaceVertexCData',atlas_SD.wssmean,'FaceColor','interp','FaceAlpha',1);colorbar;
% axis equal;axis off
% caxis([0 1.5])
% view([90 45])
% zoom(3)
%
% figure('Name','SD')
% patch('Faces',data1.F,'Vertices',data1.V,'EdgeColor','none', 'FaceVertexCData',atlas_SD.wssSD,'FaceColor','interp','FaceAlpha',1);colorbar;
% axis equal;axis off
% caxis([0 1.5])
% view([90 45])
% zoom(3)

[FILENAME, PATHNAME] = uigetfile('c:\_ensightCases\bav_tissue\','Load aorta data_done');
% PATHNAME = 'c:\_ensightCases\bav_tissue\PT279-MK\3dpc\mrstruct_subject_20130918_user022\'
% FILENAME = 'data_done'
load(strcat(PATHNAME,FILENAME))
%data.velocity_magnitude = squeeze(data.velocity_magnitude);
data2 = data; clear data
if ~isfield(data2,'PC_unaliased')
    data2.PC_unaliased = data2.PC_zeros;
end
L2 = (data2.PC_unaliased(:,:,:,1,1) ~= 0);
[x,y,z] = meshgrid((1:size(data2.PC_unaliased,2)).* data2.vox(2), ...
    (1:size(data2.PC_unaliased,1)).* data2.vox(1),(1:size(data2.PC_unaliased,3)).* data2.vox(3));
data2.x_coor_vel = x(L2);data2.y_coor_vel = y(L2);data2.z_coor_vel = z(L2);
clear x, clear y, clear z

data2.vox(3)

% figure('Name','mean_min_2SD_atlas ')
% patch('Faces',data2.F_matrix{1},'Vertices',[data2.V_matrix{1}(:,1)./data2.vox(1),data2.V_matrix{1}(:,2)./data2.vox(1),data2.V_matrix{1}(:,3)./data2.vox(1)],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);

%caxis([0 1.5])
%view([90 45])
%zoom(3)
%break

for t = 1:size(data2.PC_unaliased,5)
    vx = squeeze(data2.PC_unaliased(:,:,:,1,t));
    vy = squeeze(data2.PC_unaliased(:,:,:,2,t));
    vz = squeeze(data2.PC_unaliased(:,:,:,3,t));
    vmagn = sqrt(vx.^2 + vy.^2 + vz.^2);
    L =(vmagn ~= 0);
    mean_velo(t) = mean(vmagn(L));
end

% figure('Name','Mean velocity')
% plot(1:size(data2.PC_unaliased,5),mean_velo,'-ro','LineWidth',5,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',16);

[I,time2] = find(mean_velo==max(mean_velo));

%%% Wall shear stress coordinates for both datasets
data1.x_coor_wss = data1.V(:,1);
data1.y_coor_wss = data1.V(:,2);
data1.z_coor_wss = data1.V(:,3);

data2.x_coor_wss = data2.V_matrix{1}(:,1);
data2.y_coor_wss = data2.V_matrix{1}(:,2);
data2.z_coor_wss = data2.V_matrix{1}(:,3);

% Systolic timesteps averaged
data2.x_value_wss_t1 = data2.WSS_matrix{time2-2}(:,1);data2.y_value_wss_t1 = data2.WSS_matrix{time2-2}(:,2);data2.z_value_wss_t1 = data2.WSS_matrix{3}(:,3);
data2.x_value_wss_t2 = data2.WSS_matrix{time2-1}(:,1);data2.y_value_wss_t2 = data2.WSS_matrix{time2-1}(:,2);data2.z_value_wss_t2 = data2.WSS_matrix{4}(:,3);
data2.x_value_wss_t3 = data2.WSS_matrix{time2}(:,1);data2.y_value_wss_t3 = data2.WSS_matrix{time2}(:,2);data2.z_value_wss_t3 = data2.WSS_matrix{5}(:,3);
data2.x_value_wss_t4 = data2.WSS_matrix{time2+1}(:,1);data2.y_value_wss_t4 = data2.WSS_matrix{time2+1}(:,2);data2.z_value_wss_t4 = data2.WSS_matrix{6}(:,3);
data2.x_value_wss_t5 = data2.WSS_matrix{time2+2}(:,1);data2.y_value_wss_t5 = data2.WSS_matrix{time2+2}(:,2);data2.z_value_wss_t5 = data2.WSS_matrix{7}(:,3);

data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4 + data2.x_value_wss_t5)./5;
data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4 + data2.y_value_wss_t5)./5;
data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4 + data2.z_value_wss_t5)./5;

data2.wss_m = sqrt(data2.x_value_wss.^2 + data2.y_value_wss.^2 + data2.z_value_wss.^2);

figure('Name','To be registered')
plot3(data1.x_coor_vel,data1.y_coor_vel,data1.z_coor_vel,'r.')
hold on
plot3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,'b.')
legend('to remain the same','to be transformed')
axis equal; axis ij; axis off; view([-180 -90])  

%%% Velocities
PSF = fspecial('gaussian',5,1);
mask1_to_register = mask1;
mask1_to_register = imfilter(mask1_to_register,PSF,'conv');

mask2 = double(L2); clear L2
mask2 = imfilter(mask2,PSF,'conv');

% figure('Name','mask1')
% V = vol3d('cdata',mask1_to_register,'texture','3D','texturemap',mask1);
% axis tight; axis equal
%
% figure('Name','mask2')
% V = vol3d('cdata',mask2,'texture','3D','texturemap',mask2);
% axis tight; axis equal

disp(' ')
disp('...Busy registering...')
disp(' ')
disp('...This can take up to 5 minutes...')

% directory with flirt.exe and cygwin1.dll (use cygwin convention)
fsldir = '/cygdrive/d/research/matlabcode/matlab_registration/flirt/';

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

disp('rigid registration')
flirtcmd= [...
    'flirt -searchrx -5 5 -searchry -5 5 -searchrz -5 5 -dof 6' ...
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
worldmat
rotmat = worldmat(1:3,1:3);

%%% Velocity
yxz_coor_vel = [data2.y_coor_vel data2.x_coor_vel data2.z_coor_vel];
yxz_coor_vel(:,4) = 1;
yxz_coor_vel_new = inv(worldmat)*yxz_coor_vel'; clear YXZ_coor_vel

%%% WSS coordinates
yxz_coor_wss = [data2.y_coor_wss data2.x_coor_wss data2.z_coor_wss];
yxz_coor_wss(:,4) = 1;
yxz_coor_wss_new = inv(worldmat)*yxz_coor_wss'; clear YXZ_coor_wss

%%% WSS values
yxz_value_wss = [data2.y_value_wss data2.x_value_wss data2.z_value_wss];
yxz_value_wss_new = inv(rotmat)*yxz_value_wss'; clear YXZ_coor_wss

%%% Velocity
data2.x_coor_vel_new = yxz_coor_vel_new(2,:)';% - min(YXZ_mask2_new(1,:));
data2.y_coor_vel_new = yxz_coor_vel_new(1,:)';% - min(YXZ_mask2_new(2,:));
data2.z_coor_vel_new = yxz_coor_vel_new(3,:)';% - min(YXZ_mask2_new(3,:));
%clear yxz_coor_vel_new

%%% WSS coordinates
data2.x_coor_wss_new = yxz_coor_wss_new(2,:)';% - min(YXZ_mask2_new(1,:));
data2.y_coor_wss_new = yxz_coor_wss_new(1,:)';% - min(YXZ_mask2_new(2,:));
data2.z_coor_wss_new = yxz_coor_wss_new(3,:)';% - min(YXZ_mask2_new(3,:));
%clear yxz_coor_wss_new

%%% WSS coordinates
data2.x_value_wss_new = yxz_value_wss_new(2,:)';% - min(YXZ_mask2_new(1,:));
data2.y_value_wss_new = yxz_value_wss_new(1,:)';% - min(YXZ_mask2_new(2,:));
data2.z_value_wss_new = yxz_value_wss_new(3,:)';% - min(YXZ_mask2_new(3,:));
data2.wss_m_new = sqrt(data2.x_value_wss_new.^2 + data2.y_value_wss_new.^2 + data2.z_value_wss_new.^2);
%clear yxz_value_wss_new

figure('Name','Registered')
plot3(data1.x_coor_vel,data1.y_coor_vel,data1.z_coor_vel,'r.')
hold on
plot3(data2.x_coor_vel_new,data2.y_coor_vel_new,data2.z_coor_vel_new,'b.')
legend('remains the same','transformed')
axis equal; axis ij; axis off; view([-180 -90])  

x_vel_round = round(data2.x_coor_vel_new./data1.vox(1));
y_vel_round = round(data2.y_coor_vel_new./data1.vox(2));
z_vel_round = round(data2.z_coor_vel_new./data1.vox(3));

% This is a little trick to prevent negative coordinates.
% If there are any negative coordinates, we'll get into trouble later with negavtive indices
[I1,J1] = find(x_vel_round <= 0 | y_vel_round <= 0 | z_vel_round <= 0);
[I2,J2] = find(x_vel_round > 0 & y_vel_round > 0 & z_vel_round > 0);
if ~isempty(I1)
    x_vel_round(I1,1) = x_vel_round(I2(1),1);
    y_vel_round(I1,1) = y_vel_round(I2(1),1);
    z_vel_round(I1,1) = z_vel_round(I2(1),1);
end

% figure('Name','Registered2')
% plot3(data1.x_coor_vel./data1.vox(1),data1.y_coor_vel./data1.vox(2),data1.z_coor_vel./data1.vox(3),'r.')
% hold on
% plot3(x_vel_round,y_vel_round,z_vel_round,'b.')
% legend('remains the same','transformed')
% axis equal;% axis off

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

% figure('Name','mask2_new')
% V = vol3d('cdata',mask2_new,'texture','3D','texturemap',mask2_new);
% axis tight; axis equal

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

%L_mask1 = mask1_new;%double(mask1_new ~= 0);
L_mask2 = double(mask2_new ~= 0);

% figure('Name','mask1_filtered')
% V = vol3d('cdata',mask1,'texture','3D','texturemap',mask1);
% axis tight; axis equal
%
% figure('Name','mask2_filtered')
% V = vol3d('cdata',L_mask2,'texture','3D','texturemap',L_mask2);
% axis tight; axis equal
%
% figure('Name','difference')
difference = abs(mask1-L_mask2);
% V = vol3d('cdata',difference,'texture','3D','texturemap',difference);colorbar
% axis tight; axis equal

[I1,J] = find(mask1~=0);
[I2,J] = find(L_mask2~=0);
mean_I = (size(I1,1) + size(I2,1))/2;
[I_diff,J] = find(difference~=0);
diff_voxels = size(I_diff,1);
diff_percentage = (diff_voxels / mean_I) * 100;
disp(['Difference between aorta and atlas = ' num2str(diff_percentage)])

% figure('Name','Registered WSS')
% plot3(data1.x_coor_wss,data1.y_coor_wss,data1.z_coor_wss,'r.')
% hold on
% plot3(data2.x_coor_wss,data2.y_coor_wss,data2.z_coor_wss,'b.')
% legend('remains the same','transformed')
% axis equal; axis off

interpolation_function = TriScatteredInterp([data1.x_coor_wss data1.y_coor_wss data1.z_coor_wss],atlas_SD.wssmean,'nearest');
atlas_mean2 = interpolation_function([data2.x_coor_wss_new data2.y_coor_wss_new data2.z_coor_wss_new]);
atlas_mean2(isnan(atlas_mean2)) = 0;
atlas_mean = atlas_mean2;

interpolation_function = TriScatteredInterp([data1.x_coor_wss data1.y_coor_wss data1.z_coor_wss],atlas_SD.wssSD,'nearest');
atlas_SD2 = interpolation_function([data2.x_coor_wss_new data2.y_coor_wss_new data2.z_coor_wss_new]);
atlas_SD2(isnan(atlas_SD2)) = 0;
atlas_SD = atlas_SD2;

% figure('Name','Mean')
% patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',atlas_mean,'FaceColor','interp','FaceAlpha',1);colorbar;
% axis equal;axis off
% caxis([0 1.5])
% view([90 45])
% zoom(3)
%
% figure('Name','SD')
% patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',atlas_SD,'FaceColor','interp','FaceAlpha',1);colorbar;
% axis equal;axis off
% caxis([0 1.5])
% view([90 45])
% zoom(3)

mean_plus_2SD_atlas = atlas_mean + 1.96.*atlas_SD;
mean_minus_2SD_atlas = atlas_mean - 1.96.*atlas_SD;

% figure('Name','mean_plus_2SD_atlas ')
% patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',mean_plus_2SD_atlas ,'FaceColor','interp','FaceAlpha',1);colorbar;
% axis equal;axis off
% caxis([0 1.5])
% view([90 45])
% zoom(3)
%
% figure('Name','mean_min_2SD_atlas ')
% patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',mean_minus_2SD_atlas ,'FaceColor','interp','FaceAlpha',1);colorbar;
% axis equal;axis off
% caxis([0 1.5])
% view([90 45])
% zoom(3)
% pause

% figure('Name','Aorta shape')
% patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
% axis equal;axis off
% caxis([0 1.5])
% view([90 45])
% zoom(3)
%
% figure('Name','Aorta')
% patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss],'EdgeColor','none', 'FaceVertexCData',data2.wss_m,'FaceColor','interp','FaceAlpha',1);colorbar;
% axis equal;axis off
% caxis([0 1.5])
% view([90 45])
% zoom(3)
% pause

wss_m_BAV = data2.wss_m;%sqrt(data2.x_value_wss.^2 + data2.y_value_wss.^2 + data2.z_value_wss.^2);

heat_mapp = zeros(size(wss_m_BAV,1),1);
for i=1:size(wss_m_BAV,1)
    if wss_m_BAV(i) < mean_minus_2SD_atlas(i)
        heat_mapp(i,1) = 0;
    elseif wss_m_BAV(i) > mean_plus_2SD_atlas(i)
        heat_mapp(i,1) = 1;
    else
        heat_mapp(i,1) = 2;
    end
end

disp('done')

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

f2 = figure('Name','Higher/Lower than 1.96*SD')
patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss_new data2.y_coor_wss_new data2.z_coor_wss_new],'EdgeColor','none', 'FaceVertexCData',heat_mapp,'FaceColor','interp','FaceAlpha',1);
%colormap(jet)
colormap(color1);
%c1 = colorbar;
caxis([0 2]);
%set(c1,'YTick',[0.375 1.125 1.875 2.625]);
%set(c1,'YTickLabel',{'p>0.05, control<dilated' 'p>0.05, control>dilated' 'p<0.05, control<dilated' 'p<0.05, control>dilated' })
axis equal;axis ij; axis off
view([-180 -90])  

figure('Name','before/after registration')
plot3(data2.x_coor_wss,data2.y_coor_wss,data2.z_coor_wss,'r.')
hold on
plot3(data2.x_coor_wss_new,data2.y_coor_wss_new,data2.z_coor_wss_new,'b.')
axis ij; axis equal; view([-180 -90])  

%%% WSS coordinates
yxz_coor_wss = [data2.y_coor_wss_new data2.x_coor_wss_new data2.z_coor_wss_new];
yxz_coor_wss(:,4) = 1;
yxz_coor_wss_new = worldmat*yxz_coor_wss'; clear YXZ_coor_wss

%%% WSS coordinates
x_coor_wss = yxz_coor_wss_new(2,:)';% - min(YXZ_mask2_new(1,:));
y_coor_wss = yxz_coor_wss_new(1,:)';% - min(YXZ_mask2_new(2,:));
z_coor_wss = yxz_coor_wss_new(3,:)';% - min(YXZ_mask2_new(3,:));

% figure('Name','before/after registration')
% plot3(data2.x_coor_wss,data2.y_coor_wss,data2.z_coor_wss,'r.')
% hold on
% plot3(x_coor_wss,y_coor_wss,z_coor_wss,'b.')
% axis equal

gray_colormap = colormap(gray);
color1(1,:) = [0 0 1];
color1(2,:) = [1 0 0];
color1(3,:) = [0.5 0.5 0.5];
color1(4:64,:) = gray_colormap(4:64,:);


count = 0;
angles(1) = 0;
f1 = figure('Name','Higher/Lower than 1.96*SD2');
x = x_coor_wss./data2.vox(1);
y = y_coor_wss./data2.vox(2);
z = z_coor_wss./data2.vox(3);
p=patch('Faces',data2.F_matrix{1},'Vertices',[x y z],'EdgeColor','none', 'FaceVertexCData',heat_mapp,'FaceColor','interp','FaceAlpha',1);
colormap(color1);
caxis([0 64]);
axis equal; axis ij; axis off;
magnitude = flipdim(double(data2.mag_struct.dataAy(:,:,:,5)),3);
magnitude(magnitude == 0) = 3;
magnitude(magnitude == 1) = 3;
magnitude(magnitude == 2) = 3;
hold on
s = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3),magnitude(:,:,1),'EdgeColor','none');
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
    'Callback', {@rotater,gca});

uicontrol('Style','text',...
    'Position',[10 75 120 20],...
    'String','Slice slider')
uicontrol('Style', 'slider',...
    'Min',1,'Max',size(magnitude,3),'Value',1,...
    'Position', [10 50 120 20],...
    'SliderStep',[1/(size(magnitude,3)-1) 10/(size(magnitude,3)-1)],...
    'Callback', {@move_slice,gca});

uicontrol('Style','text',...
    'Position',[440 375 120 20],...
    'String','Show Heat Map')
uicontrol('Style','checkbox',...
     'Value',1, 'Position', [440 375 20 20], ...
     'Callback', {@show_heat_mapp,gca});

    function move_slice(hObj,event,ax)       
        sliceobj = findobj(s);
        delete(sliceobj)
        slice = round(get(hObj,'Value'));        
        s = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]).*(size(magnitude,3)-(slice-1)),magnitude(:,:,slice),'EdgeColor','none');
    end

    function show_heat_mapp(hObj,event,ax)       
        show = round(get(hObj,'Value'));      
        if show == 1
             patchobj = findobj(p);
             set(patchobj,'HandleVisibility','on','Visible','on');            
        elseif show == 0
             patchobj = findobj(p);
             set(patchobj,'HandleVisibility','off','Visible','off');
        end
    end

    function rotater(hObj,event,ax)
        count = count + 1;
        dtheta = get(hObj,'Value');     
        dphi = 0;
        
        if count == 1;
            dtheta2 = dtheta*3;
        else
            dtheta2 = (dtheta - angles(count-1))*3;         
        end            
     
        camorbit(dtheta2,dphi,'data',[0 1 0])
        angles(count) = dtheta;
    end

set(f1,'toolbar','figure');

heat_map.heat_map = heat_mapp;
heat_map.vertices = [x y z];
heat_map.faces = data2.F_matrix{1};
heat_map.color = color1;

% set up results folder
dir_orig = pwd;
dir_new = PATHNAME; cd(dir_new); cd('..')
dir_new = pwd;
mkdir('results_heatmap')
dir_new = strcat(dir_new,'\results_heatmap');

% save results in results folder
save(strcat(dir_new,'\heat_map'),'heat_map')
savefig(f2,strcat(dir_new,'\heat_map'))
print(f2,'-zbuffer', '-djpeg','-r600',strcat(dir_new,'\heat_map.jpg'));
cd(dir_orig)
end