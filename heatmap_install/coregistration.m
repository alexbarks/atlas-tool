%function coregistration(offset,strMaskSeg1,strMaskSeg2)

%%% function coregistration(strMaskSeg1,strMaskSeg2)
%
% This function coregisters two different masks
%
% 2014, Pim van Ooij, Northwestern University
%
% Input
% 1)strMaskSeg1     : Path and filename of mrstruct mask 1
% 2)strMaskSeg2     : Path and filename of mrstruct mask 2
%
% Usage
% This code is for co-registration of 2 mask mrstructs
%
% Examples:
% coregistration(strMaskSeg1,strMaskSeg2)
% coregistration('C:\mask1.mat',C:\mask2.mat)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

if nargin < 1 || isempty(offset)
    offset = 40;
end

% if ~exist(strMaskSeg1) == 2 || isempty(strMaskSeg1)
%     [mask1_name, mask1_path] = uigetfile('*.txt','Start by loading the Mimics All_greyvalues text file','Multiselect','Off');
%     strMaskSeg1 = [mask1_path,'/',mask1_name];
% end
%
% if ~exist(strMaskSeg2) == 2 || isempty(strMaskSeg2)
%     [all_name, all_path] = uigetfile('*.txt','Start by loading the Mimics All_greyvalues text file','Multiselect','Off');
%     strMaskSeg2 = [mask2_path,'/',mask2_name];
% end

strMaskSeg1 = 'C:\1_Chicago\Data\MIMICS\CEMRA\8\ce_Ao_grayvalues_mask_struct'
strMaskSeg2 = 'C:\1_Chicago\Data\MIMICS\CEMRA\8\pc_Ao_grayvalues_mask_struct'

% Load up mask1 mrstruct
[maskstruct1,path1] = mrstruct_read(strMaskSeg1);
% Load up mask2 mrstruct
[maskstruct2,path2] = mrstruct_read(strMaskSeg2);

mask1 = maskstruct1.dataAy;
mask1_vox = maskstruct1.vox;
mask2 = maskstruct2.dataAy;
mask2_vox = maskstruct2.vox;

%%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
sizes = [size(mask1,1)+offset size(mask1,2)+offset size(mask1,3)+offset];% size(mask1,4) size(mask1,5)];
mask1b = zeros(sizes);
mask1b((offset+1):size(mask1b,1),(offset+1):size(mask1b,2),(offset+1):size(mask1b,3),:,:) = mask1;
mask1 = mask1b;clear mask1b

L1 = (mask1 ~= 0);
% create all coordinates
[x,y,z] = meshgrid((1:size(mask1,2)).* mask1_vox(2), ...
    (1:size(mask1,1)).* mask1_vox(1),(1:size(mask1,3)).* mask1_vox(3));
mask1_x_coor = x(L1);mask1_y_coor = y(L1);mask1_z_coor = z(L1);
clear x, clear y, clear z

L2 = (mask2 ~= 0);
% create all coordinates
[x,y,z] = meshgrid((1:size(mask2,2)).* mask2_vox(2), ...
    (1:size(mask2,1)).* mask2_vox(1),(1:size(mask2,3)).* mask2_vox(3));
mask2_x_coor = x(L2);mask2_y_coor = y(L2);mask2_z_coor = z(L2);
clear x, clear y, clear z

figure('Name','To be registered')
plot3(mask1_x_coor,mask1_y_coor,mask1_z_coor,'r.')
hold on
plot3(mask2_x_coor,mask2_y_coor,mask2_z_coor,'b.')
legend('to remain the same','to be transformed')
axis equal; axis ij; axis off; view([-180 -90])

%%% Registration
PSF = fspecial('gaussian',1,1);
mask1_to_register = mask1;
mask1_to_register = imfilter(mask1_to_register,PSF,'conv');

mask2 = imfilter(mask2,PSF,'conv');

disp(' ')
disp('...Busy registering... Rigid registration (dof = 6)')
disp(' ')
disp('...This can take up to 5 minutes...')

tic
% directory with flirt.exe and cygwin1.dll (use cygwin convention)
fsldir = '/cygdrive/c/1_Chicago/WSSprojectWithAmsterdam/flirt/';

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
system('c:\cygwin\bin\bash runflirt.sh');

% load transformation mask
load Rotation_Translation -ascii

flirtmat = 'Rotation_Translation.mat';
src = infile;
trg = reffile;
[worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
rotmat = worldmat(1:3,1:3);

%%% Velocity coordinates
yxz_coor = [mask2_y_coor mask2_x_coor mask2_z_coor];
yxz_coor(:,4) = 1;
yxz_coor_new = inv(worldmat)*yxz_coor'; clear yxz_coor
mask2_x_coor_new = yxz_coor_new(2,:)';
mask2_y_coor_new = yxz_coor_new(1,:)';
mask2_z_coor_new = yxz_coor_new(3,:)';clear yxz_coor_new

%    if plotFlag == 1
figure('Name','Registered')
plot3(mask1_x_coor,mask1_y_coor,mask1_z_coor,'r.')
hold on
plot3(mask2_x_coor_new,mask2_y_coor_new,mask2_z_coor_new,'b.')
legend('Mask 1','Mask 2')
axis equal; axis off;view([-180 -90]); axis ij
pause(10)
% end