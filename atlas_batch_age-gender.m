clc, clear, close all

directory1 = 'X:\cv_mri\Aorta-4D_Flow\controls' % location of all patient folders

[num,txt,raw] = xlsread('Z:\Aorta_Atlas\target_enrollment_atlases_20171110','Atlas Recruitment','G43:G53')    % reads the list of patient folders from spreadsheet

for n = 1:size(txt,1)
    path = raw(n,1);
    PATHNAME{n} =  [directory1 '\' path{1} '\3dpc'];
end

% Calculate WSS
% for n = 1:size(PATHNAME,2)
% 
%     path = PATHNAME{n}
% 
% %    MrstructPath = strcat(path,'\3dpc\mrstruct\');
% %    MimicsSegPath = strcat(path,'\3dpc\results_022\');%strcat(path,'\3dpc\results_022\');
%     
%     currDir=pwd;
% 
%     cd(path)
%     folders = ls;
%     if exist('mrstruct','dir')==7
%         MrstructPath = strcat(PATHNAME{n},'\mrstruct\3dpc\');
%     else
%         for i=3:size(folders,1)
%             [a,b]=find(folders(i,1:8)=='mrstruct');
%             if sum(a)==8
%                 MrstructPath=strcat(PATHNAME{n},'\3dpc\',folders(i,:),'\');
%                 break
%             end
%         end
%     end
% 
%    MrstructPath = strcat(PATHNAME{n},'\mrstruct\')
% 
% if exist(strcat(path,'\results_022'))
%    MimicsSegPath = strcat(path,'\results_022\')
% elseif exist(strcat(path,'\results_user025'))
%    MimicsSegPath = strcat(path,'\results_user025\')
% end
% 
% if ~exist([MrstructPath 'Wss_point_cloud_aorta.mat'])
% % if  exist(MrstructPath) && exist(MimicsSegPath) && ~exist([MrstructPath 'Wss_point_cloud_aorta.mat'])
% %     tic
%     mimics_to_Wss(MrstructPath,MimicsSegPath,1,1,1,1,0,0,0);
%     toc
% end
%     close all
% end

tic
[probability_mask] = make_geometry_point_cloud(PATHNAME,1,1,2,4);
toc

% Make atlas velocity & WSS
PATHNAME_probability_mask = 'C:\Users\Emily\Documents\MR_processing\heatmap_age-matching\group1\p_mask';
[atlas] = make_atlas_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,1,1,0,1);
% save the atlas to the directory of choice
directory = uigetdir('C:\Users\Emily\Documents\MR_processing\heatmap_age-matching\group1\atlas');
disp('...saving atlas...')
save(strcat(directory,'\atlas'),'atlas')
disp(['saved to ' directory])

