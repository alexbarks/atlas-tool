clc, clear, close all

directory1 = 'L:\cv_mri\Aorta-4D_Flow\controls'
directory2 = 'L:\cv_mri\Aorta-4D_Flow\BAV'
% %directory = 'L:\cv_mri\Aorta-4D_Flow\Aneurysm'
directory = 'L:\cv_mri\Aorta-4D_Flow\BAV';
% %filenames = dir(directory);
%
%PATHNAME = filenames;
%[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\Sievers0AP','B6:B38')
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\controls_v2','C85:C94')
%
for n = 1:size(txt,1)
    path = raw(n,1);
    PATHNAME{n} =  [directory1 '\' path{1} '\3dpc'];
end
pause

% [num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\controls_v2','C85:C140')
% %
% for n = 1:size(txt,1)
%     path = raw(n,1);
%     PATHNAME1{n} =  [directory1 '\' path{1} '\3dpc'];
%     PATHNAME_total{n} =  [directory1 '\' path{1} '\3dpc'];
% end
% 
% [num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\Sievers0AP','B6:B38')
% %
% for m = 1:size(txt,1)
%     path = raw(m,1);
%     PATHNAME2{m} =  [directory2 '\' path{1} '\3dpc'];
%     PATHNAME_total{n+m} =  [directory2 '\' path{1} '\3dpc'];
% end
% 
% size(PATHNAME1,2)
% size(PATHNAME2,2)
% size(PATHNAME_total,2)
%pause
%
% %load C:\1_Chicago\Data\MIMICS\BAVcohorts\Sievers_cohorts\Sievers1_RN_LAT\Sievers1_RN_LAT
% % for n = 1:size(Sievers1_RN_LAT,1)
% % PATHNAME1{n} = [directory '\' Sievers1_RN_LAT{n} '\3dpc']
% % pause
% % end
% %
% % size(PATHNAME)
% % pause
% % load C:\1_Chicago\Data\MIMICS\BAVcohorts\1_AP\AP
% % for n = 1:size(AP,1)
% % PATHNAME2{n} = [directory '\' AP{n} '\3dpc'];
% % end

% % Calculate Diameter
% %matlabpool open
% for n = 1:size(PATHNAME,2)
%     
%     path = PATHNAME{n}%(n).name;
%     
%     %    MrstructPath = strcat(path,'\3dpc\mrstruct\');
%     %    MimicsSegPath = strcat(path,'\3dpc\');%strcat(path,'\3dpc\results_022\');
%     
%     MrstructPath = strcat(PATHNAME{n},'\mrstruct\')
%     
%     if exist(strcat(path,'\results_022'))
%         MimicsSegPath = strcat(path,'\results_022\')
%     elseif exist(strcat(path,'\results_user022'))
%         MimicsSegPath = strcat(path,'\results_user022\')        
%     elseif exist(strcat(path,'\results_user025'))
%         MimicsSegPath = strcat(path,'\results_user025\')
%     end
%     
%     if ~exist([MrstructPath 'Diameter_point_cloud_aorta.mat'])
%         % if  exist(MrstructPath) && exist(MimicsSegPath) && ~exist([MrstructPath 'Wss_point_cloud_aorta.mat'])
%         tic
%         mimics_to_diameter(MrstructPath,MimicsSegPath,1,3,1,1);
%         toc
%         close all
%     end
% end
% 
% % Calculate WSS
% for n = 1:size(PATHNAME,2)
% 
%     path = PATHNAME{n}%(n).name;
% 
% %    MrstructPath = strcat(path,'\3dpc\mrstruct\');
% %    MimicsSegPath = strcat(path,'\3dpc\');%strcat(path,'\3dpc\results_022\');
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

% tic
% [probability_mask] = make_geometry_point_cloud(PATHNAME,1,1,2,4);
% toc
% % % %

%matlabpool close

% % %Make atlas diameter
% PATHNAME_probability_mask = '';%'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80';
% [atlas_diameter] = make_atlas_diameter_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,1,1);
% %save the atlas to the directory of choice
% disp('paused')
% pause
% directory = uigetdir('C:\1_Chicago\Data\MIMICS\BAVcohorts\cohort_averaged_maps\');
% disp('...saving atlas...')
% save(strcat(directory,'\atlas_diameter'),'atlas_diameter')
% disp(['saved to ' directory])

% %Make atlas velocity & WSS
PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\age_matching\controls\19_30';%'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80';
[atlas] = make_atlas_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,1,1,0,1);
%save the atlas to the directory of choice
directory = uigetdir('C:\1_Chicago\Data\MIMICS\BAVcohorts\cohort_averaged_maps\');
disp('...saving atlas...')
save(strcat(directory,'\atlas'),'atlas')
disp(['saved to ' directory])
pause

% for n = 1:size(PATHNAME,2)
%     n
%    %    path = [directory '\' PATHNAME(n).name '\3dpc']
%
%    path = PATHNAME{n} ;
%
%     MrstructPath = strcat(path,'\mrstruct');
%     MimicsSegPath = strcat(path,'\results_022');
%
%     MrstructPath
%     MimicsSegPath
%
%   %      if  ~exist([path '\results_heatmap']) && exist(MrstructPath) && exist(MimicsSegPath) && exist([MrstructPath '\Wss_point_cloud_aorta.mat'])
%
%      %       cd(path)
%
%             filename = dir([path '\scanInfo_*.xls']);
%             [NUM,TXT,RAW]=xlsread([path '\' filename(1).name],'scanInfo','B4');
%             Age = char(TXT);
%             Age = str2num(Age(2:3))
%
%             if Age >= 19 && Age <= 35
%                 AtlasPath = 'C:\1_Chicago\Data\MIMICS\17_Followup\atlases\19_35'
%             elseif Age >= 36 && Age <= 45
%                 AtlasPath = 'C:\1_Chicago\Data\MIMICS\17_Followup\atlases\36_45'
%             elseif Age >= 46 && Age <= 55
%                 AtlasPath = 'C:\1_Chicago\Data\MIMICS\17_Followup\atlases\46_55'
%             elseif Age >= 56 && Age <= 79
%                 AtlasPath = 'C:\1_Chicago\Data\MIMICS\17_Followup\atlases\56_78'
%             end
%             %4
%       %      cd ..
%             [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,PATHNAME{n},1,0,0,1);
%
% %     AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60'
% %     PATHNAME{n}
% %     [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,PATHNAME{n},1,0,0,1);
% %     %  [traffic_light]=traffic_light_time_to_peak_systole(AtlasPath,PATHNAME{n},1,1,1);
%     pause
%     close all
% end
% %clc,close all
% %end


% PATHNAME_probability_mask = '';
% [rushhour,temperature] = make_temperature_map(PATHNAME,PATHNAME_probability_mask,1,1);

% PATHNAME_probability_mask = '';
% [p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,PATHNAME_probability_mask,1,0,0,1,0,1)
