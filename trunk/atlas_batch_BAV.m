clc, clear, close all

directory = 'L:\cv_mri\Aorta-4D_Flow\controls'
%directory = 'L:\cv_mri\Aorta-4D_Flow\BAV'
%directory = 'L:\cv_mri\Aorta-4D_Flow\Aneurysm'
%directory3 = 'L:\cv_mri\Aorta-4D_Flow\BAV';
%filenames = dir(directory);

%PATHNAME = filenames;
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\17_Followup\controls_v2','D59:D67')

for n = 1:size(txt,1)
    path = raw(n,1);
    PATHNAME{n} = [directory '\' path{1} '\3dpc'];
end

size(PATHNAME,2)

%load C:\1_Chicago\Data\MIMICS\BAVcohorts\Sievers_cohorts\Sievers1_RN_LAT\Sievers1_RN_LAT
% for n = 1:size(Sievers1_RN_LAT,1)
% PATHNAME{n} = [directory '\' Sievers1_RN_LAT{n} '\3dpc']
% pause
% end
%
% size(PATHNAME)
% pause
% load C:\1_Chicago\Data\MIMICS\BAVcohorts\1_AP\AP
% for n = 1:size(AP,1)
% PATHNAME2{n} = [directory '\' AP{n} '\3dpc'];
% end


% % Calculate WSS
% for n = 1:size(PATHNAME,1)
%
%     path = PATHNAME(n).name;
%
%     MrstructPath = strcat(path,'\3dpc\mrstruct\');
%     MimicsSegPath = strcat(path,'\3dpc\results_022\');
%
% %if ~exist([MrstructPath 'Wss_point_cloud_aorta.mat'])
% if  exist(MrstructPath) && exist(MimicsSegPath) && ~exist([MrstructPath 'Wss_point_cloud_aorta.mat'])
%     tic
%     mimics_to_Wss(MrstructPath,MimicsSegPath,1,1,1,1,0,0,0);
%     toc
% end
%     close all
% end

tic
[probability_mask] = make_geometry_point_cloud(PATHNAME,1,1,3,7);
toc
pause
% %
% %Make atlas
PATHNAME_probability_mask = '';%'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80';
[atlas] = make_atlas_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,0,1,0,1);
%save the atlas to the directory of choice
directory = uigetdir('C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\');
disp('...saving atlas...')
save(strcat(directory,'\atlas'),'atlas')
disp(['saved to ' directory])
pause

for n = 1:size(PATHNAME,2)
    n
   %    path = [directory '\' PATHNAME(n).name '\3dpc']
    
   path = PATHNAME{n} ;
   
    MrstructPath = strcat(path,'\mrstruct');
    MimicsSegPath = strcat(path,'\results_022');
    
    MrstructPath
    MimicsSegPath
    
  %      if  ~exist([path '\results_heatmap']) && exist(MrstructPath) && exist(MimicsSegPath) && exist([MrstructPath '\Wss_point_cloud_aorta.mat'])
    
     %       cd(path)
    
            filename = dir([path '\scanInfo_*.xls']);
            [NUM,TXT,RAW]=xlsread([path '\' filename(1).name],'scanInfo','B4');
            Age = char(TXT);
            Age = str2num(Age(2:3))
    
            if Age >= 19 && Age <= 35
                AtlasPath = 'C:\1_Chicago\Data\MIMICS\17_Followup\atlases\19_35'
            elseif Age >= 36 && Age <= 45
                AtlasPath = 'C:\1_Chicago\Data\MIMICS\17_Followup\atlases\36_45'
            elseif Age >= 46 && Age <= 55
                AtlasPath = 'C:\1_Chicago\Data\MIMICS\17_Followup\atlases\46_55'
            elseif Age >= 56 && Age <= 79
                AtlasPath = 'C:\1_Chicago\Data\MIMICS\17_Followup\atlases\56_78'                
            end
            %4
      %      cd ..
            [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,PATHNAME{n},1,0,0,1);
    
%     AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60'
%     PATHNAME{n}
%     [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,PATHNAME{n},1,0,0,1);
%     %  [traffic_light]=traffic_light_time_to_peak_systole(AtlasPath,PATHNAME{n},1,1,1);
    pause
    close all
end
%clc,close all
%end


% PATHNAME_probability_mask = '';
% [rushhour,temperature] = make_temperature_map(PATHNAME,PATHNAME_probability_mask,1,1);

% probability_mask = '';
% [p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,probability_mask,1,0,0,1,0,1)
