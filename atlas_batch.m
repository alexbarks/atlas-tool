clc, clear, close all

% PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\1_CAMRI-JEHAM\20140307_084523_Aera_NMH\3dpc';
% PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\2_CAMRI-LAJAM\20140313_102146_Aera_NMH\3dpc_nav80';
% PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\3_CAMRI-KIMCO\20140313_115319_Aera_NMH\3dpc_nav80';
% PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\4_CAMRI-THZUM\20140314_073534_Aera_NMH\3dpc_nav80';
% PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\5_CAMRI-MAGIM2\20140314_090452_Aera_NMH\3dpc_nav80';
% PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\6_CAMRI-WIANM\20140317_095816_Aera_NMH\3dpc_nav80';
% PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\7_CAMRI_ARNIF\20140404_082635_Aera_NMH';
% PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\8_CAMRI-FREPAM\20140404_094040_Aera_NMH\3dpc_nav80';
% PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\9_CAMRI_JUVIF\20140408_083611_Aera_NMH\3dpc_nav80';
% PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\10_CAMRI_JAKIM\20140414_085925_Aera_NMH';
% PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\11_CAMRI_VIVIM\20140416_081355_Aera_NMH';
% PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\12_CAMRI_RIJO\20140429_080804_Aera_NMH';
% PATHNAME{13} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\13_CAMRI_MASNF\20140502_075511_Aera_NMH';
% PATHNAME{14} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\14_CAMRI_ANMAM\20140506_081402_Aera_NMH';

PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\1_CAMRI-JEHAM\20140401_085557_Aera_NMH';
PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\2_CAMRI-LAJAM\20140403_094800_Aera_NMH';
PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\3_CAMRI-KIMCO\20140401_152110_Aera_NMH';
PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\4_CAMRI-THZUM\20140328_074601_Area_NMH\3dpc';
PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\5_CAMRI-MAGIM2\20140328_092618_Aera_NMH';
PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\6_CAMRI-WIANM\20140331_084210_Area_NMH\3dpc';
PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\7_CAMRI_ARNIF\20140418_080619_Area_NMH';
PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\8_CAMRI-FREPAM\20140418_090819_Aera_NMH';
PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\9_CAMRI_JUVIF\20140422_082611_Area_NMH';
PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\10_CAMRI_JAKIM\20140501_081300_Aera_NMH';
PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\11_CAMRI_VIVIM\20140430_080819_Aera_NMH';
PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\12_CAMRI_RIJO\20140513_075753_Aera_NMH';
PATHNAME{13} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\13_CAMRI_MASNF\20140516_075402_Aera_NMH';
PATHNAME{14} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\14_CAMRI_ANMAM\20140520_081054_Aera_NMH';
% 
% % % % load C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\probability_mask_baseline.mat
% probability_mask = '';
% [p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,probability_mask,1,0,0,1,0,1)

% PATHNAME{1} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\01_PT221-BJ';
% PATHNAME{2} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\02_PT23-LS';
% PATHNAME{3} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\03_PT68-DD';
% PATHNAME{4} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\04_PT47-MR';
% PATHNAME{5} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\05_PT115-AG';
% PATHNAME{6} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\06_PT151-NJ';
% PATHNAME{7} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\07_PT193-SM';
% PATHNAME{8} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\08_PT202-AP';
% PATHNAME{9} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\09_PT208-KK';
% PATHNAME{10} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\10_PT214-DT';
% PATHNAME{11} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\11_PT6-RG_followup';
% PATHNAME{12} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\12_PT242-JL';
% PATHNAME{13} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\13_PT53-JJ_followup';
% PATHNAME{14} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\14_PT279-MK';
% PATHNAME{15} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\15_PT294_NJ';
% PATHNAME{16} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\16_PT280_AK';
% PATHNAME{17} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\17_PT258-TS';
% PATHNAME{18} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\18_PT333-SH';
% PATHNAME{19} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\19_PT100-FC_follow-up';
% PATHNAME{20} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\20_PT361-SG';
% PATHNAME{21} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\21_PT128-MH_follow-up_3timesteps';
% PATHNAME{22} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\22_PT48-HK_follow-up2';
% PATHNAME{23} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\BAV\23_PT315-SC';

% PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\1_20140508_100742_Skyra_NMH';
% PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\2_20140508_084846_Skyra_NMH';
% PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\3_20140425_073703_Skyra_NMH';
% PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\4_20140430_080819_Aera_NMH';
% PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\5_20140318_084120_Aera_NMH';
% PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\6_20121206_115454_Avanto_NMH';
% PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\7_20121109_080007_Skyra_NMH';
% PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\8_20121130_124948_Skyra_NMH';
% PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\9_20130118_131920_Aera_NMH';
% PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\traffic_light\4_controls\10_20140513_075753_Aera_NMH';

% Calculate WSS
% size(PATHNAME,2)
% for n = 1:size(PATHNAME,2)
%         
%     MrstructPath = strcat(PATHNAME{n},'\mrstruct_interobserver\');
%     MimicsSegPath = strcat(PATHNAME{n},'\results_022\');
%  %   MimicsSegPath = strcat(PATHNAME{n},'\');
%     
% if ~exist([MrstructPath 'Wss_point_cloud_aorta.mat'])
%     tic
%     mimics_to_Wss(MrstructPath,MimicsSegPath,1,1,1,1,0,0,0);
%     toc 
% end
%     close all
% end
% 
% Make idealized aorta geometry
% tic
% [probability_mask] = make_geometry_point_cloud(PATHNAME,1,1);
% toc

% load C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\TR_18_30
% TRs = TR; clear TR;
% for n = 1:size(PATHNAME,2)
%     PATHNAME{n}
%     TR = TRs(n)
%     [tps_matrix]=calculate_time_to_peak_systole(PATHNAME{n},TR,1,1,1);
%     pause(10)
%     close all
% end
% break


% % 
% close all
% 

% % %Make atlas
PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\';
[atlas] = make_atlas_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,1,1,1,1);
%save the atlas to the directory of choice
directory = uigetdir('C:\1_Chicago\Data\MIMICS\18_Test-Retest\');
disp('...saving atlas...')
save(strcat(directory,'\atlas'),'atlas')
disp(['saved to ' directory])

pause

% PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60';
% [atlas_tps] = make_atlas_time_to_peak_systole(PATHNAME,PATHNAME_probability_mask,1,1,1);
% directory = uigetdir('C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\');
% disp('...saving atlas...')
% save(strcat(directory,'\atlas_tps'),'atlas_tps')
% disp(['saved to ' directory])


% for n = 1:size(PATHNAME,2)
%            
%     Age(n)
%     
%     if Age(n) >= 18 && Age(n) <= 30
%         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30'
%     elseif Age(n) >= 31 && Age(n) <= 40
%         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40'  
%     elseif Age(n) >= 41 && Age(n) <= 50
%         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50'          
%     elseif Age(n) >= 51 && Age(n) <= 60
%         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60' 
%     elseif Age(n) >= 61
%         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80' 
%     end
%     
%     tic
% %[traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,PATHNAME,plotFlag,calculateIE_Flag,calculate_area_of_higherlowerFlag,peak_systolicFlag,images_for_surgeryFlag)    
% 
%       [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,PATHNAME{n},0,0,0,1,0);
%     % [traffic_light]=traffic_light_time_to_peak_systole(AtlasPath,PATHNAME{n},1,1,1);      
% %     toc
% 
%      clc,close all 
% end

% probability_mask = '';
% [p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,probability_mask,1,0,0,1,0,1)
