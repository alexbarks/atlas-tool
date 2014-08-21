clc, clear, close all

% % PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\1_M_20140404_094040_Aera_NMH\mrstruct\';
% % PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\2_F_20140404_082635_Aera_NMH\mrstruct\';
% % PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\3_M_20120427_100153_Espree_NMH\mrstruct\';
% % PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\4_M_20120629_141637_Espree_NMH\mrstruct\';
% % PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\5_M_20120906_064323_Espree_NMH\mrstruct\';
% % PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\6_M_20121213_152738_Skyra_NMH\mrstruct\';
% % PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\7_M_20130928_132456_Avanto_NMH\mrstruct\';
% % PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\8_F_20130620_163030_Skyra_NMH\mrstruct\';
% % PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\9_F_20120705_092028_Skyra_NMH\mrstruct\';
% % PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\10_F_20121221_100044_Skyra_NMH\mrstruct\';
% % PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\11_F_20130620_144334_Skyra_NMH\mrstruct\';
% % PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\12_F_20140502_075511_Aera_NMH\mrstruct\';
% 
% % PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\1_M_20120420_132106_Espree_NMH\mrstruct\';
% % PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\2_M_20120606_140221_Skyra_NMH\mrstruct\';
% % PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\3_M_20120703_133227_Skyra_NMH\mrstruct\';
% % PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\4_M_20120705_161151_Skyra_NMH\mrstruct\';
% % PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\5_M_20120906_135536_Espree_NMH\mrstruct\';
% % PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\6_M_20121106_140410_Skyra_NMH\mrstruct\';
% % PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\7_M_20121213_110744_Skyra_NMH\mrstruct\';
% % PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\8_M_20130403_111332_Skyra_NMH\mrstruct\';
% % PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\9_F_CAMRI_BOUJOF\mrstruct\';
% % PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\10_F_20121004_132900_Skyra_NMH\mrstruct\';
% % PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\11_F_20121025_073209_Skyra_NMH\mrstruct\';
% % PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\12_F_JACKSON_0666066048\mrstruct\';
% 
% % PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\1_M_20120502_134311_Espree_NMH\mrstruct\';
% % PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\2_M_20120608_094441_Skyra_NMH\mrstruct\';
% % PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\3_M_20120627_093614_Espree_NMH\mrstruct\';
% % PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\4_M_20120702_092347_Espree_NMH\mrstruct\';
% % PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\5_M_20120822_092013_Skyra_NMH\mrstruct\';
% % PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\6_M_20121008_082605_Skyra_NMH\mrstruct\';
% % PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\7_M_20121022_101001_Skyra_NMH\mrstruct\';
% % PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\8_M_20121106_102335_Skyra_NMH\mrstruct\';
% % PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\9_M_21021112_141522_Skyra_NMH\mrstruct\';
% % PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\10_M_2012114_151338_Skyra_NMH\mrstruct\';
% % PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\11_M_20130306_071348_Aera_NMH\mrstruct\';
% % PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\12_M_20130513_134314_Aera_NMH\mrstruct\';
% % PATHNAME{13} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\14_F_20120712_112128_Skyra_NMH\mrstruct\';
% % PATHNAME{14} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\15_F_20121029_075439_Skyra_NMH\mrstruct\';
% 
% % PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\1_M_20120522_170003_Skyra_NMH\mrstruct\';
% % PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\2_M_20121130_124948_Skyra_NMH\mrstruct\';
% % PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\3_M_20130118_131920_Aera_NMH\mrstruct\';
% % PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\4_M_20130516_130638_Skyra_NMH\mrstruct\';
% % PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\5_M_20130620_084441_Aera_NMH\mrstruct\';
% % PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\6_M_20140508_100742_Skyra_NMH\mrstruct\';
% % PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\7_M_20130905_090323_Skyra_NMH\mrstruct\';
% % PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\8_M_20140313_102146_Aera_NMH\mrstruct\';
% % PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\9_M_20140414_085925_Aera_NMH\mrstruct\';
% % PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\10_M_20140421_090435_Skyra_NMH\mrstruct\';
% % PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\11_M_20140506_081402_Aera_NMH\mrstruct\';
% % PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\12_F_20120605_114952_Skyra_NMH\mrstruct\';
% % PATHNAME{13} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\13_F_20121109_080007_Skyra_NMH\mrstruct\';
% 
% PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\7_M_20140416_081355_Aera_NMH\mrstruct\';
% PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\6_M_20130904_083954_Skyra_NMH\mrstruct\';
% PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\5_M_20140421_090435_Skyra_NMH\mrstruct\';
% PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\4_M_20121128_130955_Skyra_NMH\mrstruct\';
% PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\3_M_20130620_084441_Aera_NMH\mrstruct\';
% PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\2_M_20130118_131920_Aera_NMH\mrstruct\';
% PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\1_M_20121130_124948_Skyra_NMH\mrstruct\';
% PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\11_F_20140318_084120_Aera_NMH\mrstruct\';
% PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\10_F_20121206_115454_Avanto_NMH\mrstruct\';
% PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\9_F_20120817_081853_Skyra_NMH\mrstruct\';

PATHNAME{1} = 'C:\1_Chicago\BAV_tissue\Controls\1_20120420_132106_Espree_NMH_interobserver';
PATHNAME{2} = 'C:\1_Chicago\BAV_tissue\Controls\2_20120426_132244_Espree_NMH_interobserver';
PATHNAME{3} = 'C:\1_Chicago\BAV_tissue\Controls\3_20121206_115454_Avanto_NMH';
PATHNAME{4} = 'C:\1_Chicago\BAV_tissue\Controls\4_20120522_170003_Skyra_NMH_interobserver';
PATHNAME{5} = 'C:\1_Chicago\BAV_tissue\Controls\5_20120502_134311_Espree_NMH_interobserver';
PATHNAME{6} = 'C:\1_Chicago\BAV_tissue\Controls\6_20120627_093614_Espree_NMH_interobserver';
PATHNAME{7} = 'C:\1_Chicago\BAV_tissue\Controls\7_20120702_092347_Espree_NMH_interobserver';
PATHNAME{8} = 'C:\1_Chicago\BAV_tissue\Controls\8_20120831_101148_Espree_NMH';
PATHNAME{9} = 'C:\1_Chicago\BAV_tissue\Controls\9_20121109_080007_Skyra_NMH';
PATHNAME{10} = 'C:\1_Chicago\BAV_tissue\Controls\10_20130621_122315_Aera_NMH';

% Calculate WSS
for n = 1:size(PATHNAME,2)
        
    MrstructPath = strcat(PATHNAME{n},'\mrstruct\');
    MimicsSegPath = strcat(PATHNAME{n},'\results_022\');
    
    tic
    mimics_to_Wss(MrstructPath,MimicsSegPath,1,1,1,1,0,0,0);
    toc
%    close all
end

% Make idealized aorta geometry
tic
[probability_mask] = make_geometry_point_cloud_NEW(PATHNAME,1,0);
toc

close all

% % Make atlas
[atlas] = make_atlas_point_cloud_scalars_affine_registration_NEW(PATHNAME,probability_mask,1,1,0,0);

% save the atlas to the directory of choice
directory = uigetdir('C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\');
disp('...saving atlas...')
save(strcat(directory,'\atlas'),'atlas')

disp(['saved to ' directory])
