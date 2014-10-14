clc, clear, close all

PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\1_M_20140404_094040_Aera_NMH';
PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\2_F_20140404_082635_Aera_NMH';
PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\3_M_20120427_100153_Espree_NMH';
PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\4_M_20120629_141637_Espree_NMH';
PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\5_M_20120906_064323_Espree_NMH';
PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\6_M_20121213_152738_Skyra_NMH';
PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\7_M_20130928_132456_Avanto_NMH';
PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\8_F_20130620_163030_Skyra_NMH';
PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\9_F_20120705_092028_Skyra_NMH';
PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\10_F_20121221_100044_Skyra_NMH';
PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\11_F_20130620_144334_Skyra_NMH';
PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30\12_F_20140502_075511_Aera_NMH';

% PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\1_M_20120420_132106_Espree_NMH';
% PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\2_M_20120606_140221_Skyra_NMH';
% PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\3_M_20120703_133227_Skyra_NMH';
% PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\4_M_20120705_161151_Skyra_NMH';
% PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\5_M_20120906_135536_Espree_NMH';
% PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\6_M_20121106_140410_Skyra_NMH';
% PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\7_M_20121213_110744_Skyra_NMH';
% PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\8_M_20130403_111332_Skyra_NMH';
% PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\9_F_CAMRI_BOUJOF';
% PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\10_F_20121004_132900_Skyra_NMH';
% PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\11_F_20121025_073209_Skyra_NMH';
% PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40\12_F_JACKSON_0666066048';

% PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\1_M_20120502_134311_Espree_NMH';
% PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\2_M_20120608_094441_Skyra_NMH';
% PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\3_M_20120627_093614_Espree_NMH';
% PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\4_M_20120702_092347_Espree_NMH';
% PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\5_M_20120822_092013_Skyra_NMH';
% PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\6_M_20121008_082605_Skyra_NMH';
% PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\7_M_20121022_101001_Skyra_NMH';
% PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\8_M_20121106_102335_Skyra_NMH';
% PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\9_M_21021112_141522_Skyra_NMH';
% PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\10_M_2012114_151338_Skyra_NMH';
% PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\11_M_20130306_071348_Aera_NMH';
% PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\12_M_20130513_134314_Aera_NMH';
% PATHNAME{13} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\13_M_20130516_101643_Skyra_NMH';
% PATHNAME{14} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\14_F_20120712_112128_Skyra_NMH';
% PATHNAME{15} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\15_F_20121029_075439_Skyra_NMH';
% PATHNAME{16} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\16_F_20121108_080242_Skyra_NMH';
% PATHNAME{17} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\17_F_20121203_081058_Skyra_NMH';
% PATHNAME{18} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\18_F_20140313_115319_Aera_NMH';
% PATHNAME{19} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50\19_F_20140408_083611_Aera_NMH';
% % % % % 
% PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\1_M_20120522_170003_Skyra_NMH';
% PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\2_M_20121130_124948_Skyra_NMH';
% PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\3_M_20130118_131920_Aera_NMH';
% PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\4_M_20130516_130638_Skyra_NMH';
% PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\5_M_20130620_084441_Aera_NMH';
% PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\6_M_20140508_100742_Skyra_NMH';
% PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\7_M_20130905_090323_Skyra_NMH';
% PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\8_M_20140313_102146_Aera_NMH';
% PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\9_M_20140414_085925_Aera_NMH';
% PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\10_M_20140421_090435_Skyra_NMH';
% PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\11_M_20140506_081402_Aera_NMH';
% PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\12_F_20120605_114952_Skyra_NMH';
% PATHNAME{13} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60\13_F_20121109_080007_Skyra_NMH';

% PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\1_M_20130904_083954_Skyra_NMH';
% PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\2_F_20140318_084120_Aera_NMH';
% PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\3_M_20121128_130955_Skyra_NMH';
% PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\4_M_20130620_084441_Aera_NMH';
% PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\5_M_20130118_131920_Aera_NMH';
% PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\6_M_20121130_124948_Skyra_NMH';
% PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\7_M_20140416_081355_Aera_NMH';
% PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\8_M_20140421_090435_Skyra_NMH';
% PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80\9_F_20121109_080007_Skyra_NMH';

% PATHNAME{1} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\1_20120420_132106';
% PATHNAME{2} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\2_20120426_132244';
% PATHNAME{3} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\3_20121206_115454';
% PATHNAME{4} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\4_20120522_170003';
% PATHNAME{5} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\5_20120502_134311';
% PATHNAME{6} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\6_20120627_093614';
% PATHNAME{7} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\7_20120702_092347';
% PATHNAME{8} = 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\8_20120831_101148';
% PATHNAME{9} = 'l:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\9_20121109_080007';
% PATHNAME{10}= 'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls\10_20130621_122315';
% 
% PATHNAME1{1} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\1_PT5_GD';
% PATHNAME1{2} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\2_PT9_LD';
% PATHNAME1{3} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\3_PT12_MK';
% PATHNAME1{4} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\4_PT15_PR';
% PATHNAME1{5} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\5_PT17_SD';
% PATHNAME1{6} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\6_PT20_WB';
% PATHNAME1{7} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\7_PT25_SE';
% PATHNAME1{8} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\8_PT35_KP';
% PATHNAME1{9} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\9_PT39_DL';
% PATHNAME1{10} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\10_PT42_BD';
% PATHNAME1{11} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\11_PT44_RS';
% PATHNAME1{12} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\12_PT60_KG';
% PATHNAME1{13} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RL\13_PT64_JD';
% 
% PATHNAME2{1} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\1_PT18_SE';
% PATHNAME2{2} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\2_PT152_OB';
% PATHNAME2{3} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\3_PT232_KS';
% PATHNAME2{4} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\4_PT240_MS';
% PATHNAME2{5} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\5_PT143_BT';
% PATHNAME2{6} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\6_PT150_FC';
% PATHNAME2{7} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\7_PT11_MP';
% PATHNAME2{8} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\8_PT28_AT';
% PATHNAME2{9} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\9_PT66_JF';
% PATHNAME2{10} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\10_PT10_MA';
% PATHNAME2{11} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\11_PT256_MW';
% PATHNAME2{12} = 'C:\1_Chicago\Data\MIMICS\RL_RN_comparison\RN\12_PT254_JD';

% % RL-RN age
% Age(1) = 45;
% Age(2) = 60;
% Age(3) = 57;
% Age(4) = 52;
% Age(5) = 41;
% Age(6) = 42;
% Age(7) = 43;
% Age(8) = 33;
% Age(9) = 32;
% Age(10) = 40;
% Age(11) = 29;
% Age(12) = 51;
% Age(13) = 58;
% Age(14) = 64;
% Age(15) = 62;
% Age(16) = 49;
% Age(17) = 58;
% Age(18) = 32;
% Age(19) = 48;
% Age(20) = 44;
% Age(21) = 45;
% Age(22) = 46;
% Age(23) = 33;
% Age(24) = 40;
% Age(25) = 44;
% Age(26) = 38;
% Age(27) = 34;
% Age(28) = 53;
% Age(29) = 55;
% Age(30) = 49;

% PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\1_PT76_MG_RL\20120907';
% PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\2_PT42_BD_RL\20120626';
% PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\3_PT33_GB_RL\20120425';
% PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\4_PT55_ML_RN\20120717';
% PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\5_PT24_AM_True\20120111';
% PATHNAME{11} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\6_PT13_MS_True\20111101';
% PATHNAME{13} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\7_PT11_MP_RN\20111223';
% PATHNAME{15} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\8_PT103_DN_unicuspid\20121026';
% PATHNAME{17} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\9_PT20_WB_RL\20120111';
% PATHNAME{19} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\10_PT22_YM_True\20111227';
% PATHNAME{21} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\11_PT2_DM_unicuspid\20120622';
% PATHNAME{23} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\12_PT232_KS_RN\20130606';
% PATHNAME{25} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\13_PT113_PM_RN\20120904';
% PATHNAME{27} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\14_PT19_SD_functionally_unicuspid\20120116';
% PATHNAME{29} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\15_PT166_KW_RL\20130212';
% PATHNAME{31} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\16_PT78_KK_RL\20120829';
% PATHNAME{33} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\17_PT119_HM_RL\20121127';
% PATHNAME{35} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\18_PT48_HK_RL\20120814';
% PATHNAME{37} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\19_PT12_MK_RL\20111222';
% PATHNAME{39} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\20_PT86_FE_LN\20120927';
% PATHNAME{41} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\21_PT8_JD_RL\20120120';
% PATHNAME{43} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\22_PT53_JJ_RL\20120713';
% PATHNAME{45} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\23_PT168_EP_LN_very_aliased\20130218';
% PATHNAME{47} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\24_PT172_MB_RN\20130220';
% 
% PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\1_PT76_MG_RL\20130325';
% PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\2_PT42_BD_RL\20130528';
% PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\3_PT33_GB_RL\20130529';
% PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\4_PT55_ML_RN\20130701';
% PATHNAME{10} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\5_PT24_AM_True\20130204';
% PATHNAME{12} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\6_PT13_MS_True\20121227';
% PATHNAME{14} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\7_PT11_MP_RN\20131021';
% PATHNAME{16} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\8_PT103_DN_unicuspid\20130904';
% PATHNAME{18} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\9_PT20_WB_RL\20130304';
% PATHNAME{20} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\10_PT22_YM_True\20130125';
% PATHNAME{22} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\11_PT2_DM_unicuspid\20121228';
% PATHNAME{24} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\12_PT232_KS_RN\20131210';
% PATHNAME{26} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\13_PT113_PM_RN\20130903';
% PATHNAME{28} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\14_PT19_SD_functionally_unicuspid\20130129';
% PATHNAME{30} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\15_PT166_KW_RL\20130806';
% PATHNAME{32} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\16_PT78_KK_RL\20130625';
% PATHNAME{34} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\17_PT119_HM_RL\20140116';
% PATHNAME{36} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\18_PT48_HK_RL\20130322';
% PATHNAME{38} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\19_PT12_MK_RL\20121224';
% PATHNAME{40} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\20_PT86_FE_LN\20130325';
% PATHNAME{42} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\21_PT8_JD_RL\20130220';
% PATHNAME{44} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\22_PT53_JJ_RL\20130731';
% PATHNAME{46} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\23_PT168_EP_LN_very_aliased\20130913';
% PATHNAME{48} = 'C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\24_PT172_MB_RN\20130820';

% % % Follow-up age
% Age(1) = 29;
% Age(3) = 29;
% Age(5) = 33;
% Age(7) = 34;
% Age(9) = 37;
% Age(11) = 37;
% Age(13) = 40;
% Age(15) = 41;
% Age(17) = 42;
% Age(19) = 42;
% Age(21) = 44;
% Age(23) = 44;
% Age(25) = 45;
% Age(27) = 46;
% Age(29) = 52;
% Age(31) = 52;
% Age(33) = 54;
% Age(35) = 55;
% Age(37) = 57;
% Age(39) = 58;
% Age(41) = 61;
% Age(43) = 64;
% Age(45) = 64;
% Age(47) = 76;
% Age(2) = 29;
% Age(4) = 29;
% Age(6) = 33;
% Age(8) = 34;
% Age(10) = 37;
% Age(12) = 37;
% Age(14) = 40;
% Age(16) = 41;
% Age(18) = 42;
% Age(20) = 42;
% Age(22) = 44;
% Age(24) = 44;
% Age(26) = 45;
% Age(28) = 46;
% Age(30) = 52;
% Age(32) = 52;
% Age(34) = 54;
% Age(36) = 55;
% Age(38) = 57;
% Age(40) = 58;
% Age(42) = 61;
% Age(44) = 64;
% Age(46) = 64;
% Age(48) = 76;
% 
% % % load C:\1_Chicago\Data\MIMICS\17_Followup\2_BAV\probability_mask_baseline.mat
% probability_mask = '';
% [p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,probability_mask,1,0,0,1,0,1)
% 
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
for n = 1:size(PATHNAME,2)
        
    MrstructPath = strcat(PATHNAME{n},'\mrstruct\');
    MimicsSegPath = strcat(PATHNAME{n},'\results_022\');
    
if ~exist([MrstructPath 'Wss_point_cloud_aorta.mat'])
    tic
    mimics_to_Wss(MrstructPath,MimicsSegPath,1,1,1,1,0,0,0);
    toc 
end
    close all
end
% 
% Make idealized aorta geometry

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

% tic
% [probability_mask] = make_geometry_point_cloud(PATHNAME,1,1);
% toc

% % 
% close all
% 

% % %Make atlas
% PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80';
% [atlas] = make_atlas_point_cloud_scalars_affine_registration(PATHNAME,PATHNAME_probability_mask,1,1,1,1);
% %save the atlas to the directory of choice
% directory = uigetdir('C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\');
% disp('...saving atlas...')
% save(strcat(directory,'\atlas'),'atlas')
% disp(['saved to ' directory])
% 
% pause

% PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80';
% [atlas_tps] = make_atlas_time_to_peak_systole(PATHNAME,PATHNAME_probability_mask,1,1,1);
% directory = uigetdir('C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\');
% disp('...saving atlas...')
% save(strcat(directory,'\atlas_tps'),'atlas_tps')
% disp(['saved to ' directory])


AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80';%'L:\data\NU\Aorta-4D_Flow\Results\Pim\Data\MIMICS\BAV_tissue\Controls'
for n = 1:size(PATHNAME,2)
% %         
% % 
% %     
% %     Age(n)
% %     
% %     if Age(n) >= 18 && Age(n) <= 30
% %         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_18_30'
% %     elseif Age(n) >= 31 && Age(n) <= 40
% %         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_31_40'  
% %     elseif Age(n) >= 41 && Age(n) <= 50
% %         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_41_50'          
% %     elseif Age(n) >= 51 && Age(n) <= 60
% %         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_51_60' 
% %     elseif Age(n) >= 61
% %         AtlasPath = 'C:\1_Chicago\Data\MIMICS\3_ControlsSagittalView\AgeGroups\Mixed_55_80' 
% %     end
% %     
% %     tic
% %   %  [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,MrstructPath,plotFlag,calculateRE_Flag,calculateIE_Flag,calculate_traffic_light_volumeFlag,calculate_area_of_higherlowerFlag,peak_systolicFlag)    

      [traffic_light,heat_map]=heat_map_traffic_light_scalars_affine_registration(AtlasPath,PATHNAME{n},1,1,1,1,0,1);
    % [traffic_light]=traffic_light_time_to_peak_systole(AtlasPath,PATHNAME{n},1,1,1);      
     pause
%     toc
     clc,close all 
end