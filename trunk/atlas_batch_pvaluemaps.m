
clc, clear, close all
directory1 = 'L:\cv_mri\Aorta-4D_Flow\BAV';%'L:\cv_mri\Aorta-4D_Flow\controls'
directory2 = 'L:\cv_mri\Aorta-4D_Flow\BAV'
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\AP','B2:B182')
for n = 1:size(txt,1)
    path = raw(n,1);
    PATHNAME1{n} =  [directory1 '\' path{1} '\3dpc'];
    PATHNAME_total{n} =  [directory1 '\' path{1} '\3dpc'];
end
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\LAT','B2:B52')
for m = 1:size(txt,1)
    path = raw(m,1);
    PATHNAME2{m} =  [directory2 '\' path{1} '\3dpc'];
    PATHNAME_total{n+m} =  [directory2 '\' path{1} '\3dpc'];
end
PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\BAVcohorts\p_value_maps\AP_vs_LAT';
[p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,PATHNAME_probability_mask,1,0,0,1,0,1)

break

clc, clear, close all
directory1 = 'L:\cv_mri\Aorta-4D_Flow\BAV';%'L:\cv_mri\Aorta-4D_Flow\controls'
directory2 = 'L:\cv_mri\Aorta-4D_Flow\BAV'
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\Sievers1RL','B2:B95')
for n = 1:size(txt,1)
    path = raw(n,1);
    PATHNAME1{n} =  [directory1 '\' path{1} '\3dpc'];
    PATHNAME_total{n} =  [directory1 '\' path{1} '\3dpc'];
end
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\Sievers1RN','B2:B11')
for m = 1:size(txt,1)
    path = raw(m,1);
    PATHNAME2{m} =  [directory2 '\' path{1} '\3dpc'];
    PATHNAME_total{n+m} =  [directory2 '\' path{1} '\3dpc'];
end
PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\BAVcohorts\p_value_maps\1_Sievers1RL_vs_Sievers1RN_no_stenosis';
[p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,PATHNAME_probability_mask,1,0,0,1,0,1)

break

clc, clear, close all
directory1 = 'L:\cv_mri\Aorta-4D_Flow\BAV'
directory2 = 'L:\cv_mri\Aorta-4D_Flow\BAV'
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\controls_v2','C123:C140')
for n = 1:size(txt,1)
    path = raw(n,1);
    PATHNAME1{n} =  [directory1 '\' path{1} '\3dpc'];
    PATHNAME_total{n} =  [directory1 '\' path{1} '\3dpc'];
end
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\Sievers1RL','B136:B155')
for m = 1:size(txt,1)
    path = raw(m,1);
    PATHNAME2{m} =  [directory2 '\' path{1} '\3dpc'];
    PATHNAME_total{n+m} =  [directory2 '\' path{1} '\3dpc'];
end

PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\BAVcohorts\p_value_maps\3_controls_vs_Sievers1RL_moderate_stenosis';
[p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,PATHNAME_probability_mask,1,0,0,1,0,1)

clc, clear, close all
directory1 = 'L:\cv_mri\Aorta-4D_Flow\controls'
directory2 = 'L:\cv_mri\Aorta-4D_Flow\BAV'
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\controls_v2','C85:C140')
for n = 1:size(txt,1)
    path = raw(n,1);
    PATHNAME1{n} =  [directory1 '\' path{1} '\3dpc'];
    PATHNAME_total{n} =  [directory1 '\' path{1} '\3dpc'];
end
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\Sievers1RN','B22:B28')
for m = 1:size(txt,1)
    path = raw(m,1);
    PATHNAME2{m} =  [directory2 '\' path{1} '\3dpc'];
    PATHNAME_total{n+m} =  [directory2 '\' path{1} '\3dpc'];
end

PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\BAVcohorts\p_value_maps\3_controls_vs_Sievers1RN_moderate_stenosis';
[p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,PATHNAME_probability_mask,1,0,0,1,0,1)

break

clc, clear, close all
directory1 = 'L:\cv_mri\Aorta-4D_Flow\controls'
directory2 = 'L:\cv_mri\Aorta-4D_Flow\BAV'
% Siever 0 LAT moderate stenosis
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\controls_v2','C94:C140')
for n = 1:size(txt,1)
    path = raw(n,1);
    PATHNAME1{n} =  [directory1 '\' path{1} '\3dpc'];
    PATHNAME_total{n} =  [directory1 '\' path{1} '\3dpc'];
end
[num,txt,raw] = xlsread('C:\1_Chicago\Data\MIMICS\BAVcohorts\Sievers0LAT','B31:B41')
for m = 1:size(txt,1)
    path = raw(m,1);
    PATHNAME2{m} =  [directory2 '\' path{1} '\3dpc'];
    PATHNAME_total{n+m} =  [directory2 '\' path{1} '\3dpc'];
end

PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\BAVcohorts\p_value_maps\3_controls_vs_Sievers0LAT_moderate_stenosis';
[p_value_map] = make_pvalue_map_point_cloud_scalars_affine_registration(PATHNAME1,PATHNAME2,PATHNAME_probability_mask,1,0,0,1,0,1)