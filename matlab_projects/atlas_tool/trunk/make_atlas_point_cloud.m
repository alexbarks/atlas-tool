clc, clear, close all

%function make_atlas

calculate_interplation_errorFlag = 0;

PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\';
%probability_mask = load(strcat(PATHNAME_probability_mask,'probability_mask.mat'));
load('C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\probability_mask_TEST.mat')

%%%% WSS to load
PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\1_CAMRI_FERGCRF\20120426_132244_Espree_NMH\';
PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\2_CAMRI_JOHNCHRM\20120606_140221_Skyra_NMH\';
PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\3_CAMRI_GALTJAF\20120627_153138_Espree_NMH\';
PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\4_CAMRI_HALLJEM\20120702_092347_Espree_NMH\';
PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\5_CAMRI_ROZDJEM\20120822_092013_Skyra_NMH\';
PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\6_CAMRI_MARLE\20121025_073209_Skyra_NMH\';
PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\7_CAMRI_WILWAD\20121029_075439_Skyra_NMH\';
PATHNAME{8}= 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\8_CAMRI_LUSTSCM\20121114_151338_Skyra_NMH\';
PATHNAME{9}= 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\9_CAMRI_GIRAMAM\20121128_130955_Skyra_NMH\';
PATHNAME{10}= 'C:\1_Chicago\Data\MIMICS\17_Followup\1_Controls\10_CAMRI_JAHAL\20121213_152738_Skyra_NMH\';
FILENAME = 'data_done';

mask1 = probability_mask.matrix;
mask1_vox =probability_mask.vox;clear probability_mask

%%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
sizes = [size(mask1,1)+40 size(mask1,2)+40 size(mask1,3)+40];% size(mask1,4) size(mask1,5)];
mask1b = zeros(sizes);
mask1b(41:size(mask1b,1),41:size(mask1b,2),41:size(mask1b,3),:,:) = mask1;
mask1 = mask1b;clear mask1b

L1 = (mask1 ~= 0);
% create all coordinates
[x,y,z] = meshgrid((1:size(mask1,2)).* mask1_vox(2), ...
    (1:size(mask1,1)).* mask1_vox(1),(1:size(mask1,3)).* mask1_vox(3));
geo.x_coor_vel = x(L1);geo.y_coor_vel = y(L1);geo.z_coor_vel = z(L1);
clear x, clear y, clear z
contours = zeros(size(L1));
contours(L1==0) = -1;
contours(L1==1) = 1;
[geo.F,V] = isosurface(contours,0); % make a surface from the detected contours
geo.V = V .* (ones(size(V,1),1) * mask1_vox(1:3));

WSSx = zeros(size(geo.V,1),size(PATHNAME,2));
WSSy = zeros(size(geo.V,1),size(PATHNAME,2));
WSSz = zeros(size(geo.V,1),size(PATHNAME,2));
for n = 1:size(PATHNAME,2)
    
    %     figure('Name','mask1')
    %     V = vol3d('cdata',L,'texture','3D','texturemap',L);
    %     xlabel('x');ylabel('y');zlabel('z')
    %     axis tight; axis equal;axis ij
    
    load(strcat(PATHNAME{n},FILENAME))
    data2 = data; clear data;
    
    % create coordinates
    if ~isfield(data2,'PC_unaliased')
        data2.PC_unaliased = data2.PC_zeros;
    end
    mask2 = (data2.PC_unaliased(:,:,:,1,1) ~= 0);
    
    %%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
    sizes = [size(mask2,1)+40 size(mask2,2)+40 size(mask2,3)+40];% size(mask1,4) size(mask1,5)];
    mask2b = zeros(sizes);
    mask2b(41:size(mask2b,1),41:size(mask2b,2),41:size(mask2b,3),:,:) = mask2;
    mask2 = mask2b;clear mask2b
        
    L2 = (mask2 ~= 0);
    % create all coordinates
    [x,y,z] = meshgrid((1:size(mask2,2)).* data2.vox(2), ...
        (1:size(mask2,1)).*data2.vox(1),(1:size(mask2,3)).* data2.vox(3));
    data2.x_coor_vel = x(L2);data2.y_coor_vel = y(L2);data2.z_coor_vel = z(L2);
    clear x, clear y, clear z
    contours = zeros(size(L2));
    contours(L2==0) = -1;
    contours(L2==1) = 1;
    [data2.F,V] = isosurface(contours,0); % make a surface from the detected contours
    data2.V = V .* (ones(size(V,1),1) * data2.vox(1:3));    
    
%     [x,y,z] = meshgrid((1:size(mask2,2)).* data2.vox(2), ...
%         (1:size(mask2,1)).* data2.vox(1),(1:size(mask2,3)).* data2.vox(3));
%     data2.x_coor_vel = x(L2);data2.y_coor_vel = y(L2);data2.z_coor_vel = z(L2);
%     clear x, clear y, clear z
% %   mask2 = double(L2);
    
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
    geo.x_coor_wss = geo.V(:,1);
    geo.y_coor_wss = geo.V(:,2);
    geo.z_coor_wss = geo.V(:,3);
    
    data2.x_coor_wss = data2.V_matrix{1}(:,1) + (40*data2.vox(1));
    data2.y_coor_wss = data2.V_matrix{1}(:,2) + (40*data2.vox(2));
    data2.z_coor_wss = data2.V_matrix{1}(:,3) + (40*data2.vox(3));
    
    % WSS averaged for systolic timesteps
    data2.x_value_wss_t1 = data2.WSS_matrix{time2-2}(:,1);data2.y_value_wss_t1 = data2.WSS_matrix{time2-2}(:,2);data2.z_value_wss_t1 = data2.WSS_matrix{3}(:,3);
    data2.x_value_wss_t2 = data2.WSS_matrix{time2-1}(:,1);data2.y_value_wss_t2 = data2.WSS_matrix{time2-1}(:,2);data2.z_value_wss_t2 = data2.WSS_matrix{4}(:,3);
    data2.x_value_wss_t3 = data2.WSS_matrix{time2}(:,1);data2.y_value_wss_t3 = data2.WSS_matrix{time2}(:,2);data2.z_value_wss_t3 = data2.WSS_matrix{5}(:,3);
    data2.x_value_wss_t4 = data2.WSS_matrix{time2+1}(:,1);data2.y_value_wss_t4 = data2.WSS_matrix{time2+1}(:,2);data2.z_value_wss_t4 = data2.WSS_matrix{6}(:,3);
    data2.x_value_wss_t5 = data2.WSS_matrix{time2+2}(:,1);data2.y_value_wss_t5 = data2.WSS_matrix{time2+2}(:,2);data2.z_value_wss_t5 = data2.WSS_matrix{7}(:,3);
    
    data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4 + data2.x_value_wss_t5)./5;
    data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4 + data2.y_value_wss_t5)./5;
    data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4 + data2.z_value_wss_t5)./5;
    
    data2.wss_m = sqrt(data2.x_value_wss.^2 + data2.y_value_wss.^2 + data2.z_value_wss.^2);
    
    figure('Name','data2')
    a = [8 20];
    c = [ ];
    skipfactor = 4;
    %patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss], ...
    %    'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
    %hold on
    [F2,V2,C2]=quiver3Dpatch(data2.x_coor_wss(1:skipfactor:size(data2.x_coor_wss,1)),data2.y_coor_wss(1:skipfactor:size(data2.x_coor_wss,1)),data2.z_coor_wss(1:skipfactor:size(data2.x_coor_wss,1)),data2.x_value_wss(1:skipfactor:size(data2.x_coor_wss,1)) ...
        ,data2.y_value_wss(1:skipfactor:size(data2.x_coor_wss,1)),data2.z_value_wss(1:skipfactor:size(data2.x_coor_wss,1)),c,a);
    patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
    c2=colorbar;caxis([0 max(data2.wss_m(:))])
    caxis([0 1.5])
    %T = text(min(data.V_matrix{t2}(:,1)),min(data.V_matrix{t2}(:,2)),max(data.V_matrix{t2}(:,3))+0.005,['Timestep = ' num2str(t2)]);
    %set(T,'FontSize',24)
    axis equal;axis off; axis ij
    view([0 90])
  
   if calculate_interplation_errorFlag == 1;
        wss_magn_new = sqrt(data2.x_value_wss.^2+data2.y_value_wss.^2+data2.z_value_wss.^2);       
        if ~exist(strcat(PATHNAME{n},'mask1.mat'),'file')
            for i = 1:6
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo)
                region = getPosition(polyAAo);
                %mask = inpolygon(x, y, AAo(:,1), AAo(:,2));
                
                disp('saving, pausing')
                save(strcat([PATHNAME{n} '\mask' num2str(i)]),'region');
                pause
            end
        end
        load(strcat(PATHNAME{n},'\mask1'))
        data2_mask_AAo_inner = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_inner = mean(wss_magn_new(data2_mask_AAo_inner));
        load(strcat(PATHNAME{n},'\mask2'))
        data2_mask_AAo_outer = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_outer = mean(wss_magn_new(data2_mask_AAo_outer));
        load(strcat(PATHNAME{n},'\mask3'))
        data2_mask_arch_inner = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_inner = mean(wss_magn_new(data2_mask_arch_inner));
        load(strcat(PATHNAME{n},'\mask4'))
        data2_mask_arch_outer = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_outer = mean(wss_magn_new(data2_mask_arch_outer));
        load(strcat(PATHNAME{n},'\mask5'))
        data2_mask_DAo_inner = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_inner = mean(wss_magn_new(data2_mask_DAo_inner));
        load(strcat(PATHNAME{n},'\mask6'))
        data2_mask_DAo_outer = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_outer = mean(wss_magn_new(data2_mask_DAo_outer));
        
        mean_wss_before_interpolation(n,1) = mean_wss_asc_inner;
        mean_wss_before_interpolation(n,2) = mean_wss_asc_outer;
        mean_wss_before_interpolation(n,3) = mean_wss_arch_inner;
        mean_wss_before_interpolation(n,4) = mean_wss_arch_outer;
        mean_wss_before_interpolation(n,5) = mean_wss_DAo_inner;
        mean_wss_before_interpolation(n,6) = mean_wss_DAo_outer;
        
        figure(101)
        scatter3(data2.x_coor_wss(data2_mask_AAo_inner),data2.y_coor_wss(data2_mask_AAo_inner),data2.z_coor_wss(data2_mask_AAo_inner),20,wss_magn_new(data2_mask_AAo_inner),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(102)
        scatter3(data2.x_coor_wss(data2_mask_AAo_outer),data2.y_coor_wss(data2_mask_AAo_outer),data2.z_coor_wss(data2_mask_AAo_outer),20,wss_magn_new(data2_mask_AAo_outer),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(103)
        scatter3(data2.x_coor_wss(data2_mask_arch_inner),data2.y_coor_wss(data2_mask_arch_inner),data2.z_coor_wss(data2_mask_arch_inner),20,wss_magn_new(data2_mask_arch_inner),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(104)
        scatter3(data2.x_coor_wss(data2_mask_arch_outer),data2.y_coor_wss(data2_mask_arch_outer),data2.z_coor_wss(data2_mask_arch_outer),20,wss_magn_new(data2_mask_arch_outer),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(105)
        scatter3(data2.x_coor_wss(data2_mask_DAo_inner),data2.y_coor_wss(data2_mask_DAo_inner),data2.z_coor_wss(data2_mask_DAo_inner),20,wss_magn_new(data2_mask_DAo_inner),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(106)
        scatter3(data2.x_coor_wss(data2_mask_DAo_outer),data2.y_coor_wss(data2_mask_DAo_outer),data2.z_coor_wss(data2_mask_DAo_outer),20,wss_magn_new(data2_mask_DAo_outer),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
    end    
    
    %     figure('Name','mask2')
    %     V = vol3d('cdata',L,'texture','3D','texturemap',L);
    %     xlabel('x');ylabel('y');zlabel('z')
    %     axis tight; axis equal;axis ij
    
    figure('Name',strcat('To be registered',num2str(n)))
    plot3(geo.x_coor_vel,geo.y_coor_vel,geo.z_coor_vel,'r.')
    hold on
    plot3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,'b.')
    plot3(data2.x_coor_wss,data2.y_coor_wss,data2.z_coor_wss,'g.')
    legend('to remain the same','to be transformed')
    axis equal; axis off;view([0 90]); axis ij
    pause
    
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
    cnii=make_nii(mask1_to_register,[mask1_vox(1) mask1_vox(2) mask1_vox(3)]);
    save_nii(cnii,'mask1.nii');
    mnii=make_nii(mask2,[data2.vox(1) data2.vox(2) data2.vox(3)]);
    save_nii(mnii,'mask2.nii');
    
    % run flirt (needs cygwin installed)
    infile = 'mask2.nii';
    reffile = 'mask1.nii';
    outfile = 'rmask2.nii';
    omat = 'Rotation_Translation.mat';
    
    disp('affine registration (dof = 12)')
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
    
    %%% Velocity
    yxz_coor_vel = [data2.y_coor_vel data2.x_coor_vel data2.z_coor_vel];clear x_coor_vel2, clear y_coor_vel2, clear z_coor_vel2
    yxz_coor_vel(:,4) = 1;
    yxz_coor_vel_new = inv(worldmat)*yxz_coor_vel'; clear yxz_coor_vel
    
    %%% Velocity
    x_coor_vel = yxz_coor_vel_new(2,:)';
    y_coor_vel = yxz_coor_vel_new(1,:)';
    z_coor_vel = yxz_coor_vel_new(3,:)';
    clear yxz_coor_vel_new
    
    %%% WSS coordinates
    yxz_coor_wss = [data2.y_coor_wss data2.x_coor_wss data2.z_coor_wss];clear x_coor_wss2, clear y_coor_wss2, clear z_coor_wss2
    yxz_coor_wss(:,4) = 1;
    yxz_coor_wss_new = inv(worldmat)*yxz_coor_wss'; clear yxz_coor_vel
    
    %%% WSS coordinates
    x_coor_wss = yxz_coor_wss_new(2,:)';
    y_coor_wss = yxz_coor_wss_new(1,:)';
    z_coor_wss = yxz_coor_wss_new(3,:)';
    clear yxz_coor_wss_new
    
    %%% WSS values
    yxz_value_wss = [data2.y_value_wss data2.x_value_wss data2.z_value_wss];
    yxz_value_wss_new = inv(rotmat)*yxz_value_wss'; clear YXZ_value_wss
    
    %%% WSS values
    x_value_wss = yxz_value_wss_new(2,:)';
    y_value_wss = yxz_value_wss_new(1,:)';
    z_value_wss = yxz_value_wss_new(3,:)';
    
    figure('Name','Registered')
    plot3(geo.x_coor_wss,geo.y_coor_wss,geo.z_coor_wss,'r.')
    hold on
    plot3(x_coor_wss,y_coor_wss,z_coor_wss,'b.')
    legend('Prob mask','Aorta')
    %legend('probability mask/atlas','individual aorta')
    axis equal; axis off;view([0 90]); axis ij
    pause
    
    interpolation_function = TriScatteredInterp([x_coor_wss y_coor_wss z_coor_wss],x_value_wss,'nearest');
    x_value_wss = interpolation_function([geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss]);
    x_value_wss(isnan(x_value_wss)) = 0;
    
    interpolation_function = TriScatteredInterp([x_coor_wss y_coor_wss z_coor_wss],y_value_wss,'nearest');
    y_value_wss = interpolation_function([geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss]);
    y_value_wss(isnan(y_value_wss)) = 0;
    
    interpolation_function = TriScatteredInterp([x_coor_wss y_coor_wss z_coor_wss],z_value_wss,'nearest');
    z_value_wss = interpolation_function([geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss]);
    z_value_wss(isnan(z_value_wss)) = 0;
    
    figure('Name','interpolated to atlas')
    skipfactor = 4;
    a = [8 20];
    c = [ ];
    [F2,V2,C2]=quiver3Dpatch(geo.x_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),geo.y_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),geo.z_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:), ...
        x_value_wss(1:skipfactor:size(geo.x_coor_wss,1),:),y_value_wss(1:skipfactor:size(geo.x_coor_wss,1),:),z_value_wss(1:skipfactor:size(geo.x_coor_wss,1),:),c,a);
    patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
    c2=colorbar;caxis([0 1.5])
    axis equal;axis off; axis ij
    view([0 90])
    pause    
    
    if calculate_interplation_errorFlag == 1;
        wss_magn_new = sqrt(x_value_wss.^2+y_value_wss.^2+z_value_wss.^2);    
        if ~exist(strcat(PATHNAME_probability_mask,'mask1.mat'),'file')
            for i = 1:6
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo)
                region = getPosition(polyAAo);
                %mask = inpolygon(x, y, AAo(:,1), AAo(:,2));
                
                disp('saving, pausing')
                save(strcat([PATHNAME_probability_mask '\mask' num2str(i)]),'region');
                pause
            end
        end
        load(strcat(PATHNAME_probability_mask,'\mask1'))
        geo_mask_AAo_inner = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_inner = mean(wss_magn_new(geo_mask_AAo_inner));
        load(strcat(PATHNAME_probability_mask,'\mask2'))
        geo_mask_AAo_outer = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_outer = mean(wss_magn_new(geo_mask_AAo_outer));
        load(strcat(PATHNAME_probability_mask,'\mask3'))
        geo_mask_arch_inner = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_inner = mean(wss_magn_new(geo_mask_arch_inner));
        load(strcat(PATHNAME_probability_mask,'\mask4'))
        geo_mask_arch_outer = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_outer = mean(wss_magn_new(geo_mask_arch_outer));
        load(strcat(PATHNAME_probability_mask,'\mask5'))
        geo_mask_DAo_inner = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_inner = mean(wss_magn_new(geo_mask_DAo_inner));
        load(strcat(PATHNAME_probability_mask,'\mask6'))
        geo_mask_DAo_outer = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_outer = mean(wss_magn_new(geo_mask_DAo_outer));
        
        figure(101)
        scatter3(geo.x_coor_wss(geo_mask_AAo_inner),geo.y_coor_wss(geo_mask_AAo_inner),geo.z_coor_wss(geo_mask_AAo_inner),20,wss_magn_new(geo_mask_AAo_inner),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(102)
        scatter3(geo.x_coor_wss(geo_mask_AAo_outer),geo.y_coor_wss(geo_mask_AAo_outer),geo.z_coor_wss(geo_mask_AAo_outer),20,wss_magn_new(geo_mask_AAo_outer),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(103)
        scatter3(geo.x_coor_wss(geo_mask_arch_inner),geo.y_coor_wss(geo_mask_arch_inner),geo.z_coor_wss(geo_mask_arch_inner),20,wss_magn_new(geo_mask_arch_inner),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(104)
        scatter3(geo.x_coor_wss(geo_mask_arch_outer),geo.y_coor_wss(geo_mask_arch_outer),geo.z_coor_wss(geo_mask_arch_outer),20,wss_magn_new(geo_mask_arch_outer),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(105)
        scatter3(geo.x_coor_wss(geo_mask_DAo_inner),geo.y_coor_wss(geo_mask_DAo_inner),geo.z_coor_wss(geo_mask_DAo_inner),20,wss_magn_new(geo_mask_DAo_inner),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        figure(106)
        scatter3(geo.x_coor_wss(geo_mask_DAo_outer),geo.y_coor_wss(geo_mask_DAo_outer),geo.z_coor_wss(geo_mask_DAo_outer),20,wss_magn_new(geo_mask_DAo_outer),'filled');axis equal;caxis([0 1.5])
        xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        
        mean_wss_after_interpolation(n,1) = mean_wss_asc_inner;
        mean_wss_after_interpolation(n,2) = mean_wss_asc_outer;
        mean_wss_after_interpolation(n,3) = mean_wss_arch_inner;
        mean_wss_after_interpolation(n,4) = mean_wss_arch_outer;
        mean_wss_after_interpolation(n,5) = mean_wss_DAo_inner;
        mean_wss_after_interpolation(n,6) = mean_wss_DAo_outer;
        
        mean_wss_before_interpolation
        mean_wss_after_interpolation
    end
    
    pause(1)
    close all
    
    WSSx(:,n) = x_value_wss;
    WSSy(:,n) = y_value_wss;
    WSSz(:,n) = z_value_wss;
end

atlas.xcoor_wss = geo.x_coor_wss;
atlas.ycoor_wss = geo.y_coor_wss;
atlas.zcoor_wss = geo.z_coor_wss;
atlas.xcoor_vel = geo.x_coor_vel;
atlas.ycoor_vel = geo.y_coor_vel;
atlas.zcoor_vel = geo.z_coor_vel;
atlas.vox = mask1_vox;
atlas.faces = geo.F;
atlas.vertices = geo.V;
atlas.mask = mask1;
atlas.meanx_wss = sum(WSSx,2)./size(PATHNAME,2);
atlas.meany_wss = sum(WSSy,2)./size(PATHNAME,2);
atlas.meanz_wss = sum(WSSz,2)./size(PATHNAME,2);
atlas.mean_wss = sqrt(atlas.meanx_wss.^2 + atlas.meany_wss.^2 + atlas.meanz_wss.^2);

for n = 1:size(PATHNAME,2)
    std_x(:,n) = (WSSx(:,n) - atlas.meanx_wss).^2;
    std_y(:,n) = (WSSy(:,n) - atlas.meany_wss).^2;
    std_z(:,n) = (WSSz(:,n) - atlas.meanz_wss).^2;
end
atlas.stdx_wss = sqrt(sum(std_x,2)./size(PATHNAME,2));
atlas.stdy_wss = sqrt(sum(std_y,2)./size(PATHNAME,2));
atlas.stdz_wss = sqrt(sum(std_z,2)./size(PATHNAME,2));
atlas.std_wss = sqrt(atlas.stdx_wss.^2 + atlas.stdy_wss.^2 + atlas.stdz_wss.^2);

figure('Name','Mean atlas')
a = [8 20];
c = [ ];
[F2,V2,C2]=quiver3Dpatch(geo.x_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),geo.y_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),geo.z_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),...
    atlas.meanx_wss(1:skipfactor:size(geo.x_coor_wss,1),:),atlas.meany_wss(1:skipfactor:size(geo.x_coor_wss,1),:),atlas.meanz_wss(1:skipfactor:size(geo.x_coor_wss,1),:),c,a);
patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
c2=colorbar;caxis([0 1.5])
axis equal;axis off; axis ij
view([0 90])

figure('Name','SD atlas')
patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none', 'FaceVertexCData',atlas.std_wss,'FaceColor','interp','FaceAlpha',1);colorbar;
axis equal;axis off; axis ij
caxis([0 1.5])
view([0 90])

    if calculate_interplation_errorFlag == 1;
       error_matrix = abs(mean_wss_before_interpolation-mean_wss_after_interpolation) ./ ...
           ((mean_wss_before_interpolation+mean_wss_after_interpolation)./2);
       
       disp(['Mean interpolation error inner AAo = ' num2str(mean(error_matrix(:,1),1)*100) ' +/- ' num2str(std(error_matrix(:,1),1)*100) '%'])
       disp(['Mean interpolation error outer AAo = ' num2str(mean(error_matrix(:,2),1)*100) ' +/- ' num2str(std(error_matrix(:,2),1)*100) '%'])
       disp(['Mean interpolation error inner arch = ' num2str(mean(error_matrix(:,3),1)*100) ' +/- ' num2str(std(error_matrix(:,3),1)*100) '%'])
       disp(['Mean interpolation error outer arch = ' num2str(mean(error_matrix(:,4),1)*100) ' +/- ' num2str(std(error_matrix(:,4),1)*100) '%'])
       disp(['Mean interpolation error inner DAo = ' num2str(mean(error_matrix(:,5),1)*100) ' +/- ' num2str(std(error_matrix(:,5),1)*100) '%'])
       disp(['Mean interpolation error outer DAo = ' num2str(mean(error_matrix(:,6),1)*100) ' +/- ' num2str(std(error_matrix(:,6),1)*100) '%'])   
    end

directory = uigetdir('C:\1_Chicago\Data\MIMICS\');
disp('...saving atlas...')
save(strcat(directory,'\atlas'),'atlas')
disp('done')
