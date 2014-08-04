function make_atlas_point_cloud

clc, clear, close all

calculate_interplation_errorFlag = 1;
plotFlag = 1;
offset = 40; % offset needed to prevent negative coordinates after registration

PATHNAME_probability_mask = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\';
probability_mask = 'probability_mask_TEST.mat';
load(strcat(PATHNAME_probability_mask,probability_mask));

%%%% WSS to load
PATHNAME{1} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\1_CAMRI-JEHAM\20140307_084523_Aera_NMH\3dpc\';
PATHNAME{2} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\2_CAMRI-LAJAM\20140313_102146_Aera_NMH\3dpc_nav80\';
% PATHNAME{3} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\3_CAMRI-KIMCO\20140313_115319_Aera_NMH\3dpc_nav80\';
% PATHNAME{4} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\4_CAMRI-THZUM\20140314_073534_Aera_NMH\3dpc_nav80\';
% PATHNAME{5} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\5_CAMRI-MAGIM2\20140314_090452_Aera_NMH\3dpc_nav80\';
% PATHNAME{6} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\6_CAMRI-WIANM\20140317_095816_Aera_NMH\3dpc_nav80\';
% PATHNAME{7} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\7_CAMRI_ARNIF\20140404_082635_Aera_NMH\';
% PATHNAME{8} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\8_CAMRI-FREPAM\20140404_094040_Aera_NMH\3dpc_nav80\';
% PATHNAME{9} = 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\9_CAMRI_JUVIF\20140408_083611_Aera_NMH\3dpc_nav80\';
% PATHNAME{10}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\10_CAMRI_JAKIM\20140414_085925_Aera_NMH\';
% PATHNAME{11}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\11_CAMRI_VIVIM\20140416_081355_Aera_NMH\';
% PATHNAME{12}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\12_CAMRI_RIJO\20140429_080804_Aera_NMH\';
% PATHNAME{13}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\13_CAMRI_MASNF\20140502_075511_Aera_NMH\';
% PATHNAME{14}= 'C:\1_Chicago\Data\MIMICS\18_Test-Retest\14_CAMRI_ANMAM\20140506_081402_Aera_NMH\';
FILENAME = 'data_done';

mask1 = probability_mask.matrix;
mask1_vox = probability_mask.vox;clear probability_mask

%%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
sizes = [size(mask1,1)+offset size(mask1,2)+offset size(mask1,3)+offset];% size(mask1,4) size(mask1,5)];
mask1b = zeros(sizes);
mask1b((offset+1):size(mask1b,1),(offset+1):size(mask1b,2),(offset+1):size(mask1b,3),:,:) = mask1;
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

VELx = zeros(size(geo.x_coor_vel,1),size(PATHNAME,2));
VELy = zeros(size(geo.x_coor_vel,1),size(PATHNAME,2));
VELz = zeros(size(geo.x_coor_vel,1),size(PATHNAME,2));
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
    
    L2a = (mask2 ~= 0); % Using two logicals for the second mask since the size changes after the translation which is needed for coordinates but not for velocity
    
    %%% translate the geometry away from the origin to prevent coordinates < 0 after registration, otherwise the geometry can not be transformed back to a matrix
    sizes = [size(mask2,1)+offset size(mask2,2)+offset size(mask2,3)+offset];% size(mask1,4) size(mask1,5)];
    mask2b = zeros(sizes);
    mask2b((offset+1):size(mask2b,1),(offset+1):size(mask2b,2),(offset+1):size(mask2b,3),:,:) = mask2;
    mask2 = mask2b;clear mask2b
    
    L2b = (mask2 ~= 0);
    % create velocity coordinates
    [x,y,z] = meshgrid((1:size(mask2,2)).* data2.vox(2), ...
        (1:size(mask2,1)).*data2.vox(1),(1:size(mask2,3)).* data2.vox(3));
    data2.x_coor_vel = x(L2b);data2.y_coor_vel = y(L2b);data2.z_coor_vel = z(L2b);
    clear x, clear y, clear z
    contours = zeros(size(L2b));
    contours(L2b==0) = -1;
    contours(L2b==1) = 1;
    [data2.F,V] = isosurface(contours,0); % make a surface from the detected contours
    data2.V = V .* (ones(size(V,1),1) * data2.vox(1:3));
    
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
    
    [I,time] = find(mean_velo==max(mean_velo));
    
    %%% Wall shear stress coordinates for both datasets
    geo.x_coor_wss = geo.V(:,1);
    geo.y_coor_wss = geo.V(:,2);
    geo.z_coor_wss = geo.V(:,3);
    
    data2.x_coor_wss = data2.V_matrix{1}(:,1) + (offset*data2.vox(1));
    data2.y_coor_wss = data2.V_matrix{1}(:,2) + (offset*data2.vox(2));
    data2.z_coor_wss = data2.V_matrix{1}(:,3) + (offset*data2.vox(3));
    
    % velocity averaged for systolic timesteps
    data2.x_value_vel_t1 = data2.PC_unaliased(:,:,:,1,time-2);data2.y_value_vel_t1 = data2.PC_unaliased(:,:,:,2,time-2);data2.z_value_vel_t1 = data2.PC_unaliased(:,:,:,3,time-2);
    data2.x_value_vel_t2 = data2.PC_unaliased(:,:,:,1,time-1);data2.y_value_vel_t2 = data2.PC_unaliased(:,:,:,2,time-1);data2.z_value_vel_t2 = data2.PC_unaliased(:,:,:,3,time-1);
    data2.x_value_vel_t3 = data2.PC_unaliased(:,:,:,1,time);  data2.y_value_vel_t3 = data2.PC_unaliased(:,:,:,2,time);  data2.z_value_vel_t3 = data2.PC_unaliased(:,:,:,3,time);
    data2.x_value_vel_t4 = data2.PC_unaliased(:,:,:,1,time+1);data2.y_value_vel_t4 = data2.PC_unaliased(:,:,:,2,time+1);data2.z_value_vel_t4 = data2.PC_unaliased(:,:,:,3,time+1);
    data2.x_value_vel_t5 = data2.PC_unaliased(:,:,:,1,time+2);data2.y_value_vel_t5 = data2.PC_unaliased(:,:,:,2,time+2);data2.z_value_vel_t5 = data2.PC_unaliased(:,:,:,3,time+2);
    
    data2.x_value_vel = (data2.x_value_vel_t1(L2a) + data2.x_value_vel_t2(L2a) + data2.x_value_vel_t3(L2a) + data2.x_value_vel_t4(L2a) + data2.x_value_vel_t5(L2a))./5;
    data2.y_value_vel = (data2.y_value_vel_t1(L2a) + data2.y_value_vel_t2(L2a) + data2.y_value_vel_t3(L2a) + data2.y_value_vel_t4(L2a) + data2.y_value_vel_t5(L2a))./5;
    data2.z_value_vel = (data2.z_value_vel_t1(L2a) + data2.z_value_vel_t2(L2a) + data2.z_value_vel_t3(L2a) + data2.z_value_vel_t4(L2a) + data2.z_value_vel_t5(L2a))./5;
    
    data2.vel_m = sqrt(data2.x_value_vel.^2 + data2.y_value_vel.^2 + data2.z_value_vel.^2);
    
    % WSS averaged for systolic timesteps
    data2.x_value_wss_t1 = data2.WSS_matrix{time-2}(:,1);data2.y_value_wss_t1 = data2.WSS_matrix{time-2}(:,2);data2.z_value_wss_t1 = data2.WSS_matrix{time-2}(:,3);
    data2.x_value_wss_t2 = data2.WSS_matrix{time-1}(:,1);data2.y_value_wss_t2 = data2.WSS_matrix{time-1}(:,2);data2.z_value_wss_t2 = data2.WSS_matrix{time-1}(:,3);
    data2.x_value_wss_t3 = data2.WSS_matrix{time}(:,1);  data2.y_value_wss_t3 = data2.WSS_matrix{time}(:,2);  data2.z_value_wss_t3 = data2.WSS_matrix{time}(:,3);
    data2.x_value_wss_t4 = data2.WSS_matrix{time+1}(:,1);data2.y_value_wss_t4 = data2.WSS_matrix{time+1}(:,2);data2.z_value_wss_t4 = data2.WSS_matrix{time+1}(:,3);
    data2.x_value_wss_t5 = data2.WSS_matrix{time+2}(:,1);data2.y_value_wss_t5 = data2.WSS_matrix{time+2}(:,2);data2.z_value_wss_t5 = data2.WSS_matrix{time+2}(:,3);
    
    data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4 + data2.x_value_wss_t5)./5;
    data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4 + data2.y_value_wss_t5)./5;
    data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4 + data2.z_value_wss_t5)./5;
    
    data2.wss_m = sqrt(data2.x_value_wss.^2 + data2.y_value_wss.^2 + data2.z_value_wss.^2);
    
    if plotFlag == 1
        figure('Name','data2 velocity')
        a = [8 20];
        c = [ ];
        skipfactor = 4;
        patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss], ...
            'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.1);
        hold on
        [F2,V2,C2]=quiver3Dpatch(data2.x_coor_vel(1:skipfactor:size(data2.x_coor_vel,1)),data2.y_coor_vel(1:skipfactor:size(data2.x_coor_vel,1)),data2.z_coor_vel(1:skipfactor:size(data2.x_coor_vel,1)),data2.x_value_vel(1:skipfactor:size(data2.x_coor_vel,1)) ...
            ,data2.y_value_vel(1:skipfactor:size(data2.x_coor_vel,1)),data2.z_value_vel(1:skipfactor:size(data2.x_coor_vel,1)),c,a);
        patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
        colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
        
        figure('Name','data2 WSS')
        a = [8 20];
        c = [ ];
        skipfactor = 4;
        patch('Faces',data2.F_matrix{1},'Vertices',[data2.x_coor_wss data2.y_coor_wss data2.z_coor_wss], ...
            'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
        hold on
        [F2,V2,C2]=quiver3Dpatch(data2.x_coor_wss(1:skipfactor:size(data2.x_coor_wss,1)),data2.y_coor_wss(1:skipfactor:size(data2.x_coor_wss,1)),data2.z_coor_wss(1:skipfactor:size(data2.x_coor_wss,1)),data2.x_value_wss(1:skipfactor:size(data2.x_coor_wss,1)) ...
            ,data2.y_value_wss(1:skipfactor:size(data2.x_coor_wss,1)),data2.z_value_wss(1:skipfactor:size(data2.x_coor_wss,1)),c,a);
        patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
        colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    end
    
    if calculate_interplation_errorFlag == 1;
        %         wss_magn_new = sqrt(data2.x_value_wss.^2+data2.y_value_wss.^2+data2.z_value_wss.^2);
        %         wss_magn_new = sqrt(data2.x_value_wss.^2+data2.y_value_wss.^2+data2.z_value_wss.^2);
        if ~exist(strcat(PATHNAME{n},'mask6_before_registration\mask1.mat'),'file')
            for i = 1:6
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo)
                region = getPosition(polyAAo);
                %mask = inpolygon(x, y, AAo(:,1), AAo(:,2));
                
                disp('saving, pausing')
                mkdir(PATHNAME{n},'mask6_before_registration')
                save(strcat([PATHNAME{n} 'mask6_before_registration\mask' num2str(i)]),'region');
                pause
            end
        end
        load(strcat(PATHNAME{n},'mask6_before_registration\mask1'))
        data2_mask_AAo_inner_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_inner = mean(data2.vel_m(data2_mask_AAo_inner_vel));
        data2_mask_AAo_inner_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_inner = mean(data2.wss_m(data2_mask_AAo_inner_wss));
        load(strcat(PATHNAME{n},'mask6_before_registration\mask2'))
        data2_mask_AAo_outer_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_outer = mean(data2.vel_m(data2_mask_AAo_outer_vel));
        data2_mask_AAo_outer_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_outer = mean(data2.wss_m(data2_mask_AAo_outer_wss));
        load(strcat(PATHNAME{n},'mask6_before_registration\mask3'))
        data2_mask_arch_inner_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_inner = mean(data2.vel_m(data2_mask_arch_inner_vel));
        data2_mask_arch_inner_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_inner = mean(data2.wss_m(data2_mask_arch_inner_wss));
        load(strcat(PATHNAME{n},'mask6_before_registration\mask4'))
        data2_mask_arch_outer_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_outer = mean(data2.vel_m(data2_mask_arch_outer_vel));
        data2_mask_arch_outer_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_outer = mean(data2.wss_m(data2_mask_arch_outer_wss));
        load(strcat(PATHNAME{n},'mask6_before_registration\mask5'))
        data2_mask_DAo_inner_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_inner = mean(data2.vel_m(data2_mask_DAo_inner_vel));
        data2_mask_DAo_inner_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_inner = mean(data2.wss_m(data2_mask_DAo_inner_wss));
        load(strcat(PATHNAME{n},'mask6_before_registration\mask6'))
        data2_mask_DAo_outer_vel = inpolygon(data2.x_coor_vel, data2.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_outer = mean(data2.vel_m(data2_mask_DAo_outer_vel));
        data2_mask_DAo_outer_wss = inpolygon(data2.x_coor_wss, data2.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_outer = mean(data2.wss_m(data2_mask_DAo_outer_wss));
        
        mean_vel_before_interpolation(n,1) = mean_vel_asc_inner;
        mean_vel_before_interpolation(n,2) = mean_vel_asc_outer;
        mean_vel_before_interpolation(n,3) = mean_vel_arch_inner;
        mean_vel_before_interpolation(n,4) = mean_vel_arch_outer;
        mean_vel_before_interpolation(n,5) = mean_vel_DAo_inner;
        mean_vel_before_interpolation(n,6) = mean_vel_DAo_outer;
        
        mean_wss_before_interpolation(n,1) = mean_wss_asc_inner;
        mean_wss_before_interpolation(n,2) = mean_wss_asc_outer;
        mean_wss_before_interpolation(n,3) = mean_wss_arch_inner;
        mean_wss_before_interpolation(n,4) = mean_wss_arch_outer;
        mean_wss_before_interpolation(n,5) = mean_wss_DAo_inner;
        mean_wss_before_interpolation(n,6) = mean_wss_DAo_outer;
        
        if plotFlag == 1
            figure(95)
            scatter3(data2.x_coor_vel(data2_mask_AAo_inner_vel),data2.y_coor_vel(data2_mask_AAo_inner_vel),data2.z_coor_vel(data2_mask_AAo_inner_vel),20,data2.vel_m(data2_mask_AAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(96)
            scatter3(data2.x_coor_vel(data2_mask_AAo_outer_vel),data2.y_coor_vel(data2_mask_AAo_outer_vel),data2.z_coor_vel(data2_mask_AAo_outer_vel),20,data2.vel_m(data2_mask_AAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(97)
            scatter3(data2.x_coor_vel(data2_mask_arch_inner_vel),data2.y_coor_vel(data2_mask_arch_inner_vel),data2.z_coor_vel(data2_mask_arch_inner_vel),20,data2.vel_m(data2_mask_arch_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(98)
            scatter3(data2.x_coor_vel(data2_mask_arch_outer_vel),data2.y_coor_vel(data2_mask_arch_outer_vel),data2.z_coor_vel(data2_mask_arch_outer_vel),20,data2.vel_m(data2_mask_arch_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(99)
            scatter3(data2.x_coor_vel(data2_mask_DAo_inner_vel),data2.y_coor_vel(data2_mask_DAo_inner_vel),data2.z_coor_vel(data2_mask_DAo_inner_vel),20,data2.vel_m(data2_mask_DAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(100)
            scatter3(data2.x_coor_vel(data2_mask_DAo_outer_vel),data2.y_coor_vel(data2_mask_DAo_outer_vel),data2.z_coor_vel(data2_mask_DAo_outer_vel),20,data2.vel_m(data2_mask_DAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(101)
            scatter3(data2.x_coor_wss(data2_mask_AAo_inner_wss),data2.y_coor_wss(data2_mask_AAo_inner_wss),data2.z_coor_wss(data2_mask_AAo_inner_wss),20,data2.wss_m(data2_mask_AAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(102)
            scatter3(data2.x_coor_wss(data2_mask_AAo_outer_wss),data2.y_coor_wss(data2_mask_AAo_outer_wss),data2.z_coor_wss(data2_mask_AAo_outer_wss),20,data2.wss_m(data2_mask_AAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(103)
            scatter3(data2.x_coor_wss(data2_mask_arch_inner_wss),data2.y_coor_wss(data2_mask_arch_inner_wss),data2.z_coor_wss(data2_mask_arch_inner_wss),20,data2.wss_m(data2_mask_arch_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(104)
            scatter3(data2.x_coor_wss(data2_mask_arch_outer_wss),data2.y_coor_wss(data2_mask_arch_outer_wss),data2.z_coor_wss(data2_mask_arch_outer_wss),20,data2.wss_m(data2_mask_arch_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(105)
            scatter3(data2.x_coor_wss(data2_mask_DAo_inner_wss),data2.y_coor_wss(data2_mask_DAo_inner_wss),data2.z_coor_wss(data2_mask_DAo_inner_wss),20,data2.wss_m(data2_mask_DAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
            figure(106)
            scatter3(data2.x_coor_wss(data2_mask_DAo_outer_wss),data2.y_coor_wss(data2_mask_DAo_outer_wss),data2.z_coor_wss(data2_mask_DAo_outer_wss),20,data2.wss_m(data2_mask_DAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');axis ij;view([-180 -90]);
        end
    end
    
    %     figure('Name','mask2')
    %     V = vol3d('cdata',L,'texture','3D','texturemap',L);
    %     xlabel('x');ylabel('y');zlabel('z')
    %     axis tight; axis equal;axis ij
    
    if plotFlag == 1
        figure('Name',strcat('To be registered',num2str(n)))
        plot3(geo.x_coor_vel,geo.y_coor_vel,geo.z_coor_vel,'r.')
        hold on
        % plot3(geo.x_coor_wss,geo.y_coor_wss,geo.z_coor_wss,'g.')
        plot3(data2.x_coor_vel,data2.y_coor_vel,data2.z_coor_vel,'b.')
        % plot3(data2.x_coor_wss,data2.y_coor_wss,data2.z_coor_wss,'k.')
        legend('to remain the same','to be transformed')
        axis equal; axis off;view([-180 -90]); axis ij
    end
    
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
    
    disp('rigid registration (dof = 6)')
    flirtcmd= [...
        'flirt -searchrx -0 0 -searchry -0 0 -searchrz -0 0 -dof 6' ...
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
    yxz_coor_vel = [data2.y_coor_vel data2.x_coor_vel data2.z_coor_vel];clear x_coor_vel2, clear y_coor_vel2, clear z_coor_vel2
    yxz_coor_vel(:,4) = 1;
    yxz_coor_vel_new = inv(worldmat)*yxz_coor_vel'; clear yxz_coor_vel
    x_coor_vel = yxz_coor_vel_new(2,:)';
    y_coor_vel = yxz_coor_vel_new(1,:)';
    z_coor_vel = yxz_coor_vel_new(3,:)';clear yxz_coor_vel_new
    
    %%% Velocity values
    yxz_value_vel = [data2.y_value_vel data2.x_value_vel data2.z_value_vel];
    yxz_value_vel_new = inv(rotmat)*yxz_value_vel'; clear yxz_value_vel
    x_value_vel = yxz_value_vel_new(2,:)';
    y_value_vel = yxz_value_vel_new(1,:)';
    z_value_vel = yxz_value_vel_new(3,:)';clear yxz_value_vel_new
    
    %%% WSS coordinates
    yxz_coor_wss = [data2.y_coor_wss data2.x_coor_wss data2.z_coor_wss];clear x_coor_wss2, clear y_coor_wss2, clear z_coor_wss2
    yxz_coor_wss(:,4) = 1;
    yxz_coor_wss_new = inv(worldmat)*yxz_coor_wss'; clear yxz_coor_wss
    x_coor_wss = yxz_coor_wss_new(2,:)';
    y_coor_wss = yxz_coor_wss_new(1,:)';
    z_coor_wss = yxz_coor_wss_new(3,:)';clear yxz_coor_wss_new
    
    %%% WSS values
    yxz_value_wss = [data2.y_value_wss data2.x_value_wss data2.z_value_wss];
    yxz_value_wss_new = inv(rotmat)*yxz_value_wss'; clear yxz_value_wss
    x_value_wss = yxz_value_wss_new(2,:)';
    y_value_wss = yxz_value_wss_new(1,:)';
    z_value_wss = yxz_value_wss_new(3,:)'; clear yxz_value_wss_new
    
    if plotFlag == 1
        figure('Name','Registered')
        plot3(geo.x_coor_vel,geo.y_coor_vel,geo.z_coor_vel,'r.')
        hold on
        %  plot3(geo.x_coor_wss,geo.y_coor_wss,geo.z_coor_wss,'g.')
        plot3(x_coor_vel,y_coor_vel,z_coor_vel,'b.')
        %  plot3(x_coor_wss,y_coor_wss,z_coor_wss,'k.')
        legend('Prob mask','Aorta')
        %legend('probability mask/atlas','individual aorta')
        axis equal; axis off;view([-180 -90]); axis ij
        pause(10)
        
        figure('Name','transformed velocity')
        skipfactor = 4;
        a = [8 20];
        c = [ ];
        [F2,V2,C2]=quiver3Dpatch(x_coor_vel(1:skipfactor:size(x_coor_vel,1),:),y_coor_vel(1:skipfactor:size(x_coor_vel,1),:),z_coor_vel(1:skipfactor:size(x_coor_vel,1),:), ...
            x_value_vel(1:skipfactor:size(x_coor_vel,1),:),y_value_vel(1:skipfactor:size(x_coor_vel,1),:),z_value_vel(1:skipfactor:size(x_coor_vel,1),:),c,a);
        patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
        c2=colorbar;caxis([0 1.5])
        axis equal;axis off; axis ij
        view([-180 -90])
        pause(10)
        
        figure('Name','transformed WSS')
        skipfactor = 4;
        a = [8 20];
        c = [ ];
        [F2,V2,C2]=quiver3Dpatch(x_coor_wss(1:skipfactor:size(x_coor_wss,1),:),y_coor_wss(1:skipfactor:size(x_coor_wss,1),:),z_coor_wss(1:skipfactor:size(x_coor_wss,1),:), ...
            x_value_wss(1:skipfactor:size(x_coor_wss,1),:),y_value_wss(1:skipfactor:size(x_coor_wss,1),:),z_value_wss(1:skipfactor:size(x_coor_wss,1),:),c,a);
        patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
        c2=colorbar;caxis([0 1.5])
        axis equal;axis off; axis ij
        view([-180 -90])
        pause(10)
    end
    
    interpolation_function = TriScatteredInterp([x_coor_vel y_coor_vel z_coor_vel],x_value_vel,'nearest');
    x_value_vel = interpolation_function([geo.x_coor_vel geo.y_coor_vel geo.z_coor_vel]);
    x_value_vel(isnan(x_value_vel)) = 0;
    
    interpolation_function = TriScatteredInterp([x_coor_vel y_coor_vel z_coor_vel],y_value_vel,'nearest');
    y_value_vel = interpolation_function([geo.x_coor_vel geo.y_coor_vel geo.z_coor_vel]);
    y_value_vel(isnan(y_value_vel)) = 0;
    
    interpolation_function = TriScatteredInterp([x_coor_vel y_coor_vel z_coor_vel],z_value_vel,'nearest');
    z_value_vel = interpolation_function([geo.x_coor_vel geo.y_coor_vel geo.z_coor_vel]);
    z_value_vel(isnan(z_value_vel)) = 0;
    
    interpolation_function = TriScatteredInterp([x_coor_wss y_coor_wss z_coor_wss],x_value_wss,'nearest');
    x_value_wss = interpolation_function([geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss]);
    x_value_wss(isnan(x_value_wss)) = 0;
    
    interpolation_function = TriScatteredInterp([x_coor_wss y_coor_wss z_coor_wss],y_value_wss,'nearest');
    y_value_wss = interpolation_function([geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss]);
    y_value_wss(isnan(y_value_wss)) = 0;
    
    interpolation_function = TriScatteredInterp([x_coor_wss y_coor_wss z_coor_wss],z_value_wss,'nearest');
    z_value_wss = interpolation_function([geo.x_coor_wss geo.y_coor_wss geo.z_coor_wss]);
    z_value_wss(isnan(z_value_wss)) = 0;
    
    if plotFlag == 1
        figure('Name','interpolated to atlas velocity')
        skipfactor = 4;
        a = [8 20];
        c = [ ];
        [F2,V2,C2]=quiver3Dpatch(geo.x_coor_vel(1:skipfactor:size(geo.x_coor_vel,1),:),geo.y_coor_vel(1:skipfactor:size(geo.x_coor_vel,1),:),geo.z_coor_vel(1:skipfactor:size(geo.x_coor_vel,1),:), ...
            x_value_vel(1:skipfactor:size(geo.x_coor_vel,1),:),y_value_vel(1:skipfactor:size(geo.x_coor_vel,1),:),z_value_vel(1:skipfactor:size(geo.x_coor_vel,1),:),c,a);
        patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
        c2=colorbar;caxis([0 1.5])
        axis equal;axis off; axis ij
        view([-180 -90])
        pause(10)
        
        figure('Name','interpolated to atlas WSS')
        skipfactor = 4;
        a = [8 20];
        c = [ ];
        [F2,V2,C2]=quiver3Dpatch(geo.x_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),geo.y_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),geo.z_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:), ...
            x_value_wss(1:skipfactor:size(geo.x_coor_wss,1),:),y_value_wss(1:skipfactor:size(geo.x_coor_wss,1),:),z_value_wss(1:skipfactor:size(geo.x_coor_wss,1),:),c,a);
        patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
        c2=colorbar;caxis([0 1.5])
        axis equal;axis off; axis ij
        view([-180 -90])
        pause(10)
    end
    
    if calculate_interplation_errorFlag == 1;
        data2.vel_m = sqrt(x_value_vel.^2+y_value_vel.^2+z_value_vel.^2);
        data2.wss_m = sqrt(x_value_wss.^2+y_value_wss.^2+z_value_wss.^2);
        if ~exist(strcat(PATHNAME_probability_mask,'mask6\mask1.mat'),'file')
            for i = 1:6
                %Polygon and mask for AAo
                polyAAo = impoly;
                wait(polyAAo)
                region = getPosition(polyAAo);
                %mask = inpolygon(x, y, AAo(:,1), AAo(:,2));
                
                disp('saving, pausing')
                mkdir(PATHNAME_probability_mask,'mask6')
                save(strcat([PATHNAME_probability_mask 'mask6\mask' num2str(i)]),'region');
                pause
            end
        end
        load(strcat(PATHNAME_probability_mask,'mask6\mask1'))
        geo_mask_AAo_inner_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_inner = mean(data2.vel_m(geo_mask_AAo_inner_vel));
        geo_mask_AAo_inner_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_inner = mean(data2.wss_m(geo_mask_AAo_inner_wss));
        load(strcat(PATHNAME_probability_mask,'mask6\mask2'))
        geo_mask_AAo_outer_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_asc_outer = mean(data2.vel_m(geo_mask_AAo_outer_vel));
        geo_mask_AAo_outer_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_asc_outer = mean(data2.wss_m(geo_mask_AAo_outer_wss));
        load(strcat(PATHNAME_probability_mask,'mask6\mask3'))
        geo_mask_arch_inner_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_inner = mean(data2.vel_m(geo_mask_arch_inner_vel));
        geo_mask_arch_inner_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_inner = mean(data2.wss_m(geo_mask_arch_inner_wss));
        load(strcat(PATHNAME_probability_mask,'mask6\mask4'))
        geo_mask_arch_outer_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_arch_outer = mean(data2.vel_m(geo_mask_arch_outer_vel));
        geo_mask_arch_outer_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_arch_outer = mean(data2.wss_m(geo_mask_arch_outer_wss));
        load(strcat(PATHNAME_probability_mask,'mask6\mask5'))
        geo_mask_DAo_inner_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_inner = mean(data2.vel_m(geo_mask_DAo_inner_vel));
        geo_mask_DAo_inner_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_inner = mean(data2.wss_m(geo_mask_DAo_inner_wss));
        load(strcat(PATHNAME_probability_mask,'mask6\mask6'))
        geo_mask_DAo_outer_vel = inpolygon(geo.x_coor_vel, geo.y_coor_vel, region(:,1), region(:,2));
        mean_vel_DAo_outer = mean(data2.vel_m(geo_mask_DAo_outer_vel));
        geo_mask_DAo_outer_wss = inpolygon(geo.x_coor_wss, geo.y_coor_wss, region(:,1), region(:,2));
        mean_wss_DAo_outer = mean(data2.wss_m(geo_mask_DAo_outer_wss));
        
        if plotFlag == 1
            figure(95)
            scatter3(geo.x_coor_vel(geo_mask_AAo_inner_vel),geo.y_coor_vel(geo_mask_AAo_inner_vel),geo.z_coor_vel(geo_mask_AAo_inner_vel),20,data2.vel_m(geo_mask_AAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(96)
            scatter3(geo.x_coor_vel(geo_mask_AAo_outer_vel),geo.y_coor_vel(geo_mask_AAo_outer_vel),geo.z_coor_vel(geo_mask_AAo_outer_vel),20,data2.vel_m(geo_mask_AAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(97)
            scatter3(geo.x_coor_vel(geo_mask_arch_inner_vel),geo.y_coor_vel(geo_mask_arch_inner_vel),geo.z_coor_vel(geo_mask_arch_inner_vel),20,data2.vel_m(geo_mask_arch_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(98)
            scatter3(geo.x_coor_vel(geo_mask_arch_outer_vel),geo.y_coor_vel(geo_mask_arch_outer_vel),geo.z_coor_vel(geo_mask_arch_outer_vel),20,data2.vel_m(geo_mask_arch_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(99)
            scatter3(geo.x_coor_vel(geo_mask_DAo_inner_vel),geo.y_coor_vel(geo_mask_DAo_inner_vel),geo.z_coor_vel(geo_mask_DAo_inner_vel),20,data2.vel_m(geo_mask_DAo_inner_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(100)
            scatter3(geo.x_coor_vel(geo_mask_DAo_outer_vel),geo.y_coor_vel(geo_mask_DAo_outer_vel),geo.z_coor_vel(geo_mask_DAo_outer_vel),20,data2.vel_m(geo_mask_DAo_outer_vel),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(101)
            scatter3(geo.x_coor_wss(geo_mask_AAo_inner_wss),geo.y_coor_wss(geo_mask_AAo_inner_wss),geo.z_coor_wss(geo_mask_AAo_inner_wss),20,data2.wss_m(geo_mask_AAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(102)
            scatter3(geo.x_coor_wss(geo_mask_AAo_outer_wss),geo.y_coor_wss(geo_mask_AAo_outer_wss),geo.z_coor_wss(geo_mask_AAo_outer_wss),20,data2.wss_m(geo_mask_AAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(103)
            scatter3(geo.x_coor_wss(geo_mask_arch_inner_wss),geo.y_coor_wss(geo_mask_arch_inner_wss),geo.z_coor_wss(geo_mask_arch_inner_wss),20,data2.wss_m(geo_mask_arch_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(104)
            scatter3(geo.x_coor_wss(geo_mask_arch_outer_wss),geo.y_coor_wss(geo_mask_arch_outer_wss),geo.z_coor_wss(geo_mask_arch_outer_wss),20,data2.wss_m(geo_mask_arch_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(105)
            scatter3(geo.x_coor_wss(geo_mask_DAo_inner_wss),geo.y_coor_wss(geo_mask_DAo_inner_wss),geo.z_coor_wss(geo_mask_DAo_inner_wss),20,data2.wss_m(geo_mask_DAo_inner_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
            figure(106)
            scatter3(geo.x_coor_wss(geo_mask_DAo_outer_wss),geo.y_coor_wss(geo_mask_DAo_outer_wss),geo.z_coor_wss(geo_mask_DAo_outer_wss),20,data2.wss_m(geo_mask_DAo_outer_wss),'filled');axis equal;caxis([0 1.5])
            xlabel('x'),ylabel('y'),zlabel('z');view([0 -90]);
        end
        
        mean_vel_after_interpolation(n,1) = mean_vel_asc_inner;
        mean_vel_after_interpolation(n,2) = mean_vel_asc_outer;
        mean_vel_after_interpolation(n,3) = mean_vel_arch_inner;
        mean_vel_after_interpolation(n,4) = mean_vel_arch_outer;
        mean_vel_after_interpolation(n,5) = mean_vel_DAo_inner;
        mean_vel_after_interpolation(n,6) = mean_vel_DAo_outer;
        mean_wss_after_interpolation(n,1) = mean_wss_asc_inner;
        mean_wss_after_interpolation(n,2) = mean_wss_asc_outer;
        mean_wss_after_interpolation(n,3) = mean_wss_arch_inner;
        mean_wss_after_interpolation(n,4) = mean_wss_arch_outer;
        mean_wss_after_interpolation(n,5) = mean_wss_DAo_inner;
        mean_wss_after_interpolation(n,6) = mean_wss_DAo_outer;
        
        mean_vel_before_interpolation
        mean_vel_after_interpolation
        mean_wss_before_interpolation
        mean_wss_after_interpolation
    end
    
    pause(10)
    close all
    
    VELx(:,n) = x_value_vel;
    VELy(:,n) = y_value_vel;
    VELz(:,n) = z_value_vel;
    WSSx(:,n) = x_value_wss;
    WSSy(:,n) = y_value_wss;
    WSSz(:,n) = z_value_wss;
end

atlas.x_coor_vel = geo.x_coor_vel;
atlas.y_coor_vel = geo.y_coor_vel;
atlas.z_coor_vel = geo.z_coor_vel;
atlas.meanx_vel = sum(VELx,2)./size(PATHNAME,2);
atlas.meany_vel = sum(VELy,2)./size(PATHNAME,2);
atlas.meanz_vel = sum(VELz,2)./size(PATHNAME,2);
atlas.x_coor_wss = geo.x_coor_wss;
atlas.y_coor_wss = geo.y_coor_wss;
atlas.z_coor_wss = geo.z_coor_wss;
atlas.meanx_wss = sum(WSSx,2)./size(PATHNAME,2);
atlas.meany_wss = sum(WSSy,2)./size(PATHNAME,2);
atlas.meanz_wss = sum(WSSz,2)./size(PATHNAME,2);
atlas.vox = mask1_vox;
atlas.faces = geo.F;
atlas.vertices = geo.V;
atlas.mask = mask1;
atlas.meanvel = sqrt(atlas.meanx_vel.^2 + atlas.meany_vel.^2 + atlas.meanz_vel.^2);
atlas.mean_wss = sqrt(atlas.meanx_wss.^2 + atlas.meany_wss.^2 + atlas.meanz_wss.^2);

for n = 1:size(PATHNAME,2)
    std_x_vel(:,n) = (VELx(:,n) - atlas.meanx_vel).^2;
    std_y_vel(:,n) = (VELy(:,n) - atlas.meany_vel).^2;
    std_z_vel(:,n) = (VELz(:,n) - atlas.meanz_vel).^2;
    std_x_wss(:,n) = (WSSx(:,n) - atlas.meanx_wss).^2;
    std_y_wss(:,n) = (WSSy(:,n) - atlas.meany_wss).^2;
    std_z_wss(:,n) = (WSSz(:,n) - atlas.meanz_wss).^2;
end
atlas.stdx_vel = sqrt(sum(std_x_vel,2)./size(PATHNAME,2));
atlas.stdy_vel = sqrt(sum(std_y_vel,2)./size(PATHNAME,2));
atlas.stdz_vel = sqrt(sum(std_z_vel,2)./size(PATHNAME,2));
atlas.std_vel = sqrt(atlas.stdx_vel.^2 + atlas.stdy_vel.^2 + atlas.stdz_vel.^2);
atlas.stdx_wss = sqrt(sum(std_x_wss,2)./size(PATHNAME,2));
atlas.stdy_wss = sqrt(sum(std_y_wss,2)./size(PATHNAME,2));
atlas.stdz_wss = sqrt(sum(std_z_wss,2)./size(PATHNAME,2));
atlas.std_wss = sqrt(atlas.stdx_wss.^2 + atlas.stdy_wss.^2 + atlas.stdz_wss.^2);

if plotFlag == 1
    figure('Name','Mean atlas velocity')
    a = [8 20];
    c = [ ];
    [F2,V2,C2]=quiver3Dpatch(atlas.x_coor_vel(1:skipfactor:size(atlas.x_coor_vel,1),:),atlas.y_coor_vel(1:skipfactor:size(atlas.x_coor_vel,1),:),atlas.z_coor_vel(1:skipfactor:size(atlas.x_coor_vel,1),:),...
        atlas.meanx_vel(1:skipfactor:size(atlas.x_coor_vel,1),:),atlas.meany_vel(1:skipfactor:size(atlas.x_coor_vel,1),:),atlas.meanz_vel(1:skipfactor:size(atlas.x_coor_vel,1),:),c,a);
    patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
    colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    
    figure('Name','Mean atlas WSS')
    a = [8 20];
    c = [ ];
    [F2,V2,C2]=quiver3Dpatch(geo.x_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),geo.y_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),geo.z_coor_wss(1:skipfactor:size(geo.x_coor_wss,1),:),...
        atlas.meanx_wss(1:skipfactor:size(geo.x_coor_wss,1),:),atlas.meany_wss(1:skipfactor:size(geo.x_coor_wss,1),:),atlas.meanz_wss(1:skipfactor:size(geo.x_coor_wss,1),:),c,a);
    patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
    c2=colorbar;caxis([0 1.5]);axis equal;axis off; axis ij;view([-180 -90])
    
    figure('Name','SD atlas velocity')
    patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none', 'FaceVertexCData',atlas.std_vel,'FaceColor','interp','FaceAlpha',0.1);colorbar;
    axis equal;axis off; axis ij;caxis([0 1.5]);view([-180 -90])
    
    figure('Name','SD atlas wss')
    patch('Faces',geo.F,'Vertices',geo.V,'EdgeColor','none', 'FaceVertexCData',atlas.std_wss,'FaceColor','interp','FaceAlpha',1);colorbar;
    axis equal;axis off; axis ij;caxis([0 1.5]);view([-180 -90])
end

if calculate_interplation_errorFlag == 1;
    error_matrix_vel = abs(mean_vel_before_interpolation-mean_vel_after_interpolation) ./ ...
        ((mean_vel_before_interpolation+mean_vel_after_interpolation)./2);
    error_matrix_wss = abs(mean_wss_before_interpolation-mean_wss_after_interpolation) ./ ...
        ((mean_wss_before_interpolation+mean_wss_after_interpolation)./2);
    
    disp(['Velocity: Mean interpolation error inner AAo = ' num2str(mean(error_matrix_vel(:,1),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,1),1)*100) '%'])
    disp(['Velocity: Mean interpolation error outer AAo = ' num2str(mean(error_matrix_vel(:,2),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,2),1)*100) '%'])
    disp(['Velocity: Mean interpolation error inner arch = ' num2str(mean(error_matrix_vel(:,3),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,3),1)*100) '%'])
    disp(['Velocity: Mean interpolation error outer arch = ' num2str(mean(error_matrix_vel(:,4),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,4),1)*100) '%'])
    disp(['Velocity: Mean interpolation error inner DAo = ' num2str(mean(error_matrix_vel(:,5),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,5),1)*100) '%'])
    disp(['Velocity: Mean interpolation error outer DAo = ' num2str(mean(error_matrix_vel(:,6),1)*100) ' +/- ' num2str(std(error_matrix_vel(:,6),1)*100) '%'])
    disp(' ')
    disp(['WSS: Mean interpolation error inner AAo = ' num2str(mean(error_matrix_wss(:,1),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,1),1)*100) '%'])
    disp(['WSS: Mean interpolation error outer AAo = ' num2str(mean(error_matrix_wss(:,2),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,2),1)*100) '%'])
    disp(['WSS: Mean interpolation error inner arch = ' num2str(mean(error_matrix_wss(:,3),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,3),1)*100) '%'])
    disp(['WSS: Mean interpolation error outer arch = ' num2str(mean(error_matrix_wss(:,4),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,4),1)*100) '%'])
    disp(['WSS: Mean interpolation error inner DAo = ' num2str(mean(error_matrix_wss(:,5),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,5),1)*100) '%'])
    disp(['WSS: Mean interpolation error outer DAo = ' num2str(mean(error_matrix_wss(:,6),1)*100) ' +/- ' num2str(std(error_matrix_wss(:,6),1)*100) '%'])
end

directory = uigetdir('C:\1_Chicago\Data\MIMICS\');
disp('...saving atlas...')
save(strcat(directory,'\atlas'),'atlas')
disp('done')