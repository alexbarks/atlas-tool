clc, clear, close all

load d:\research\MR_GT_WSSvalidation\phantom_tube\pulsatile\70BPM_3LMIN\mrstruct_subject_20150324_user009\mask_struct_aorta.mat
mask2 = mrstruct_mask.dataAy;
mask2_vox = mrstruct_mask.vox;

L2 = (mask2 ~= 0);
contours = zeros(size(L2));
contours(L2==0) = -1;
contours(L2==1) = 1;
[F,V] = isosurface(contours,0); % make a surface from the detected contours
V = V .* (ones(size(V,1),1) * mask2_vox(1:3));

load d:\research\MR_GT_WSSvalidation\phantom_tube\pulsatile\70BPM_3LMIN\mrstruct_subject_20150324_user009\Wss_point_cloud_aorta
WSS_all = Wss_point_cloud; clear Wss_point_cloud
WSS = WSS_all{3};

figure('Name','data2 WSS vectors')
a = [2 15];
c = [ ];
patch('Faces',F,'Vertices',V, ...
    'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
hold on
[F2,V2,C2]=quiver3Dpatch(V(:,1),V(:,2),V(:,3),WSS(:,1),WSS(:,2),WSS(:,3),c,a);
patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
c2=colorbar;caxis([0 1.5])
axis equal;axis off; axis ij
view([-180 -90])