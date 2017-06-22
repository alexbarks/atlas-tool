function view_WSS_vectors(TimeFlag)

global mrstruct_mask
global Wss_point_cloud
global mrStruct

clc%, clear%, close all
currDir = pwd;
[FileName,MrstructPath,FilterIndex] = uigetfile(currDir,'Select the mag_struct file of your patient');
try
    cd(MrstructPath);
catch
    warndlg('File not found!');
    return;
end

maskStructFiles=ls('mask_struct_*');
if size(maskStructFiles,1)>1
    [FileName,FilePath,FilterIndex] = uigetfile(MrstructPath,'There are several ''mask_struct'' files, please select the right one');
    load(FileName);        % Load mask
else
    if (exist('mask_struct_aorta.mat') == 2)
        load mask_struct_aorta.mat
    else
        folders = ls;
        for i=3:size(folders,1)
            [a,b]=find(folders(i,1:12)=='mask_struct_');
            if sum(a)==12
                load(folders(i,:));
                break
            end
        end
        if ~(exist(folders(i,:)) == 2)
            [all_name, all_path] = uigetfile('*.mat','Start by loading the mask_struct_ file','Multiselect','Off');
            load(strcat(all_path,all_name));
        end
        clear folders a all_path all_name
    end
end

WssPtCloudFiles=ls('Wss_point_cloud_*');
if size(WssPtCloudFiles,1)>1
    [FileName,FilePath,FilterIndex] = uigetfile(MrstructPath,'There are several ''Wss_point_cloud'' files, please select the right one');
    load(FileName);        % Load WSS
else
    if (exist('Wss_point_cloud_aorta.mat') == 2)
        load Wss_point_cloud_aorta.mat
    else
        folders = ls;
        for i=3:size(folders,1)
            [a,b]=find(folders(i,1:16)=='Wss_point_cloud_');
            if sum(a)==16
                load(folders(i,:));
                break
            end
        end
        if ~(exist(folders(i,:)) == 2)
            [all_name, all_path] = uigetfile('*.mat','Please load the Wss_point_cloud_ file','Multiselect','Off');
            load(strcat(all_path,all_name));
        end
    end
end

mask2 = mrstruct_mask.dataAy;
mask2_vox = mrstruct_mask.vox;

L2 = (mask2 ~= 0);
contours = zeros(size(L2));
contours(L2==0) = -1;
contours(L2==1) = 1;
[F,V] = isosurface(contours,0); % make a surface from the detected contours
V = V .* (ones(size(V,1),1) * mask2_vox(1:3));

load vel_struct
velocity = mrStruct.dataAy; clear mrStruct
for t = 1:size(velocity,5)
    vx = squeeze(velocity(:,:,:,1,t));
    vy = squeeze(velocity(:,:,:,2,t));
    vz = squeeze(velocity(:,:,:,3,t));
    vmagn = sqrt(vx.^2 + vy.^2 + vz.^2);
    mean_velo(t) = mean(vmagn(L2));
end

h_meanVel=figure('Name','Mean velocity');
plot(1:size(velocity,5),mean_velo,'-ko','LineWidth',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',14);
ylabel('Mean velocity (m/s)');
xlabel('Time frame #');

[I,time] = find(mean_velo==max(mean_velo));

WSS_all = Wss_point_cloud; clear Wss_point_cloud

if TimeFlag == 0  % Peak systole only
    figure(h_meanVel)
    hold on, plot(time,mean_velo(time),'-ko','LineWidth',4,...
        'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
elseif TimeFlag==1  % Averaged systolic WSS
    if size(WSS_all,2) == 1
        warndlg('WSS was previously calculated only at peak systole!');
        return;
    elseif time == 2    % second time frame is peak systole: averaging over 4 timesteps
        figure(h_meanVel)
        hold on, plot(time-1:time+2,mean_velo(time-1:time+2),'-ko','LineWidth',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    elseif time == 1 % first time frame is peak systole: averaging over 3 timesteps
        figure(h_meanVel)
        hold on, plot(time:time+2,mean_velo(time:time+2),'-ko','LineWidth',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    else
        figure(h_meanVel)
        hold on, plot(time-2:time+2,mean_velo(time-2:time+2),'-ko','LineWidth',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    end
end

for t = 1:size(WSS_all,2)
    f=figure('Name',['Wss ' num2str(t)])
    a = [4 20];
    c = [ ];
    patch('Faces',F,'Vertices',V, ...
        'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
    hold on
    [F2,V2,C2]=quiver3Dpatch(V(:,1),V(:,2),V(:,3),WSS_all{t}(:,1),WSS_all{t}(:,2),WSS_all{t}(:,3),c,a);
    patch('Faces',F2,'Vertices',V2,'CData',C2,'FaceColor','flat','EdgeColor','none','FaceAlpha',1);
    c2=colorbar;caxis([0 1.5])
    axis equal;axis off; axis ij
    view([-180 -90])
    cameratoolbar(f)
    cameratoolbar(f,'setmode','orbit')
    cameratoolbar(f,'SetCoordSys','none')
end

end