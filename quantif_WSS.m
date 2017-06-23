function quantif_WSS(TimeFlag)

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
if TimeFlag==2
    choice = questdlg('Do you want to extract WSS values calculated...', ...
        'WSS regional quantification', ...
        'at peak systole?','while averaging up to 5 systolic timesteps?','at peak systole?');
    % Handle response
    switch choice
        case 'at peak systole?'
            TimeFlag=0;
        case 'while averaging up to 5 systolic timesteps?'
            TimeFlag=1;
    end
end
if TimeFlag==0
    % Peak systolic WSS
    figure(h_meanVel)
    hold on, plot(time,mean_velo(time),'-ko','LineWidth',4,...
        'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    if size(WSS_all,2) > 5
        WSS = WSS_all{time};
    elseif size(WSS_all,2) == 5
        WSS = WSS_all{3};
    elseif size(WSS_all,2) == 4
        WSS = WSS_all{2};
    elseif size(WSS_all,2) == 3
        WSS = WSS_all{1};
    elseif size(WSS_all,2) == 1    % Emilie: calculated at only one time (peak systole)
        WSS = WSS_all{1};
    end
    wss_m = sqrt(WSS(:,1).^2 + WSS(:,2).^2 + WSS(:,3).^2);
elseif TimeFlag==1  % Averaged systolic WSS
    if size(WSS_all,2) == 1
        warndlg('WSS was previously calculated only at peak systole!');
        return;
        % Velocity averaged over x systolic time frames
    elseif time == 2    % second time frame is peak systole: averaging over 4 timesteps
        disp('AVERAGE OVER 4 TIME FRAMES!')
        data2.x_value_wss_t1 = WSS_all{time-1}(:,1);data2.y_value_wss_t1 = WSS_all{time-1}(:,2);data2.z_value_wss_t1 = WSS_all{time-1}(:,3);
        data2.x_value_wss_t2 = WSS_all{time}(:,1);  data2.y_value_wss_t2 = WSS_all{time}(:,2);  data2.z_value_wss_t2 = WSS_all{time}(:,3);
        data2.x_value_wss_t3 = WSS_all{time+1}(:,1);data2.y_value_wss_t3 = WSS_all{time+1}(:,2);data2.z_value_wss_t3 = WSS_all{time+1}(:,3);
        data2.x_value_wss_t4 = WSS_all{time+2}(:,1);data2.y_value_wss_t4 = WSS_all{time+2}(:,2);data2.z_value_wss_t4 = WSS_all{time+2}(:,3);
        data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4)./4;
        data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4)./4;
        data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4)./4;
        figure(h_meanVel)
        hold on, plot(time-1:time+2,mean_velo(time-1:time+2),'-ko','LineWidth',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    elseif time == 1 % first time frame is peak systole: averaging over 3 timesteps
        disp('AVERAGE OVER 3 TIME FRAMES!')
        data2.x_value_wss_t1 = WSS_all{time}(:,1);  data2.y_value_wss_t1 = WSS_all{time}(:,2);  data2.z_value_wss_t1 = WSS_all{time}(:,3);
        data2.x_value_wss_t2 = WSS_all{time+1}(:,1);data2.y_value_wss_t2 = WSS_all{time+1}(:,2);data2.z_value_wss_t2 = WSS_all{time+1}(:,3);
        data2.x_value_wss_t3 = WSS_all{time+2}(:,1);data2.y_value_wss_t3 = WSS_all{time+2}(:,2);data2.z_value_wss_t3 = WSS_all{time+2}(:,3);
        data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3)./3;
        data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3)./3;
        data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3)./3;
        figure(h_meanVel)
        hold on, plot(time:time+2,mean_velo(time:time+2),'-ko','LineWidth',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
    else % timestep > 2 is peak systole: averaging over 5 timesteps
        if size(WSS_all,2) > 5
            data2.x_value_wss_t1 = WSS_all{time-2}(:,1);data2.y_value_wss_t1 = WSS_all{time-2}(:,2);data2.z_value_wss_t1 = WSS_all{time-2}(:,3);
            data2.x_value_wss_t2 = WSS_all{time-1}(:,1);data2.y_value_wss_t2 = WSS_all{time-1}(:,2);data2.z_value_wss_t2 = WSS_all{time-1}(:,3);
            data2.x_value_wss_t3 = WSS_all{time}(:,1);  data2.y_value_wss_t3 = WSS_all{time}(:,2);  data2.z_value_wss_t3 = WSS_all{time}(:,3);
            data2.x_value_wss_t4 = WSS_all{time+1}(:,1);data2.y_value_wss_t4 = WSS_all{time+1}(:,2);data2.z_value_wss_t4 = WSS_all{time+1}(:,3);
            data2.x_value_wss_t5 = WSS_all{time+2}(:,1);data2.y_value_wss_t5 = WSS_all{time+2}(:,2);data2.z_value_wss_t5 = WSS_all{time+2}(:,3);
            figure(h_meanVel)
            hold on, plot(time-2:time+2,mean_velo(time-2:time+2),'-ko','LineWidth',4,...
                'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
        elseif size(WSS_all,2) == 5
            data2.x_value_wss_t1 = WSS_all{1}(:,1);data2.y_value_wss_t1 = WSS_all{1}(:,2);data2.z_value_wss_t1 = WSS_all{1}(:,3);
            data2.x_value_wss_t2 = WSS_all{2}(:,1);data2.y_value_wss_t2 = WSS_all{2}(:,2);data2.z_value_wss_t2 = WSS_all{2}(:,3);
            data2.x_value_wss_t3 = WSS_all{3}(:,1);  data2.y_value_wss_t3 = WSS_all{3}(:,2);  data2.z_value_wss_t3 = WSS_all{3}(:,3);
            data2.x_value_wss_t4 = WSS_all{4}(:,1);data2.y_value_wss_t4 = WSS_all{4}(:,2);data2.z_value_wss_t4 = WSS_all{4}(:,3);
            data2.x_value_wss_t5 = WSS_all{5}(:,1);data2.y_value_wss_t5 = WSS_all{5}(:,2);data2.z_value_wss_t5 = WSS_all{5}(:,3);
            figure(h_meanVel)
            hold on, plot(1:5,mean_velo(1:5),'-ko','LineWidth',4,...
            'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',14);
        end
        data2.x_value_wss = (data2.x_value_wss_t1 + data2.x_value_wss_t2 + data2.x_value_wss_t3 + data2.x_value_wss_t4 + data2.x_value_wss_t5)./5;
        data2.y_value_wss = (data2.y_value_wss_t1 + data2.y_value_wss_t2 + data2.y_value_wss_t3 + data2.y_value_wss_t4 + data2.y_value_wss_t5)./5;
        data2.z_value_wss = (data2.z_value_wss_t1 + data2.z_value_wss_t2 + data2.z_value_wss_t3 + data2.z_value_wss_t4 + data2.z_value_wss_t5)./5;
    end
    wss_m = sqrt(data2.x_value_wss.^2 + data2.y_value_wss.^2 + data2.z_value_wss.^2);
end

mask2_vox = mrstruct_mask.vox;
x = V(:,1)/mask2_vox(1);
y = V(:,2)/mask2_vox(2);
z = V(:,3)/mask2_vox(3);
figure, patch('Faces',F,'Vertices',[x y z], ...
    'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',1);
axis equal;axis off; axis ij
view([-180 -90])

gray_colormap = colormap(gray);
color3(1,:) = [0 0 1];
color3(2,:) = [1 0 0];
color3(3,:) = [0.5 0.5 0.5];
color3(4:64,:) = gray_colormap(4:64,:);
colormap(color3);
caxis([0 64]);
axis equal; axis ij; axis off;
aspectRatio = 1./mask2_vox;
% set(gca,'dataaspectRatio',aspectRatio(1:3))
view([-180 -90]);
load mag_struct
magnitude = flipdim(double(mrStruct.dataAy(:,:,:,3)),3);
magnitude(magnitude == 0) = 3;
magnitude(magnitude == 1) = 3;
magnitude(magnitude == 2) = 3;
hold on
s4 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]) .* size(magnitude,3)/2,magnitude(:,:,size(magnitude,3)/2),'EdgeColor','none');
% title({'Draw 1)proximal AAo inner, 2)proximal AAo outer, 3)distal AAo inner, 4)distal AAo outer,';'5)arch inner, 6)arch outer, 7)proximal DAo inner, 8)proximal DAo outer,';'9)distal DAo inner and 10)distal DAo outer regions';'then double-click and press space'})
title({'Please keep in mind that regions of interest will be numbered in the same order you draw them';'once you''re done drawing an ROI, please double-click to validate and press space to move on to the next one'})

uicontrol('Style','text',...
    'Position',[10 200 120 70],...
    'String','Please choose using the slider the magnitude slice')
uicontrol('Style','text',...
    'Position',[10 75 120 20],...
    'String','Slice slider')
sl1 = uicontrol('Style', 'slider',...
    'Min',1,'Max',size(magnitude,3),'Value',size(magnitude,3)/2,...
    'Position', [10 50 120 20],...
    'SliderStep',[1/(size(magnitude,3)-1) 10/(size(magnitude,3)-1)],...
    'Callback', {@move_slice3,gca});

mkdir(strcat(MrstructPath, '..'),'regional_masks')

choice = questdlg('Do you want to draw new ROIs or load a previous one?', ...
    'ROIs creation', ...
    'Draw new','Load existing','Draw new');
% Handle response
switch choice
    case 'Draw new'
        
        nbROIs = inputdlg('How many ROIs do you need?','Draw ROIs',1,{'10'});
        nbROIs = (round(str2num(nbROIs{1})));
        
        for i = 1:nbROIs
            %Polygon and mask
            polyAAo = impoly;
            wait(polyAAo);
            region = getPosition(polyAAo);
            disp('saving, pausing')
            save(strcat([MrstructPath '..' '\regional_masks\mask' num2str(i)]),'region');
            mask_wss = inpolygon(x,y, region(:,1), region(:,2));
            % compute WSS in the region
            wss_mask{i} = wss_m(mask_wss);
            clear mask_wss region
            pause
        end
        
        h1 = msgbox('ROIs drawn, WSS quantification in progress...');
                
        if TimeFlag==0
            mat_file = strcat(MrstructPath,'..','\regional_masks\wss_values_syst');
        elseif TimeFlag==1
            mat_file = strcat(MrstructPath,'..','\regional_masks\wss_values_avg');
        end
        save(mat_file,'wss_mask');
        
        % compute quantitative indices
        indices{1,1} = '';
        indices{2,1} = 'mean';
        indices{3,1} = 'median';
        indices{4,1} = 'max';
        indices{5,1} = 'min';
        indices{6,1} = 'std';
        indices{7,1} = 'max5percent';
        indices{8,1} = 'max2percent';
        indices{9,1} = 'min5percent';
        indices{10,1} = 'min2percent';
        
        for i=1:nbROIs
            
            mask_wss=wss_mask{i};
            indices{1,i+1} = strcat(['region' num2str(i)]);
            indices{2,i+1} = mean(mask_wss(~isnan(mask_wss)));
            indices{3,i+1} = median(mask_wss(~isnan(mask_wss)));
            indices{4,i+1} = max(mask_wss);
            indices{5,i+1} = min(mask_wss);
            indices{6,i+1} = std(mask_wss(~isnan(mask_wss)));
            WSS_sorted = sort(mask_wss(~isnan(mask_wss)));
            indices{7,i+1} = mean(WSS_sorted(end-5/100*ceil(length(WSS_sorted)):end));
            indices{8,i+1} = mean(WSS_sorted(end-2/100*ceil(length(WSS_sorted)):end));
            indices{9,i+1} = mean(WSS_sorted(1:5/100*ceil(length(WSS_sorted))));
            indices{10,i+1} = mean(WSS_sorted(1:2/100*ceil(length(WSS_sorted))));
            
            clear mask_wss WSS_sorted
        end
        
        currDir=pwd;
        cd(strcat(MrstructPath,'..','\regional_masks'))
        % save in an Excel sheet
        if TimeFlag==0
            xls_file = 'wss_indices_syst.xls';
        elseif TimeFlag==1
            xls_file = 'wss_indices_avg.xls';
        end
        xlswrite(xls_file,indices);
        cd(currDir)
        close(h1)
        h1 = msgbox('WSS quantification done, results are saved in the regional_masks folder');
        
    case 'Load existing'
        
        choice = questdlg('Do you want to...', ...
            'Load existing ROIs', ...
            'Recompute WSS in all existing ROIs','Modify one ROI','Modify all ROIs','Modify one ROI');
        % Handle response
        switch choice
            case 'Modify one ROI'
                [FileName,MrstructPath,FilterIndex] = uigetfile(currDir,'Select the ROI you want to load');
                currentDir=pwd;
                cd(MrstructPath);
                masks=ls('mask*');
                for i=1:size(masks,1)
                    load(strcat(MrstructPath,masks(i,:)));
                    hold on, plot([region(:,1);,region(1,1)],[region(:,2);region(1,2)])
                    clear region
                end
                cd(currentDir);
                
                load(strcat(MrstructPath,FileName));
                
                polyAAo = impoly(gca,region);
                wait(polyAAo);
                region = getPosition(polyAAo);
                disp('saving, pausing')
                save(strcat(MrstructPath,FileName),'region');
                pause
                
                h1 = msgbox('ROI changed, updated WSS quantification in progress...');
                
                mask_wss = inpolygon(x,y, region(:,1), region(:,2));
                wss_mask = wss_m(mask_wss);
                if TimeFlag==0
                    mat_file = strcat([MrstructPath 'wss_values_syst_new' FileName]);
                elseif TimeFlag==1
                    mat_file = strcat([MrstructPath 'wss_values_avg_new' FileName]);
                end
                save(mat_file,'wss_mask');
                
                new_indices{1,1} = mean(wss_mask(~isnan(wss_mask)));
                new_indices{2,1} = median(wss_mask(~isnan(wss_mask)));
                new_indices{3,1} = max(wss_mask);
                new_indices{4,1} = min(wss_mask);
                new_indices{5,1} = std(wss_mask(~isnan(wss_mask)));
                WSS_sorted = sort(wss_mask(~isnan(wss_mask)));
                new_indices{6,1} = mean(WSS_sorted(end-5/100*ceil(length(WSS_sorted)):end));
                new_indices{7,1} = mean(WSS_sorted(end-2/100*ceil(length(WSS_sorted)):end));
                new_indices{8,1} = mean(WSS_sorted(1:5/100*ceil(length(WSS_sorted))));
                new_indices{9,1} = mean(WSS_sorted(1:2/100*ceil(length(WSS_sorted))));
                
                currDir=pwd;
                cd(MrstructPath)
                % save in an Excel sheet
                %         col = char(str2num(FileName(end-4))+'A');
                indMask = strfind(FileName, 'mask');
                indMat = strfind(FileName, '.mat');
                col = char(str2num(FileName(indMask+4:indMat-1))+'A');
                xlRange = strcat([col '2:' col '10']);
                if TimeFlag==0
                    xls_file = 'wss_indices_syst.xls';
                elseif TimeFlag==1
                    xls_file = 'wss_indices_avg.xls';
                end
                xlswrite(xls_file,new_indices,xlRange);
                cd(currDir)
                close(h1);
                h1 = msgbox('WSS quantification in the new ROI was updated and saved in the regional_masks folder');
                
            case 'Recompute WSS in all existing ROIs'
                
                currentDir=pwd;
                cd(strcat([MrstructPath '..' '\regional_masks']));
                masks=ls('mask*');
                for i=1:size(masks,1)
                    load(strcat([MrstructPath '..' '\regional_masks\mask' num2str(i)]));
                    hold on, plot([region(:,1);,region(1,1)],[region(:,2);region(1,2)])
                    mask_wss = inpolygon(x,y, region(:,1), region(:,2));
                    wss_mask{i} = wss_m(mask_wss);
                    clear region mask_wss
                end
                cd(currentDir);
                
                if TimeFlag==0
                    mat_file = strcat(MrstructPath,'..','\regional_masks\wss_values_syst');
                elseif TimeFlag==1
                    mat_file = strcat(MrstructPath,'..','\regional_masks\wss_values_avg');
                end
                save(mat_file,'wss_mask');
                
                % compute quantitative indices
                indices{1,1} = '';
                indices{2,1} = 'mean';
                indices{3,1} = 'median';
                indices{4,1} = 'max';
                indices{5,1} = 'min';
                indices{6,1} = 'std';
                indices{7,1} = 'max5percent';
                indices{8,1} = 'max2percent';
                indices{9,1} = 'min5percent';
                indices{10,1} = 'min2percent';
                
                for i=1:size(masks,1)
                    
                    mask_wss=wss_mask{i};
                    indices{1,i+1} = strcat(['region' num2str(i)]);
                    indices{2,i+1} = mean(mask_wss(~isnan(mask_wss)));
                    indices{3,i+1} = median(mask_wss(~isnan(mask_wss)));
                    indices{4,i+1} = max(mask_wss);
                    indices{5,i+1} = min(mask_wss);
                    indices{6,i+1} = std(mask_wss(~isnan(mask_wss)));
                    WSS_sorted = sort(mask_wss(~isnan(mask_wss)));
                    indices{7,i+1} = mean(WSS_sorted(end-5/100*ceil(length(WSS_sorted)):end));
                    indices{8,i+1} = mean(WSS_sorted(end-2/100*ceil(length(WSS_sorted)):end));
                    indices{9,i+1} = mean(WSS_sorted(1:5/100*ceil(length(WSS_sorted))));
                    indices{10,i+1} = mean(WSS_sorted(1:2/100*ceil(length(WSS_sorted))));
                    
                    clear mask_wss WSS_sorted
                end
                
                currDir=pwd;
                cd(strcat(MrstructPath,'..','\regional_masks'))
                % save in an Excel sheet
                if TimeFlag==0
                    xls_file = 'wss_indices_syst.xls';
                elseif TimeFlag==1
                    xls_file = 'wss_indices_avg.xls';
                end
                xlswrite(xls_file,indices);
                cd(currDir)
                h1 = msgbox('WSS indices were recalculated and saved in the regional_masks folder');
                
            case 'Modify all ROIs'
                
                currentDir=pwd;
                cd(strcat([MrstructPath '..' '\regional_masks']));
                masks=ls('mask*');
                for i=1:size(masks,1)
                    load(strcat([MrstructPath '..' '\regional_masks\mask' num2str(i)]));
                    hold on, plot([region(:,1);,region(1,1)],[region(:,2);region(1,2)])
                    clear region
                end
                for i=1:size(masks,1)
                    load(strcat([MrstructPath '..' '\regional_masks\mask' num2str(i)]));
                    %Polygon and mask
                    polyAAo = impoly(gca,region);
                    wait(polyAAo);
                    region = getPosition(polyAAo);
                    disp('saving, pausing')
                    save(strcat([MrstructPath '..' '\regional_masks\mask' num2str(i)]),'region');
                    mask_wss = inpolygon(x,y, region(:,1), region(:,2));
                    % compute WSS in the region
                    wss_mask{i} = wss_m(mask_wss);
                    clear mask_wss region
                    pause
                end
                cd(currentDir);
                h1 = msgbox('ROI changed, updated WSS quantification in progress...');
                
                if TimeFlag==0
                    mat_file = strcat(MrstructPath,'..','\regional_masks\wss_values_syst');
                elseif TimeFlag==1
                    mat_file = strcat(MrstructPath,'..','\regional_masks\wss_values_avg');
                end
                save(mat_file,'wss_mask');
                
                % compute quantitative indices
                indices{1,1} = '';
                indices{2,1} = 'mean';
                indices{3,1} = 'median';
                indices{4,1} = 'max';
                indices{5,1} = 'min';
                indices{6,1} = 'std';
                indices{7,1} = 'max5percent';
                indices{8,1} = 'max2percent';
                indices{9,1} = 'min5percent';
                indices{10,1} = 'min2percent';
                
                for i=1:size(masks,1)
                    
                    mask_wss=wss_mask{i};
                    indices{1,i+1} = strcat(['region' num2str(i)]);
                    indices{2,i+1} = mean(mask_wss(~isnan(mask_wss)));
                    indices{3,i+1} = median(mask_wss(~isnan(mask_wss)));
                    indices{4,i+1} = max(mask_wss);
                    indices{5,i+1} = min(mask_wss);
                    indices{6,i+1} = std(mask_wss(~isnan(mask_wss)));
                    WSS_sorted = sort(mask_wss(~isnan(mask_wss)));
                    indices{7,i+1} = mean(WSS_sorted(end-5/100*ceil(length(WSS_sorted)):end));
                    indices{8,i+1} = mean(WSS_sorted(end-2/100*ceil(length(WSS_sorted)):end));
                    indices{9,i+1} = mean(WSS_sorted(1:5/100*ceil(length(WSS_sorted))));
                    indices{10,i+1} = mean(WSS_sorted(1:2/100*ceil(length(WSS_sorted))));
                    
                    clear mask_wss WSS_sorted
                end
                
                currDir=pwd;
                cd(strcat(MrstructPath,'..','\regional_masks'))
                % save in an Excel sheet
                if TimeFlag==0
                    xls_file = 'wss_indices_syst.xls';
                elseif TimeFlag==1
                    xls_file = 'wss_indices_avg.xls';
                end
                xlswrite(xls_file,indices);
                cd(currDir)
                close(h1)
                h1 = msgbox('WSS indices were calculated and saved in the regional_masks folder');

        end
end

    function move_slice3(hObj,event,ax) % Emilie: for manual interaction to chose magnitude slice on which RPA can be visualized
        sliceobj = findobj(s4);
        delete(sliceobj)
        slice = round(get(hObj,'Value'));
        s4 = surf(1:size(magnitude,2),1:size(magnitude,1),ones([size(magnitude,1) size(magnitude,2)]).*(size(magnitude,3)-(slice-1)),magnitude(:,:,slice),'EdgeColor','none');
    end

end