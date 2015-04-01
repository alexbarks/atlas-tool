function mimics_to_diameter(MrstructPath,MimicsSegPath,DiameterFraction,DiameterThresh,plotFlag,saveFlag)

% [mrstruct_diameter, WssStruct] = mimics_to_3Ddiameter(MrstructPath,MimicsSegPath,plotFlag,saveFlag)
%
% Reads Mimics segmentation to calculate 3D diameter and to generate diameter mrStruct
%
% The diameter is calculated for the Mimics segmentation, which is smoothed using a Laplacian filter.
% 2015, by Pim van Ooij, Northwestern University, Academic Medical Center Amsterdam
%
% Input
% 1)MrstructPath    : path for mrstruct for dataset (both mag_struct and vel_struct file, they are always in the same location)
% 2)MimicsSegPath   : path for Mimics All mask and segmentation grayvalue text file (they are always in the same location)
% 3)DiameterFraction     : fraction of diameter to be used for analysis: must be within ]0,1], default = 1
% 4)DiameterTresh        : dianeter threshold for calculation of % of data > WssTresh (incidence), default = 1
% 5)plotFlag        : plot peak systolic WSS and resulting histogram, default = 1
% 6)saveFlag        : 1 = save 3D diameter map to disk, default = 1

% Output
% 1)mrstruct_diameter_*segmentation name*      = mrstruct that contains the masked diameter values
% 2)histWssStruct_*segmentation_name*   = structure with mean & std Wss and median Wss
%  histWssStruct_*segmentation_name* =
%               Wss: array with WSS data used to generate histogram
%            median: median Wss
%              mean: mean Wss
%             stdev: standard deviation of Wss
% %
% Saved
% If saveFlag = 1:
% results are saved at the same location as the mrstruct
%
% Usage
% This code is for using a Mimics segmentation file to calculate a 3D diameter map
% Examples:
% [mrstruct_diameter, DiameterStruct] = mimics_to_Wss('All.txt','Aorta.txt','mag_struct.mat','vel_struct.mat',1,1,1,1,1)
% [mrstruct_diameter, DiameterStruct] = mimics_to_Wss() - load data using gui
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact

% % Initial error check
% mrstruct_path = '';
% mrstruct_maskedVelMag.dataAy = 0;
% histStruct.histVel = 0;

if nargin < 2 || isempty(MrstructPath)
    [MrstructPath] = uigetdir('','Select folder with mag_struct and vel_struct');
    MrstructPath = strcat(MrstructPath,'\');
end

if nargin < 3 || isempty(MimicsSegPath)
    [MimicsSegPath] = uigetdir(MrstructPath,'Select folder with all.txt and aorta.txt');
    MimicsSegPath = strcat(MimicsSegPath,'\');
end

magName = 'mag_struct.mat';
strMrstructMag = strcat(MrstructPath,magName)
allName = 'all.txt';
strAllGrayValues = strcat(MimicsSegPath,allName)
segName = 'aorta.txt';
strMimicsSeg = strcat(MimicsSegPath,segName)

% % load Allgreyvalues  file
% if ~exist(strAllGrayValues) == 2 || isempty(strAllGrayValues)
%     [all_name, all_path] = uigetfile('*.txt','Start by loading the Mimics All_greyvalues text file','Multiselect','Off');
%     strAllGrayValues = [all_path,'/',all_name];
% end
%
% % load grey value text files for Mimics segmentation
% if ~exist(strMimicsSeg) == 2 || isempty(strMimicsSeg)
%     [segName, segPath] = uigetfile('*.txt','Load your segment(s) as Mimics greyvalue text file(s)','Multiselect','Off');
%     strMimicsSeg = [segPath,'/',segName];
% end
%
% % load corresponding mag_struct
% if ~exist(strMrstructMag) == 2 || isempty(strMrstructMag)
%     [magName, magPath] = uigetfile('*.mat','Load the mag mrstruct file of your 4D flow data','Multiselect','Off');
%     strMrstructMag = [magPath,'/',magName];
% end


if nargin < 3 || isempty(DiameterFraction)
    DiameterFraction = 1;
end

if DiameterFraction <= 0;
    disp('Warning: DiameterFraction must be > 0, using DiameterFraction = 1 instead');
    DiameterFraction = 1;
end

if DiameterFraction > 1;
    disp('Warning: DiameterFraction must be <= 1, using DiameterFraction = 1 instead');
    DiameterFraction = 1;
end

if nargin < 4 || isempty(DiameterThresh)
    DiameterThresh = 3;
end

if nargin < 5 || isempty(plotFlag)
    plotFlag = 1;
end

if nargin < 6 || isempty(saveFlag)
    saveFlag = 1;
end

% Convert Mimics segmentation (grayvalues files) to mrstuct_mask
[mimicsName, mimicsPath, mrstruct_path, mrstruct_mask] = mimics_to_mrstruct(strAllGrayValues,strMimicsSeg,strMrstructMag,0,0);

% Create coordinates for the entire field of view (we want meters so: .*0.001)
[x_,y_,z_] = meshgrid((1:size(mrstruct_mask.dataAy,2)).*mrstruct_mask.vox(1).*0.001,(1:size(mrstruct_mask.dataAy,1)).*mrstruct_mask.vox(2).*0.001, ...
    (1:size(mrstruct_mask.dataAy,3)).*mrstruct_mask.vox(3).*0.001);

% But use only the ones in the segmentation
%mrstruct_mask.dataAy

% Create the faces and vertices for the lumen of the vessel from the segmentation
contours = zeros(size(mrstruct_mask.dataAy));
contours(mrstruct_mask.dataAy==0) = 1;
contours(mrstruct_mask.dataAy==1) = -1;
[F,V] = isosurface(x_,y_,z_,contours,0); % make a surface from the mask
% Smooth the lumen of the vessel (the vessel wall)
[surface_faces,surface_vertices] = SmoothLaplacian(F,V,15); %laplacian smoothing for surface (Kevin Moerman)

n = patchnormals(struct('faces',surface_faces,'vertices',surface_vertices));

Diameter = getMaximumDiameter(surface_faces,surface_vertices,n);
Diameter = Diameter*100; % make centimeters

disp(['mean Diameter = ' num2str(nanmean(Diameter))])
disp(['max Diameter = ' num2str(max(Diameter))])
disp(['min Diameter = ' num2str(min(Diameter))])

disp(['Average Diameter: ' num2str(nanmean(Diameter)) '+/-' num2str(nanstd(Diameter)) ' centimeters'])

if plotFlag
    figure('Name','Diameter');
    patch('faces',F,'vertices',V,'EdgeColor','none','FaceVertexCData',Diameter,'FaceColor','interp','faceAlpha',1);colorbar
    axis equal, axis ij, view([180 -90]), axis off
    caxis([0 6])
    axis equal;
    hold off
end

% Now create a new mask to be able to save WSS in matrix format
% We do this by dilating the original mask and then subtract the original mask
% from the dilated mask to get a mask that just has pixels at the edges (the lumen)
original_mask = double(mrstruct_mask.dataAy==1);
se = strel(ones(3,3,3));
dilated_mask = imdilate(original_mask,se);
new_mask = dilated_mask-original_mask;

% Create coordinates for these pixels;
[x_new_mask,y_new_mask,z_new_mask] = meshgrid((1:size(new_mask,2)).*mrstruct_mask.vox(1).*0.001,(1:size(new_mask,1)).*mrstruct_mask.vox(2).*0.001, ...
    (1:size(new_mask,3)).*mrstruct_mask.vox(3).*0.001);
x_coor_new = x_new_mask(new_mask==1);
y_coor_new = y_new_mask(new_mask==1);
z_coor_new = z_new_mask(new_mask==1);

% Now interpolate the Wss values from the aorta lumen coordinate cloud to the new_mask matrix so that we have Wss in matrix form
Diameter_matrix = zeros([size(mrstruct_mask.dataAy,1) size(mrstruct_mask.dataAy,2) size(mrstruct_mask.dataAy,3)]);

disp('...busy interpolating...')
interpolation_function = TriScatteredInterp([V(:,1) V(:,2) V(:,3)],Diameter,'nearest');
D_new_ = interpolation_function([x_coor_new y_coor_new z_coor_new]);
D_new_(isnan(D_new_)) = 0;
Diameter_matrix(new_mask~=0) = D_new_;

Diameter_point_cloud.Diameter = Diameter;
Diameter_point_cloud.vertices = V.*1000; %.*1000 to make it compatible with make_atlas_diameter
Diameter_point_cloud.faces = F;

% Write Diameter point cloud and other results to disk
if saveFlag
    new_name_D_point_cloud = fullfile(MrstructPath,strcat('Diameter_point_cloud_',segName(1:end-4),'.mat'));
 %   new_name_D_matrix = fullfile(MrstructPath,strcat('Diameter_struct_',segName(1:end-4),'.mat'));
 %   new_name_mask = fullfile(MrstructPath,strcat('mask_struct_',segName(1:end-4),'.mat'));
    save(new_name_D_point_cloud,'Diameter_point_cloud');
 %   save(new_name_D_matrix,'Diameter_matrix');
 %   save(new_name_mask,'mrstruct_mask');
    if isempty(MrstructPath)
        disp('Could not write Diameter_point_cloud');
    end
    disp('saved')
end

% Analyze masked Diameter data and generate histogram
% get size of data and indices along time direction to extract user selected % of highest Wss
sz      = size(Diameter_matrix);
if size(sz,2)<4   % Emilie: only one time (peak systole)
    sz(4) = 1;
end
Diameter_ind = round( (1-DiameterFraction)*sz(4) ) + 1;
if Diameter_ind > sz(4)
    Diameter_ind = sz(4);
end

fDiameter=[]; % variable to collect data for histogram analysis

for zi=1:sz(1)
    %disp('.');
    %drawnow;
    for xi=1:sz(2)
        for yi=1:sz(3)
            if(new_mask(zi,xi,yi))
                tmp=squeeze(Diameter_matrix(zi,xi,yi,:));
                stmp=sort(tmp);
                fDiameter=[fDiameter; stmp(Diameter_ind:end)];
            end
        end
    end
end

histStruct.histDiameter = double(fDiameter);
% histogram analysis
maxDiameter  = max(histStruct.histDiameter);
histStruct.median    = median(histStruct.histDiameter);
histStruct.mean      = mean(histStruct.histDiameter);
histStruct.stdev     = std(histStruct.histDiameter);
f = find(histStruct.histDiameter >= DiameterThresh);
histStruct.incidence = length(f)./length(histStruct.histDiameter)*100;

% plot results
if plotFlag
    
    figure;
    clf;
    
    % plot velcoity MIP
    subplot(1,2,1);
    imagesc(squeeze(max(max(Diameter_matrix,[],3),[],4)));
    colorbar;caxis([0 6])
    axis off
    title('4D Diameter MIP','FontSize',18)
    
    % plot normalized velocity histogram
    subplot(1,2,2);
    binLocs = 0:maxDiameter/200:maxDiameter;
    numDiameterperBin         = hist(histStruct.histDiameter,binLocs);
    ylim=1.25*max(numDiameterperBin./length(histStruct.histDiameter));
    
    bar(binLocs,numDiameterperBin./length(histStruct.histDiameter))
    axis([0 maxDiameter 0 ylim])
    line([histStruct.mean histStruct.mean],[0 0.8*ylim],'Color',[0 0 0])
    text(double(histStruct.mean),0.82*ylim,['mean Diameter: ',num2str(histStruct.mean,3)],'FontSize',14)
    line([histStruct.median histStruct.median],[0 0.9*ylim],'Color',[0 0 0])
    text(double(histStruct.median),0.92*ylim,['median Diameter: ',num2str(histStruct.median,3)],'FontSize',14)
    line([DiameterThresh DiameterThresh],[0 0.7*ylim],'Color',[0 0 0])
    text(double(DiameterThresh),0.72*ylim,['Percent above ',num2str(DiameterThresh),'cm: ',num2str(histStruct.incidence,3)],'FontSize',14)
    xlabel('Diameter [cm]','FontSize',18)
    set(gca,'FontSize',16)
    title('Diameter Histogram')
    
end

%----------------------------------------------------------------------------------------------------------

function D = getMaximumDiameter(F,V,N)
% D = getMaximumDiameter(F,V,N)
% diameters D [X 1]
% Faces     F [Y 3] (connection of edge points)
% Vertices  V [X 3] (edge points of triangles)
% Normals   N [X 3] (inward normals of points)
% uses intersectLineTriangle3d by David Legland, 2011
% uses patchnormals by D. Croon
% (Both obtained through MathWorks FileExchange)
%
% Code was edited by Wouter Potters, 2012, Pim van Ooij, 2015

%t1 = tic;
D = nan(size(V,1),1);

% operation takes +/-420 seconds on 1 CPU on Mac OS X laptop with ...
% points

% mlps = matlabpool('size'); % matlab pool size // number of cpu's available
% if mlps == 0; mlps = 1; end
% total_time = num2str(ceil((0.027*numel(V(:,1))) / (mlps*60))+1);
% currenttime = sprintf('%.0f%.0f%.0f:%.0f,%.0f',((fix(clock) .* [0 0 1 1 1 1])'));
% if mlps > 1, cputxt = '''s'; else cputxt =''; end
% disp(['busy calculating ray triangle intersections.... ' total_time ' minutes left on ' num2str(mlps) ' cpu' cputxt ' starting from: ' currenttime])

tic
disp('...busy calculating...')
parfor i_v = 1:length(V)
    [~, posOfIntersection, ~] = intersectLineMesh3d([V(i_v,:) N(i_v,:)], V, F);
    
    posOfIntersection(posOfIntersection < 10^-10) = [];
    
    if ~isempty(posOfIntersection)
        D(i_v) = abs(posOfIntersection(find( abs(posOfIntersection) == min(abs(posOfIntersection)), 1, 'first' )))
    end
end
toc
%toc(t1)
%disp(['finished calculating ' length(V)^2 ' ray-triangle intersections! At time: ' currenttime])
disp('Finished!')