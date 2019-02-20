isSim = 0;

[dataFile, dataPath] = uigetfile({'*.mat';'*.*'},'MultiSelect', 'on');


% plotPhase = 1;
% if plotPhase == 1
%     [dataFile2, dataPath2] = uigetfile({'*.mat';'*.*'},'Open oufti pre process info');
% load ([dataPath2,dataFile2]);
% %     [whiteLightFile, whiteLightPath] = uigetfile({'*.tiff';'*.*'},'Open phase image');
% %     whiteLightInfo = imfinfo([whiteLightPath whiteLightFile]);
% %     whiteLight = zeros(whiteLightInfo(1).Height, whiteLightInfo(1).Width);
% %     whiteLight = whiteLight + double(imread([whiteLightPath whiteLightFile],'Info', whiteLightInfo));
% %     
% %     %set pixel values
% %     wlnmPerPixel = 110;
% %     nmPerPixel = 108;
%     
% %     %select channel
% %     dlg_title = 'Input info';
% %     prompt = {'Channel (r/g)','Filename (ex. S1R1)'};
% %     def = {'g','S1R1'};
% %     num_lines = 1;
% %     inputInfo = inputdlg(prompt,dlg_title,num_lines,def);
% %     channel = inputInfo{1};
% %     filename = inputInfo{2}; 
% %     
% %     %resize and rotate phase image
% %     whiteLight = imresize(whiteLight,wlnmPerPixel/nmPerPixel);
% %     WLHeight = size(whiteLight,1);
% %     WLWidth = size(whiteLight,2);
% %     whiteLight = rot90(whiteLight);
% %     wlCenterX = (length(whiteLight(1,:))*wlnmPerPixel)/2;
% %     wlCenterY = (length(whiteLight(:,1))*wlnmPerPixel)/2;
% %     if channel == 'r'
% %         wlShiftX = 133 * nmPerPixel;
% %         wlShiftY = 393 * nmPerPixel;
% %     elseif channel == 'g'
% %         wlShiftX = 150 * nmPerPixel;
% %         wlShiftY = 179 * nmPerPixel;
% %         whiteLight = flipud(whiteLight);
% %     end
%     
% whiteLightNew = zeros(WLWidth,WLHeight);
% whiteLightNew(wlYLim1:wlYLim2,wlXLim1:wlXLim2) = double(whiteLightCrop);
% cropShiftX = wlXLim1*nmPerPixel;
% cropShiftY = wlYLim1*nmPerPixel;
% 
% wlCenterX = (length(whiteLight(1,:))*wlnmPerPixel)/2;
% wlCenterY = (length(whiteLight(:,1))*wlnmPerPixel)/2;
% if channel == 'r'
%         wlShiftX = 133 * nmPerPixel; 
%         wlShiftY = 393 * nmPerPixel; 
% elseif channel == 'g'
%         wlShiftX = 150 * nmPerPixel;
%         wlShiftY = 179 * nmPerPixel;
% end
% 
% %plot white light, localizations, and cell outlines
% figure;
% [xWL, yWL] = meshgrid(0:WLHeight,0:WLWidth);
% xRange = xWL(1,:)*nmPerPixel-wlShiftX;
% yRange = yWL(:,1)*nmPerPixel-wlShiftY;
% imagesc(xRange,yRange,whiteLightNew);
% colormap gray;
% hold on;   
% end


%user inputs
dlg_title = 'Input info';
prompt = {'Exposure Time','Match Distance','Max Frames to Skip'};
def = {'25','2500','1'};
% def = {'25','3000','1'};
num_lines = 1;
inputInfo = inputdlg(prompt,dlg_title,num_lines,def);
time_per_frame = str2num(inputInfo{1});
matchdist = str2num(inputInfo{2});
frameSkip = str2num(inputInfo{3});

if iscell(dataFile)
% fileNum = 1:length(dataFile);
fileNum = length(dataFile);

else
    fileNum = 1;
end

for j = 1:fileNum

    if iscell(dataFile)
        load ([dataPath,dataFile{1,j}]);
    else
        load ([dataPath,dataFile]);
    end
    
    
% figure;    
tracks3 = [];   
tic
numCells = length(fiberData);
totalDiff = [];
totalDiff2 = [];
totTrackLength = [];
frameRange = fiberData(1).FrameRange;
corrected_xyz = [];
for c = 1:numCells
corrected_xyz(:,1) = fiberData(c).xLoc;
corrected_xyz(:,2) = fiberData(c).yLoc;
corrected_xyz(:,3) = fiberData(c).zLoc;
% corrected_xyz(:,1) = fiberData(c).xMean;
% corrected_xyz(:,2) = fiberData(c).yMean;
% corrected_xyz(:,3) = fiberData(c).zMean;
frameNum = fiberData(c).Frame;
points = [];
c
for a = 1:frameRange(1,2)
points{a,1} = corrected_xyz(frameNum(:,1)==a,:);
end


% tracks = Tracking(points,time_per_frame,matchdist,frameSkip,frameRange);

% % for tracking of simulated points
% time = (0:length(points)-1)*time_per_frame;
% temp = ~cellfun(@isempty,points);
% tempTime = time(temp);
% tempTracks1 = points(temp);
% tempTracks = [];
% for j = 1:length(tempTime)
% tempTracks = vertcat(tempTracks,tempTracks1{j});
% end
% tracks = cell(1);
% tracks{1} = zeros(length(tempTime),4);
% tracks{1}(:,1) = tempTime;
% tracks{1}(:,2:4) = tempTracks;
% % end tracking of simulated points

if isSim == 1
    for j = 1:length(corrected_xyz)/6 %number of tracks (tracklength of 6)
    tracks{j} = zeros(6,4);
    tracks{j}(:,1) = [25:25:150];
    tracks{j}(:,2:4) = corrected_xyz((j-1)*6+1:j*6,:);
    end
%     tracks{1} = zeros(6,4);
%     tracks{1}(:,1) = [25:25:150];
%     tracks{1}(:,2:4) = corrected_xyz;
else
[tracks, tracks_no_t_c, tracks_time_coord] = simpletracker(points,time_per_frame, 'Method','NearestNeighbor','MaxLinkingDistance',matchdist,'MaxGapClosing',1,'Debug',false);
% [tracks, tracks_no_t_c, tracks_time_coord] = simpletracker_4steps(points,time_per_frame, 'Method','NearestNeighbor','MaxLinkingDistance',matchdist,'MaxGapClosing',1,'Debug',false);
end


% simpletracker(points, time_per_frame,'Method','Hungarian','MaxLinkingDistance',2000,'MaxGapClosing',1,'Debug',false);
% load('trajectories.mat');

% %start first method to caclculate MSD
% ma = msdanalyzer(3,'nm','ms');
% ma = ma.addAll(tracks);
% ma = ma.computeMSD;

% ma.fitMeanMSD;

% r2Thresh = 0.70;
% % allMSD = ma.fitMSD;
% % allMSD = allMSD.lfit.a(allMSD.lfit.r2fit >= r2Thresh,1);
% % % allMSD = allMSD.lfit.a;
% % allDiff = allMSD/10^3;
% % 
% % totalDiff = vertcat(totalDiff,allDiff);

% figure;
% % edges = 0:0.1:ceil(max(allMSD));
% edges = 0:0.1:10;
% histogram(allMSD,edges);
% title(['Diffusion Coefficients: ' num2str(length(allMSD)) ' Tracks']);
% xlabel('Diffusion Coefficient (µm^2/s)');
% ylabel('Number of Occurances');

% diffThresh = 0.1;
% figure; 
% ha = gca;
% title({[num2str(length(ma.tracks)) ' Tracks Total'],[num2str(sum(allMSD > diffThresh)) ' Tracks Mobile - red'],...
%     [num2str(sum(allMSD <= diffThresh)) ' Tracks Stationary - blue']});
% 
% n_tracks = numel(allMSD);
% colors = jet(n_tracks);
% 
% hold(ha, 'on');
% for i = 1 : n_tracks
%     
%     track = ma.tracks{i};
%     
%     x = track(:,2);
%     y = track(:,3);
%     z = track(:,4);
%         
%     if allMSD(i,1) > diffThresh
%     plot3(ha, x, y, z, 'Color', 'r');
%     else
%     plot3(ha, x, y, z, 'Color', 'b');
%     end
%     axis vis3d equal;
% end
%end second method to calculate MSD

%    %plot cell outlines
%    hold on; 
%    scatter(fiberData(c).meshPtsX, fiberData(c).meshPtsY,'.','g','LineWidth',2);
%    hold on;
%    plot3(mean(fiberData(c).meshPtsX), mean(fiberData(c).meshPtsY), -118,'.k','MarkerSize',15);
%    %end plot cell outlines

%start second (correct) method for calculating MSD
% tracks = [];
% tracks = tracks;
% tracks2 = ma.tracks;
% single_frame_sq_displacements = [];
diff2 = [];
trackLength = zeros(length(tracks),1);
% figure;
for e = 1:length(tracks)
    single_frame_sq_displacements = [];
%     for f = 1:length(tracks2{e}) - 1
    for f = 1:length(tracks{e}(:,1)) - 1

        single_frame_sq_displacements(f,1) = sum((tracks{e}(f,2:4) - tracks{e}(f+1,2:4)).^2,2);
    end
    diff2(e,1) = mean(single_frame_sq_displacements)*10^-6/((2*3)*(time_per_frame*10^-3));
    
    
%   %plot tracks
%     x = tracks2{e}(:,2);
%     y = tracks2{e}(:,3);
%     z = tracks2{e}(:,4);
%     hold on
%     plot3( x, y, z, 'Color', 'r');
%     
% %     hold on;
% %     if diff2(e,1) > 0.25
% %     plot3( x, y, z, 'Color', 'r');
% %     else
% %     plot3( x, y, z, 'Color', 'b');
% %     end
% %     axis vis3d equal;
%     %end plot tracks
    
trackLength(e) = length(tracks{e});

end
totalDiff2 = vertcat(totalDiff2,diff2);
totTrackLength = vertcat(totTrackLength,trackLength);
%end second method


tracks3{c} = tracks;
clear corrected_xyz points ma allMSD allDiff x y z track colors
% c
end

totalDiff = totalDiff/(2*3);
totalMeanDiff = mean(totalDiff);

totalMeanDiff2 = mean(totalDiff2)


% figure;
% hold on;
% title(['Diffusion Coefficients: ' num2str(length(totalDiff2)) ' Tracks, ' num2str(numCells) ' Cells ']);
% xlabel('Diffusion Coefficient (µm^2/s)');
% ylabel('Probability Density');
% 
% edges = 0:0.05:10;
% binSpacing = edges(2)-edges(1);
% 
% histogram_counts = histcounts(totalDiff2,(edges+binSpacing/2));
% pdf_distribution = histogram_counts./(sum(histogram_counts)*binSpacing);
% scale = sum(histogram_counts)*binSpacing;
% bar1 = bar(edges(2:end),pdf_distribution);
% 
% k = 1.56;
% % initGuess = [3 4];                    %1 species
% initGuess = [0.1 0.1 1.5];              %2 species
% % initGuess = [0.1 0.6 0.1 2 4 4 4 4];  %3 species
% number_of_species =2;
% x_values = edges(2:end);
% hold on;
% hist_fitting(pdf_distribution,x_values,k, initGuess, number_of_species);
% totalMeanDiff2 = mean(totalDiff2);
% 
toc
%save necessary info
diffusionCoefficients = totalDiff2;

%  if iscell(dataFile)
%      fileName = strsplit(dataFile{1,j}, '.');
%      fileName = fileName{1};
%  else
%      fileName = strsplit(dataFile, '.');
%      fileName = fileName{1};
%  end
% outputFilePrefix = [num2str(d), datestr(now,' yyyymmdd HHMM')];
outputFilePrefix =  datestr(now,' yyyymmdd HHMM');
if isfield(fiberData,'dVal')
    if isstr(fiberData.dVal)
        save([dataPath 'd' fiberData.dVal ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','totTrackLength','numCells','tracks3'); 
    else
        save([dataPath 'd' num2str(fiberData.dVal) ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','totTrackLength','numCells','tracks3'); 
    end
else
%     save([dataPath ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','totTrackLength','numCells','tracks3'); 
% save([dataPath dataFile{j}(1:6) ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','totTrackLength','numCells','tracks3'); 
% save([dataPath fileName ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','totTrackLength','numCells','tracks3'); 
% save([dataPath filename ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','totTrackLength','numCells','tracks3'); 
save([dataPath ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','totTrackLength','numCells','tracks3'); 

end
j
end



