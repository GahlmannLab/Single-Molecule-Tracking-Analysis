% Tracking MSD
%   This function takes the localizations output from a single-molecule
%   experiment and links the localizations into trajectories. The
%   mean-squared-displacement (MSD) is calculated for each trajectory.

% Set isSim to 1 only if the input data has been generated from direct
% simulation of trajectories. Do not set to 1 if the simulation was carried
% out by simulation of DHPSF images.
isSim = 0;

%Set dimensionality for 2D or 3D tracking
dim = 3;


% Load in the files with the localizations sorted into their individual
% cells. Multiple files may be loaded in for data sets with multiple
% field's of view (FOV).
[dataFile, dataPath] = uigetfile({'*.mat';'*.*'},'MultiSelect', 'on');

% User inputs. Set the exposure time of the camera, the distance threshold
% for linking localizations, and the max number of frames to skip for
% localization linking.
dlg_title = 'Input info';
prompt = {'Exposure Time (ms)','Match Distance (nm)','Max Frames to Skip'};
def = {'25','2200','1'};
num_lines = 1;
inputInfo = inputdlg(prompt,dlg_title,num_lines,def);
time_per_frame = str2num(inputInfo{1});
matchdist = str2num(inputInfo{2});
frameSkip = str2num(inputInfo{3});

%Multiple files with cell localizations can be input
if iscell(dataFile)
    fileNum = length(dataFile);
else
    fileNum = 1;
end

for j = 1:fileNum
    
    %Load in file with localizations
    if iscell(dataFile)
        load ([dataPath,dataFile{1,j}]);
    else
        load ([dataPath,dataFile]);
    end
    
    %Initialize variables
    numCells = length(fiberData);
    tracks = cell(numCells,1);
    diffusionCoefficients = [];
    totTrackLength = [];
    frameRange = fiberData(1).FrameRange;
    corrected_xyz = [];
    
    %Loop through each cell. The localization linking is performed one cell
    %at a time to prevent linking localizations in adjacent cells.
    for c = 1:numCells
        %The fiberData variable is created in outfti. It contains
        %information on the cell outlines and localizations within the cell
        %outlines.
        corrected_xyz(:,1) = fiberData(c).xLoc;
        corrected_xyz(:,2) = fiberData(c).yLoc;
        corrected_xyz(:,3) = fiberData(c).zLoc;
        frameNum = fiberData(c).Frame;
        points = cell(frameRange(1,2),1);
        for a = 1:frameRange(1,2)
            points{a,1} = corrected_xyz(frameNum(:,1)==a,:);
        end
        
        %Do the localization tracking
        if isSim == 1
            for j = 1:length(corrected_xyz)/6 %number of tracks (simulated tracklength of 6)
                tracks{j} = zeros(6,4);
                tracks{j}(:,1) = [25:25:150];
                tracks{j}(:,2:4) = corrected_xyz((j-1)*6+1:j*6,:);
            end
        else
            if dim == 3
            %The simpletracker function was modified from the original code
            %by Jean-Yves Tinevez.
            [temptracks, tracks_no_t_c, tracks_time_coord] = simpletracker(points,time_per_frame, 'Method','NearestNeighbor','MaxLinkingDistance',matchdist,'MaxGapClosing',1,'Debug',false);
            elseif dim == 2
                [tracks, tracks_no_t_c, tracks_time_coord] = simpletracker_2D(points,time_per_frame, 'Method','NearestNeighbor','MaxLinkingDistance',matchdist,'MaxGapClosing',1,'Debug',false);

            end
        end
        
        
        %Comput the diffusion coefficient for each trajectory by
        %calculating the mean squared displacement (MSD) between localizations
        tempdiff = zeros(length(temptracks),1);
        trackLength = zeros(length(temptracks),1);
        for e = 1:length(temptracks)
            single_frame_sq_displacements = [];
            for f = 1:length(temptracks{e}(:,1)) - 1
                single_frame_sq_displacements(f,1) = sum((temptracks{e}(f,2:dim+1) - temptracks{e}(f+1,2:dim+1)).^2,2);
            end
            tempdiff(e,1) = mean(single_frame_sq_displacements)*10^-6/((2*dim)*(time_per_frame*10^-3));
            trackLength(e) = length(temptracks{e});
        end
        diffusionCoefficients = vertcat(diffusionCoefficients,tempdiff);
        totTrackLength = vertcat(totTrackLength,trackLength);
        tracks{c} = temptracks;
        c
        clear corrected_xyz points tracks_no_t_c tracks_time_coord single_frame_sq_displacements
    end
    
    
    %Delete tracks and diffusion coefficients that have 2 or more points in
    %a cell during the same frame
    diffusionCoefficientsAll = diffusionCoefficients;
    tracksAll = tracks;
    frameRange = cell(length(tracks),1);
    cropIdxLoc = cell(length(tracks),1);
    cropIdxLocTot = [];
    for c = 1:length(tracks)
        frameRange{c} = zeros(length(tracks{c}),2);
        cropIdxLoc{c} = ones(length(tracks{c}),1);
        for l = 1:length(tracks{c})
            frameRange{c}(l,1) = tracks{c}{l}(1,1)/time_per_frame + 1;
            frameRange{c}(l,2) = tracks{c}{l}(end,1)/time_per_frame + 1;
            tempRange = frameRange{c}(l,1):1:frameRange{c}(l,2);
            for m = tempRange(1):tempRange(end)
                if sum(fiberData(c).Frame == m) > 1
                    cropIdxLoc{c}(l,1) = 0;
                end
            end
        end
        cropIdxLocTot = vertcat(cropIdxLocTot,cropIdxLoc{c});
        tracks{c} = tracks{c}(logical(cropIdxLoc{c}));
    end
    diffusionCoefficients = diffusionCoefficientsAll(cropIdxLocTot == 1);
    
    
    
    %Create filename and save the relevent variables
    outputFilePrefix =  datestr(now,' yyyymmdd HHMM');
    if isfield(fiberData,'dVal')
        if isstr(fiberData.dVal)
            save([dataPath 'd' fiberData.dVal ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','diffusionCoefficientsAll','numCells','tracks','tracksAll');
        else
            save([dataPath 'd' num2str(fiberData.dVal) ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','diffusionCoefficientsAll','numCells','tracks','tracksAll');
        end
    else
        save([dataPath ' Diffusion Coefficients ' outputFilePrefix '.mat'],'diffusionCoefficients','diffusionCoefficientsAll','numCells','tracks','tracksAll');
    end
end



