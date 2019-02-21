% Copyright (c)2013, The Board of Trustees of The Leland Stanford Junior
% University. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
% Neither the name of the Leland Stanford Junior University nor the names
% of its contributors may be used to endorse or promote products derived
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [outputFilePrefix,logFile,logPath] = ...
    f_fitSMs_V7(dataFile,dataPath,calFile,calBeadIdx,templateFrames,peakThreshold,...
    darkFile,darkPath,logFile,logPath,boxRadius,channel, sigmaBounds,gaussianFilterSigma,minDistBetweenSMs,...
    lobeDistBounds,conversionGain,nmPerPixel,EMGain,templateLocs,threshFile,ROI,nhaData, fittingMethod)
% f_fitSMs is a module in easy_dhpsf that finds the positions of likely
% DH-PSF profiles by matching to a series of templates generated in
% f_calDHPSF and prepared in f_calSMidentification. These are then more
% precisely localized using a double gaussian fit and corrected for drift
% using the results from f_trackFiducials.
if exist('threshFile')
    load(threshFile,'usePolyROI','medianFilter','windowSize','medianBlurSigma','dispRaw','interpVal','medBlurFilt')
    meanCFocPos = zeros(120);
    if ~exist(calFile)
        [calFile, calPath] = uigetfile({'*.mat'},...
            'calFile not found! Please reselect path:');
        calFile = [calPath calFile];
    end
    
    load(calFile);
    usePolyROI = usePolyROI;
    medianFilter = medianFilter;
    windowSize = windowSize;
    medianBlurSigma = medianBlurSigma;
    medBlurFilt = medBlurFilt;
    dispRaw = dispRaw;
    interpVal = interpVal;
    dataFile_stringname='';
    
    if usePolyROI
        load(threshFile,'FOVmask','FOVmask1');
    end
end
if nhaData
    peakThreshold=peakThreshold*10000; %cancels out the call to divide in easy_dhpsf
end

printOutputFrames = 0;
if printOutputFrames == 1 % this will save all process images (correlation image, raw data, and reconstruction)
    mkdir('output images');
end

% initialize parameters
scrsz = get(0,'ScreenSize');
conversionFactor = conversionGain/EMGain;

ampRatioLimit = 0.5;
sigmaRatioLimit = 0.4;

% Options for lsqnonlin
options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
%    'FinDiffType','central','DerivativeCheck','on');
fileInfo = '';



%% preprocessing for files
%dataFile1 = dataFile;
%% need a better way to deal with dataFile1?
if isstr(dataFile)
    isTif = strsplit(dataFile, '.');
    strIsTif = isTif{end};
    isDcimg =strcmpi(strIsTif, 'dcimg');
else
    isTif = strsplit(dataFile{1}, '.');
    strIsTif = isTif{end};
    isDcimg =strcmpi(strIsTif, 'dcimg');
end

dataFile_stringname = dataFile{1}((length(dataPath)+1):length(dataFile{1}));
outputFilePrefix = cell(1,length(dataFile_stringname));

% can comment the code below in if it is desirable to select specific files
% see f_calSMidentification for how to do this

% if length(dataFile) > 1
%     dlg_title = 'Select Files';
%     prompt = {  'Choose files for thresholding' };
%     def = {num2str(1:length(dataFile))};
%     num_lines = 1;
%     inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
%     if isempty(inputdialog)
%         error('User cancelled the program')
%     end
%     selectedFiles = str2num(inputdialog{1});
% else
%     selectedFiles = 1:length(dataFile);
% end

selectedFiles = 1:length(dataFile);

% allows user to select subset of frames, at the beginning, for all files
%dlg_title = 'Select Frames (default: use all frames)';

%%This is not the best way to do this, but I am being lazy and want
%%to make sure everything else works;
absLocs = 0;
goodFit_b = zeros(1,120);
goodFit_f = zeros(1,120);
meanAmp1= zeros(1,120);
meanAmp2= zeros(1,120);
meanAmpRatio= zeros(1,120);
meanAngles= zeros(1,120);
meanInterlobeDistance= zeros(1,120);
meanPhotons= zeros(1,120);
meanX= zeros(1,120);
meanX_absolute= zeros(1,120);
meanY= zeros(1,120);
meanY_absolute= zeros(1,120);
stdAmpRatio= zeros(1,120);
stddevAngles= zeros(1,120);
stddevPhotons= zeros(1,120);
stdInterlobeDistance= zeros(1,120);
stdX= zeros(1,120);
stdY= zeros(1,120);
xAngleZero =0;
yAngleZero =0;
z= zeros(1,120);
zAngleZero =0;
ROIsteps = [];
ROIstepsize = [];
framesatstep = [];

%%%%%%%%%%%%%%% this comment and the one below can be uncommented to allow
%%%%%%%%%%%%%%% selection of specific frames: however, this is *slow*
% populates 'fileInfo' and 'numFrames' for all files, and generates the
% fields for the frame selection dlg
fileInfoAll=cell(length(selectedFiles),1);
numFramesAll = zeros(length(selectedFiles),1);

def = cell(length(selectedFiles)+4,1);
prompt = cell(length(selectedFiles)+4,1);
dlg_title = 'Input:';
num_lines = 1;

for i = 1:length(selectedFiles)
    %     fileInfoAll{i} = imfinfo([dataPath dataFile{selectedFiles(i)}]);
    %     numFramesAll(i) = length(fileInfoAll{i});
    if isDcimg
        dcimgfile = fullfile(dataPath, dataFile_stringname);
        dcimgfile = strrep(dcimgfile, '\', '\\');
        [framedata,totalframes]= dcimgmatlab(0, dcimgfile);
        totalframes = double(totalframes);
        framedatatrans = transpose (framedata);
    else
        %         fileInfoAll{i} = imfinfo([dataPath dataFile{selectedFiles(i)}]);
        fileInfoAll{i} = imfinfo([dataFile{selectedFiles(i)}]);
        numFramesAll(i) = length(fileInfoAll{i});
        totalframes = numFramesAll(i);
    end
    def{i} = ['[1:' num2str(totalframes) ']'];
    prompt{i} = ['Frames for ' dataFile{selectedFiles(i)}];
end
def{end-3} = 'No';
prompt{end-3} = 'Do you want to restart a run?';
def{end-2} = 'No';
prompt{end-2} = 'Do you want to estimate the laser profile?';
def{end-1} = 'No';
prompt{end-1} = 'Do you want to try to find true lobe sigmas?';
def{end} = 'No';
prompt{end} = 'Would you like to do the spatial correction?';
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);


%%%%%%%%%%%%%%%** this must be uncommented if selecting frames **
for i = 1:length(selectedFiles)
    framesAll{i} = str2num(inputdialog{i});
end
reloading = ~strcmp(inputdialog{end-3},'No');
findLaserInt = ~strcmp(inputdialog{end-2},'No');
findTrueSigma = ~strcmp(inputdialog{end-1},'No');
spatialCorrection = ~strcmp(inputdialog{end},'No');

if spatialCorrection
    %function [PSFfits] = spatialcorrection(PSFfits,calFile,type,calBeadIdx)
    [correctionFile, correctionPath] = uigetfile({'*.mat';'*.*'},'MultiSelect','on','Open .mat spatial calibration file');
    if correctionFile == 0
        error('Spatial correction failed because no NHA data was selected');
    else
        %         load ([correctionPath,correctionFile{1}]);
        load ([correctionPath,correctionFile]);
        
        % The NHA data we load here has been edited. Nan values,
        % shifting angle-z curves are deleted. Only
        % x,y,z,angle,x shifts and y shifts of NHA are
        % saved in the loaded file.
    end
end

sigSearchBounds = [0.9 4];
sigOptions = optimset('FunValCheck','on','Diagnostics','off','Jacobian','off', 'Display', 'off');

if findTrueSigma % this uses constrained sigmas for the fit of location, etc., and then separately finds the sigmas
    sigSearchBounds = [0.9 4];
    sigOptions = optimset('FunValCheck','on','Diagnostics','off','Jacobian','off', 'Display', 'off');
end
% This information should be passed from the earlier execution of
% this code in f_calSMidentification, but we get another chance to
% select it if it wasn't specified
if isequal(logFile,0)
    [logFile, logPath] = uigetfile({'*.dat';'*.*'},...
        'Open sequence log file(s) corresponding to image stack(s) (optional: hit cancel to skip)',...
        'MultiSelect', 'on');
    if isequal(logPath,0)
        logFile = 'not specified';
    end
    if ischar(logFile)
        logFile = cellstr(logFile);
    end
end

if reloading
    [valsFile, valsPath] = uigetfile({'*.mat';'*.*'},...
        'Would you like to restart a fitSM run?');
    if ~isequal(valsPath,0)
        load([valsPath valsFile],'PSFvals', 'meanBGvals', 'meanSignalvals');
    end
end
% sets absolute frame number for relating data to sequence log
absFrameNum = 1;

%profile on;

%% begin fitting loop over files
for stack = selectedFiles % = 1:length(dataFile)
    
    %%% this is an efficient way to set these variables if frames were
    %%% selected above
    fileIdx = find(selectedFiles == stack);
    %         fileInfo = fileInfoAll{fileIdx};
    frames = framesAll{fileIdx};
    numFrames = double(totalframes);
    
    %%% these variables must be set this way if frames were not selected
    %     fileInfo = imfinfo([dataPath dataFile{stack}]);
    %     numFrames = length(fileInfo);
    %     frames = 1:numFrames;
    if isDcimg
        [imgHeight,imgWidth] = size(framedatatrans);
    else
        dataFileInfo = fileInfoAll{fileIdx};
        imgHeight = dataFileInfo.Height;
        imgWidth = dataFileInfo.Width;
    end
    %% create output log filenames
    
    % saves in labeled directory if a channel is selected
    if channel == '0'
        if isDcimg
            outputFilePrefix{stack} = [dataFile{stack}(1:length(dataFile{stack})-6) filesep 'molecule fits  ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        else
            outputFilePrefix{stack} = [dataFile{stack}(1:length(dataFile{stack})-6) filesep 'molecule fits  ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        end
    else
        if isDcimg
            outputFilePrefix{stack} = [dataFile{stack}(1:length(dataFile{stack})-6) filesep channel(1) ' molecule fits  ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        else
            outputFilePrefix{stack} = [dataFile{stack}(1:length(dataFile{stack})-6) filesep channel(1) ' molecule fits  ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        end
    end
    mkdir(outputFilePrefix{stack});
    
    iROI = [1,1];
    if channel == 'g'
        %     ROI = imrect(gca,[1 1 64 64]);
        load('gainGreen.mat', 'gain');
        %load('offsetGreen.mat', 'BG');
        load('VarGreen.mat', 'ReadN_DI');
    elseif channel == 'r'
        %     ROI = imrect(gca,[234 234 64 64]);
        %load('offsetRed.mat', 'BG');
        load('gainRed.mat', 'gain');
        load('VarRed.mat', 'ReadN_DI');
    else
        %     ROI = imrect(gca,[1 1 64 64]);
        load('gainRed.mat', 'gain');
        %load('offsetRed.mat', 'BG');
        load('VarGreen.mat', 'ReadN_DI');
        
        %     ROI = imrect(gca,[700 700 600 600]);
    end
    %     gain = slope;
    %     offset = BG;
    gain = gain;
%     variance = ReadN_DI;
    variance = ReadN_DI.^2;

    if ~isequal(size(gain),[imgHeight imgWidth])
        warning('Gain and data image stack are not the same size. Please update initial coordinates.');
        
        dlg_title = 'Enter initial coordinates';
        prompt = { 'XO','YO'};
        def = { '1','1'};
        num_lines = 1;
        inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
        
        iROI(1) = str2double(inputdialog{1});
        iROI(2) = str2double(inputdialog{2});
        
        %         if or(size(gain)< [iROI(1)+ imgWidth-1,iROI(2)+ imgHeight-1])
        %             error('Out of bounds!')
        %         end
        gain = gain(iROI(2):iROI(2)+imgHeight-1,iROI(1):iROI(1)+imgWidth-1);
        variance = variance(iROI(2):iROI(2)+imgHeight-1,iROI(1):iROI(1)+imgWidth-1);
    end
    if stack == 1 %%% begin init
        %% miscellaneous bookkeeping on the darkfile / template / data
        % Compute darkAvg counts
        
        
        if ~isequal(darkFile,0)
            % Computes average of darkAvg frames for background subtraction
            if isDcimg
                darkFile1 = [darkPath darkFile];
            else
                darkFile1 = darkFile;
                darkFileInfo = imfinfo([darkPath darkFile]);
            end
            % Computes average of darkAvg frames for background subtraction
            %             darkFileInfo = imfinfo(darkFile);
            if isDcimg
                dcimgfile_D = fullfile(darkPath, darkFile);
                
                dcimgfile_D = strrep(dcimgfile_D, '\', '\\');
                
                [framedata_D,totalframes_D]= dcimgmatlab(0, dcimgfile_D);
                framedatatrans_D = transpose(framedata_D);
                numDarkFrames = totalframes_D;
                
                [darkHeight, darkWidth] = size(framedatatrans_D);
            else
                numDarkFrames = length(darkFileInfo);
                darkHeight= darkFileInfo(1).Height;
                darkWidth = darkFileInfo(1).Width;
            end
            darkAvg = zeros(darkHeight, darkWidth);
            for frame = 1:numDarkFrames
                if isDcimg
                    [framedata_D,totalframes_D]= dcimgmatlab(frame-1, dcimgfile_D);
                    framedatatrans_D = transpose (framedata_D);
                    darkAvg = darkAvg + double(framedatatrans_D);
                else
                    darkAvg = darkAvg + double(imread([darkPath darkFile],frame,'Info', darkFileInfo));
                end
            end
            darkAvg = darkAvg/double(numDarkFrames);
            if ~isequal(size(darkAvg),[imgHeight imgWidth])
                if isequal(size(darkAvg),size(gain))
                    warning('Dark count image and data image stack are not the same size. Chopping out region...');
                    darkAvg = darkAvg(iROI(2):iROI(2)+imgHeight-1,iROI(1):iROI(1)+imgWidth-1);
                else
                    warning('Dark count image and data image stack are not the same size. Resizing dark count image...');
                    darkAvg = imresize(darkAvg,[imgHeight imgWidth]);
                end
            end
        else
            darkAvg = 0;
        end
        if ~isDcimg
            clear darkFileInfo;
            fileInfo=imfinfo(dataFile{stack});
        end
        % Compute average image
        
        %         avgImg = zeros(imgHeight,imgWidth);
        %         for a = 1:200
        %             avgImg = avgImg + double(imread([dataPath dataFile{stack}],a,'Info',fileInfo)) - darkAvg;
        %             %    avgImg = avgImg + double(imread([dataPath dataFile],a)) - darkAvg;
        %             %Deleted 'Info' -AC 6/21
        %         end
        %         avgImg = avgImg/200;
        
        %define variables related to the templates
        %         if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
        %             templateInfo = imfinfo(templateFile);
        %             if templateInfo(1).Height ~= templateInfo(1).Width
        %                 error('Template is not square');
        %             end
        %             templateSize = templateInfo(1).Height;
        %         else
        load(threshFile,'template');
        templateSize = size(template,2);
        %         end
        numTemplates = length(templateFrames);
        templateColors = jet(numTemplates);
        
        % make sure ROI is an even number of pixels: should also be done in
        % f_calSMidentification
        if mod(ROI(3),2)==1
            ROI(3) = ROI(3)-1;
        end
        if mod(ROI(4),2)==1
            ROI(4) = ROI(4)-1;
        end
        cropWidth = ROI(3);
        cropHeight = ROI(4);
        
        %% trace out FoV if computing laser profile
        
        if findLaserInt
            temp = inputdlg({'What was the laser power at the objective? (in mW)'},...
                'Input laser power',...
                1,...
                {'0'});
            powerAtObjective = str2double(temp{1})/1000;
            
        end
        
        if findLaserInt == 1;
            avgImg = zeros(imgHeight,imgWidth);
            avgImgFrames = min(200,length(frames));
            for a = 1:avgImgFrames
                if isDcimg
                    dcimgfile = fullfile(dataPath, dataFile_stringname);
                    dcimgfile = strrep(dcimgfile, '\', '\\');
                    
                    [framedata,totalframes]= dcimgmatlab(a-1, dcimgfile);
                    totalframes = double(totalframes);
                    framedatatrans = transpose (framedata);
                    laserBkgnd =  f_waveletBackground(framedatatrans);
                    %                     avgImg = avgImg + double(framedatatrans) - darkAvg;
                    avgImg = avgImg + laserBkgnd - darkAvg;
                    
                    
                else
                    dcimgfile = fullfile(dataPath, dataFile);
                    avgImg = avgImg + double(imread([dataFile{stack}],frames(a),'Info', fileInfo)) - darkAvg;
                end
            end
            avgImg = avgImg/avgImgFrames;
            
            
            
            %             hFOVmaskFig=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
            %             imagesc(avgImg(ROI(2):ROI(2)+ROI(4)-1, ...
            %                 ROI(1):ROI(1)+ROI(3)-1),[0 300]);
            %             axis image;colorbar;colormap hot;
            %             title('Trace out the field of view boundaries by clicking at vertices');
            %             FOVmask = roipoly;
            %             close(hFOVmaskFig);
        end
        
        %% prepare template for template matching
        
        % pad template to same size as input
        templatePad = zeros(numTemplates,cropHeight,cropWidth);
        templateFT = zeros(numTemplates,cropHeight,cropWidth);
        for a=1:numTemplates
            %
            %             if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
            %                 templatePad(a,:,:) = padarray(squeeze(template(a,:,:)),...
            %                     [(cropHeight-size(template,2))/2 ...
            %                     (cropWidth-size(template,3))/2],min(min(template(a,:,:))));
            %             else
            templatePad(a,:,:) = padarray(squeeze(template(templateFrames(a),:,:)),...
                [(cropHeight-size(template,2))/2 ...
                (cropWidth-size(template,3))/2],min(min(template(templateFrames(a),:,:))));
            %             end
            
            % multiplying by conjugate of template in FT domain is equivalent
            % to flipping the template in the real domain
            templateFT(a,:,:) = conj(fft2(squeeze(templatePad(a,:,:))));
        end
        clear templatePad;
        
        % apply Gaussian filter to phase correlation data to weight low frequencies
        % more heavily since SNR is higher there
        gaussianFilter = abs(fft2(fspecial('gaussian', [cropHeight cropWidth], gaussianFilterSigma)));
        
        
        
    end %%% end init
    if nhaData
        minDistBetweenSMs = 25;
        %         boxRadius = 10;
    end
    
    %% Identify frames to analyze based on the sequence log
    % load data and register sequence log to data frames
    if ~isequal(logPath,0)
        if length(logFile) == length(dataFile) % if one sif file/data file
            sifLogData =  importdata([logPath logFile{stack}]);
            sifLogData = sifLogData(1:numFrames,:);
        else % i.e., if only one sif file for all data files
            sifLogData =  importdata([logPath logFile{1}]);
            sifLogData = sifLogData(absFrameNum:absFrameNum+numFrames-1,:);
            absFrameNum = absFrameNum + numFrames;
        end
        
        dlg_title = 'Filtering control points?';
        prompt = {  'If control points: how many stationary frames for each position? (Cancel if not CP)',...
            };
        def = {    '20', ...
            };
        num_lines = 1;
        inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
        if ~isempty(inputdialog)
            maxNumMeasurement = str2double(inputdialog{1});
            filterCPs = true;
        else
            filterCPs = false;
        end
        if filterCPs   %%% Option 1: Select only stationary steps (for CP)
            % works by selecting periods between the '-1' markers that signify motion
            % if the length of the period != step dwell time (i.e., = translation period), skip
            %
            motionFrames = find(sifLogData(:,1)==-1);  % This finds -1 entries in the shutters correspoding to moving frames
            validPeriod = diff(motionFrames)==maxNumMeasurement+1;
            validFrames = [];
            validMotionFrames = motionFrames(validPeriod);
            
            for i = 1:length(validMotionFrames)
                initFrame = validMotionFrames(i);
                temp = (initFrame+3:initFrame+maxNumMeasurement)';
                validFrames = [validFrames; temp];
            end
            frames = intersect(validFrames,frames);
        else   %%% Option 2: Select only frames that are illuminated (SMACM 2-color)
            if channel == 'g' %&& ~filterCPs   % use an intersect in case frames are limited by user
                frames = intersect(find(sifLogData(:,2) == 1),frames);
            elseif channel == 'r' %&& ~filterCPs
                frames = intersect(find(sifLogData(:,3) == 1),frames);
            end
        end
    end
    
    %% do template matching
    
    
    hSMFits=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    %totalPSFfits=[frame#, loc# w/in frame, xLocation yLocation matchingTemplateNumber matchConfidence (6),
    %amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2 bkgndMean(15)
    %totalFitError goodFit xCenter yCenter angle numPhotons
    %interlobeDistance, amplitude ratio, sigma ratio(24)
    %corrected X, corrected Y, corrected Z (27)]
    %     totalPSFfits = zeros(12000000, 6+18+3);
    numPSFfits = 0;
    startTime = tic;
    logFlag = 0;
    cropWidth = ROI(3);
    cropHeight = ROI(4);
    bkgndImgTotal = zeros(length(ROI(2):ROI(2)+ROI(4)-1),...
        length(ROI(1):ROI(1)+ROI(3)-1));
    %     bkgndImgTotalLaser = zeros(length(ROI(2):ROI(2)+ROI(4)-1),...
    %         length(ROI(1):ROI(1)+ROI(3)-1));
    numbkgndImg = 0;
    
    if size(frames,1) > 1 % make sure it will work in the for loop (need 1xn)
        frames = frames';
    end
    
    
    % set up mean BG for laser spot calculation
    meanBG = nan(length(frames),1);
    meanSignal = nan(length(frames),1);
    if ~exist('PSFvals')
        PSFvals = cell(length(frames),1);
        meanBGvals = cell(length(frames),1);
        meanSignalvals = cell(length(frames),1);
    end
    %     meanBGvals = cell(length(frames),1);
    %     meanSignalvals = cell(length(frames),1);
    
    bkgndImgvals = cell(idivide(int32(length(frames)),int32(interpVal))+1,1);
    %    bkgndImgvals = cell(length(frames),1);
    %%makeshift, there is a more clever way to do this.
    usePolyROI = usePolyROI;
    FOVmask = FOVmask;
    FOVmask1 = FOVmask1;
    medianFilter = medianFilter;
    windowSize = windowSize;
    medianBlurSigma = medianBlurSigma;
    dispRaw = dispRaw;
    interpVal = interpVal;
    trueFOVmask =0;
    bkgndImgvals = cell(length(frames),1);
    if exist('FOVmask')
        trueFOVmask =1;
    end
    
    %     if (imgHeight < 2048) || (imgWidth < 2048)
    %                         prompt = {'x coordinate?','y coordinate?'};
    %                         dlg_title = 'FOV coordinates';
    %                         num_lines = 1;
    %                         def = {'0','0'};
    %                         coordinates = inputdlg(prompt,dlg_title,num_lines,def);
    %                         xcoordinate = str2num(coordinates{1});
    %                         ycoordinate = str2num(coordinates{2});
    %                         xFOV(1:2) = [xcoordinate, xcoordinate + imgWidth-1];
    %                         yFOV(1:2) = [ycoordinate, ycoordinate + imgHeight-1];
    %                         gain = gain(yFOV(1):yFOV(2),xFOV(1):xFOV(2));
    %     end
    
    dimMask = size(FOVmask);bufMask = 5;
    sMask = false(dimMask(1)+2*bufMask,dimMask(2)+2*bufMask,4);
    sMask(1+2*bufMask:dimMask(1)+2*bufMask,1+bufMask:dimMask(2)+bufMask,1) = FOVmask;
    sMask(1:dimMask(1),1+bufMask:dimMask(2)+bufMask,2) = FOVmask;
    sMask(1+bufMask:dimMask(1)+bufMask,1+2*bufMask:dimMask(2)+2*bufMask,3) = FOVmask;
    sMask(1+bufMask:dimMask(1)+bufMask,1:dimMask(2),4) = FOVmask;
    FTFOVmask = FOVmask & sMask(1+bufMask:dimMask(1)+bufMask,1+bufMask:dimMask(2)+bufMask,1) & sMask(1+bufMask:dimMask(1)+bufMask,1+bufMask:dimMask(2)+bufMask,2) & sMask(1+bufMask:dimMask(1)+bufMask,1+bufMask:dimMask(2)+bufMask,3) & sMask(1+bufMask:dimMask(1)+bufMask,1+bufMask:dimMask(2)+bufMask,4);
    
    
    datavar = double((variance./gain)./gain);
    
    datavarCrop = datavar(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
    gain = gain(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
    ROIinitial = ROI;
    if ~isempty(frames)
        
        if medianFilter
            bkgndImgvals = cell(idivide(int32(length(frames)), int32(interpVal)) + 1,1);
            % set up median filter
            dataWindow = nan([2*floor(windowSize/2)+1,size(bkgndImgTotal)]);
            if medianBlurSigma~=0
                medBlurFilt = fspecial('gaussian',100,medianBlurSigma);
            else
                medBlurFilt=1;
            end
            lastFrame = nan;
            finalframe = min(100, length(frames));
            for c = frames(finalframe:-interpVal:1)
                currIdx = find(frames==c);
                if isnan(lastFrame)
                    if isDcimg
                        [bkgndImgMed, dataWindow] = f_medianFilter_DCImg(dcimgfile, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                    else
                        [bkgndImgMed, dataWindow] = f_medianFilter(dataFile{stack}, fileInfo, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                    end
                    lastFrame = c;
                    bkgndImgvals{idivide(int32(c), int32(interpVal))+1}= bkgndImgMed;
                    
                    bkgndImgvals{idivide(int32(c), int32(interpVal))}= bkgndImgMed;
                else
                    if isDcimg
                        [bkgndImgMed, dataWindow] = f_medianFilter_DCImg(dcimgfile, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                    else
                        [bkgndImgMed, dataWindow] = f_medianFilter([dataPath dataFile{stack}], fileInfo, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                    end
                    lastFrame =c;
                    bkgndImgvals{idivide(int32(c), int32(interpVal))}= bkgndImgMed;
                end
            end
            if isDcimg
                clear mex;
            end
        end
    end
    
    for c = 1:min(100,length(frames))
        if ~isempty(PSFvals{c,1})
            continue
        end
        if isDcimg
            [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dcimgfile,medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
        else
            [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dataFile{stack},medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
        end
        bkgndImgTotal = bkgndImgTotal + bkgndImg;
        numbkgndImg = numbkgndImg +1;
        
        dataFT = fft2(data,cropHeight,cropWidth);
        [PSFLocs,numPSFLocs,maxPeakImg] = GenerateLocs(cropHeight,cropWidth,gaussianFilter,dataFT,templateFT,FOVmask,peakThreshold(stack,:),numTemplates,minDistBetweenSMs,nhaData,trueFOVmask,FTFOVmask);
        if isempty(PSFLocs)
            continue
        end
        meanSignalvals{c,1} = mean(data(FOVmask));
        meanBGvals{c,1} = mean(bkgndImg(FOVmask));
        PSFfits = FittingLocs(PSFLocs,FOVmask1,ROI,boxRadius,fittingMethod,cropWidth,cropHeight,templateLocs,templateSize,data,datavarCrop,bkgndImg,sigmaBounds,conversionFactor,nmPerPixel,sigSearchBounds,sigmaRatioLimit,lobeDistBounds,ampRatioLimit,findTrueSigma,iROI);
        PlotImage(PSFfits,PSFLocs,cropWidth,cropHeight,ROI,data,bkgndImg,dispRaw,printOutputFrames,numPSFLocs,peakThreshold,maxPeakImg,templateColors,frames(c),stack,hSMFits,iROI);
        PSFvals{c,1} = [frames(c)*ones(numPSFLocs,1), (1:numPSFLocs)', PSFLocs, PSFfits];
        numPSFfits = numPSFfits+numPSFLocs;
        
        
        if isDcimg
            clear mex;
        end
    end
    

    restartframe = 101;
    if reloading
        tfr= frames(~cellfun(@isempty,PSFvals));
        restartframe = round(tfr(end)/100)*100+1;
    end
    
    if ~isDcimg
        dcimgfile = [];
    end
    if length(frames) > 100 && isempty(PSFvals{end,1})
        if (floor(length(frames)/100)*100) >= 100
            for n = restartframe:100:length(frames)
                p = n + 99;
                p = min(p, length(frames));
                
                if medianFilter
                    bkgndImgvals = cell(idivide(int32(length(frames)), int32(interpVal))+1,1);
                    dataWindow = nan([2*floor(windowSize/2)+1,size(bkgndImgTotal)]);
                    
                    lastFrame = nan;
                    for c = frames(p:-interpVal:n)
                        currIdx = find(frames==c);
                        if isnan(lastFrame)
                            if isDcimg
                                [bkgndImgMed, dataWindow] = f_medianFilter_DCImg(dcimgfile, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                            else
                                [bkgndImgMed, dataWindow] = f_medianFilter([dataFile{stack}], fileInfo, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                            end
                            lastFrame = c;
                            bkgndImgvals{idivide(int32(c), int32(interpVal))+1}= bkgndImgMed;
                            
                            bkgndImgvals{idivide(int32(c), int32(interpVal))}= bkgndImgMed;
                        else
                            if isDcimg
                                [bkgndImgMed, dataWindow] = f_medianFilter_DCImg(dcimgfile, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                            else
                                [bkgndImgMed, dataWindow] = f_medianFilter([dataFile{stack}], fileInfo, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                            end
                            lastFrame =c;
                            bkgndImgvals{idivide(int32(c), int32(interpVal))}= bkgndImgMed;
                        end
                    end
                    if isDcimg
                        clear mex;
                    end
                end
                parfor c = n:p
                    if isDcimg
                        [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dcimgfile,medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
                    else
                        [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dataFile{stack},medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
                    end
                    bkgndImgTotal = bkgndImgTotal + bkgndImg;
                    numbkgndImg = numbkgndImg +1;
                    dataFT = fft2(data,cropHeight,cropWidth);
                    [PSFLocs,numPSFLocs,~] = GenerateLocs(cropHeight,cropWidth,gaussianFilter,dataFT,templateFT,FOVmask,peakThreshold(stack,:),numTemplates,minDistBetweenSMs,nhaData,trueFOVmask,FTFOVmask);
                    if isempty(PSFLocs)
                        continue
                    end
                    meanSignalvals{c,1} = mean(data(FOVmask));
                    meanBGvals{c,1} = mean(bkgndImg(FOVmask));
                    PSFfits = FittingLocs(PSFLocs,FOVmask1,ROI,boxRadius,fittingMethod,cropWidth,cropHeight,templateLocs,templateSize,data,datavarCrop,bkgndImg,sigmaBounds,conversionFactor,nmPerPixel,sigSearchBounds,sigmaRatioLimit,lobeDistBounds,ampRatioLimit,findTrueSigma,iROI);
                    %                         PlotImage(PSFfits,PSFLocs,cropWidth,cropHeight,ROI,data,bkgndImg,dispRaw,printOutputFrames,fittingMethod,numPSFLocs,peakThreshold,maxPeakImg,templateColors,frames(c),stack,hSMFits);
                    PSFvals{c,1} = [frames(c)*ones(numPSFLocs,1), (1:numPSFLocs)', PSFLocs, PSFfits];
                    numPSFfits = numPSFfits+numPSFLocs;
                end
                if isDcimg
                    clear mex;
                end
                pcalled = ['Frame Number: ' num2str(p)];
                if rem(p, 1000) == 0
                    save('-v7.3', [outputFilePrefix{stack} 'partial molecule fits.mat'],'PSFvals', 'meanBGvals', 'meanSignalvals');
                end
                disp(pcalled);
            end
        end
    end
    elapsedTime = toc(startTime);
    ROI = ROIinitial;
    datavarCrop = datavar(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
    %         totalPSFfit1 = []; %Initializing this to a cell array takes slightly longer to complete processing.
    totalPSFfits = [];
    meanBG = zeros(length(meanBGvals),1);
    meanSignal = zeros(length(meanSignalvals),1);
    totalPSFfits = cell2mat(PSFvals);
    totalPSFfits(:,25:27) = zeros(size(totalPSFfits,1),3);
    meanBG(~cellfun(@isempty, meanBGvals)) = cell2mat(meanBGvals);
    meanSignal(~cellfun(@isempty, meanSignalvals)) = cell2mat(meanSignalvals);
    
    %         meanBG(:,1) = meanBG1(:,2);
    %         meanSignal(:,1) = meanSignal1(:,2);
    
    %         totalPSFfits = totalPSFfits(1:numPSFfits,:);
    clear data bkgnd residual fileInfo maxPeakImg reconstructImg xIdx yIdx temp totalPSFfit2 meanBG1 meanSignal1;
    %     clear data bkgnd residual fileInfo maxPeakImg reconstructImg templateFT xIdx yIdx temp;
    
    if logFlag ~= 0
        sprintf('There were %d frames in which the peakImg matrix contained complex numbers or NaNs. See log file "peakImg log.txt" for more details',logFlag)
        logFlag = 0;
    end
    close(hSMFits)
    fps = length(frames)/elapsedTime;
    moleculesPerSec = numPSFfits/elapsedTime;
    display(['Frames per second: ' num2str(fps) '  Molecules per second: ' num2str(moleculesPerSec)]);
    %movie2avi(M,'output_v1.avi','fps',10,'Compression','FFDS');
    
    %         %% Measure the Gaussian Laser Intensity Distribution
    %         if findLaserInt == 1;
    %             % This takes the average bkgndImg found using f_waveletBackground
    %             % and tries to extract a Gaussian laser profile from it, assuming the
    %             % background intensity is uniformly proportional to laser intensity
    %
    %             %This is measured in counts
    % %             bkgndImg_avg = bkgndImgTotal/numbkgndImg;
    %               avgImg = avgImg./gain;
    %             [laser_x_nm, laser_y_nm ,sigma_x_nm, sigma_y_nm, theta, peakIntensity, waist]...
    %                 = f_findGaussianLaserProfile...
    %                 (bkgndImg_avg, FOVmask, nmPerPixel, powerAtObjective, ROI);
    
    %             laserProfile=avgImg*powerAtObjective/sum(sum(avgImg))...
    %                 /(nmPerPixel^2)*10^14;
    %             %             laserProfile=bkgndImg_avg*powerAtObjective/sum(sum(bkgndImg_avg))...
    %             %                 /(nmPerPixel^2)*10^14; %multiply image by power measured
    %divided by sum of all intensities
    %divided by area of a pixel
    %             laserProfile_interp = vertcat(laserProfile_interp,laserProfile(e,:)');
    
    
end

%% Translate angle into corrected x,y positions and z position

load(calFile);
z = meanCFocPos - meanCFocPos(60);
goodFit_forward = logical((goodFit_f(calBeadIdx,:)));
totalPSFfits(:,25) = totalPSFfits(:,18) ...
    - interp1((meanAngles(calBeadIdx,goodFit_forward)),...
    (meanX(calBeadIdx,goodFit_forward)),totalPSFfits(:,20),'spline',NaN);
totalPSFfits(:,26) = totalPSFfits(:,19) ...
    - interp1((meanAngles(calBeadIdx,goodFit_forward)),...
    (meanY(calBeadIdx,goodFit_forward)),totalPSFfits(:,20),'spline',NaN);
totalPSFfits(:,27) = interp1((meanAngles(calBeadIdx,goodFit_forward)),...
    (z(goodFit_forward)),totalPSFfits(:,20),'spline',NaN);
totalPSFfits(totalPSFfits(:,20) < min(meanAngles(calBeadIdx,goodFit_forward)) | totalPSFfits(:,20) > max(meanAngles(calBeadIdx,goodFit_forward)),25:27) = NaN;
%% NHA based correction


%         if spatialCorrection
%                 dataFrames = unique(totalPSFfits(:,1));
%                 corFrames = length(dataFrames);
%
%                 % Calculate distance and find the Nearest Neighbour
%                 Dis_NN_frame = cell(1,corFrames);
%                 val_Dis_NN = cell(1,corFrames);
%                 indx_Dis_NN = cell(1,corFrames);
%                 Corrected = cell(corFrames,1);
%                 %         Corrected_x = cell(corFrames,1);
%                 %         Corrected_y = cell(corFrames,1);
%                 %         Corrected_z = cell(corFrames,1);
%
%                 for k = 1:corFrames % Column 24 25 26 (uncorrected x y and angle) in totalPSFfit1 are 18 19 20 in totalPSFfits
%                     temp = totalPSFfits(ismember(totalPSFfits(:,1),dataFrames(k)),18:20);
%                     Dis_NN_frame{k} = sqrt((bsxfun(@minus,x_NHA,temp(:,1)')).^2 ...
%                         + (bsxfun(@minus,y_NHA,temp(:,2)').^2));
%                     [val_Dis_NN{k},indx_Dis_NN{k}] = min(Dis_NN_frame{k});
%
%                     for SMs = 1:length(indx_Dis_NN{k})
%                         Corrected{k}(SMs,1) = temp(SMs,1) - interp1(angles_NHA_good(indx_Dis_NN{k}(SMs),:),...
%                             meanX_good(indx_Dis_NN{k}(SMs),:),temp(SMs,3),'spline');
%                         Corrected{k}(SMs,2) = temp(SMs,2) - interp1(angles_NHA_good(indx_Dis_NN{k}(SMs),:),...
%                             meanY_good(indx_Dis_NN{k}(SMs),:),temp(SMs,3),'spline');
%                         Corrected{k}(SMs,3) = interp1(angles_NHA_good(indx_Dis_NN{k}(SMs),:),...
%                             -1*meanCFocPos_NHA_rel_max,temp(SMs,3),'spline');
%                     end
%                     %
%                     %             Corrected_x{k} = totalPSFfit1{k}(:,31);
%                     %             Corrected_y{k} = totalPSFfit1{k}(:,32);
%                     %             Corrected_z{k} = totalPSFfit1{k}(:,33);
%
%                 end
%                 %         totalPSFfits(:,28:30) = [cell2mat(Corrected_x),cell2mat(Corrected_y),cell2mat(Corrected_z)];
%                 totalPSFfits(:,28:30) = cell2mat(Corrected);
%                 %         totalPSFfits(:,28) = cell2mat(Corrected_x);
%                 %         totalPSFfits(:,29) = cell2mat(Corrected_y);
%                 %         totalPSFfits(:,30) = cell2mat(Corrected_z);
%         end

%% New Trial Correction Code
% if spatialCorrection
%     dataFrames = unique(totalPSFfits(:,1));
%     corFrames = length(dataFrames);
%
%     % Calculate distance and find the Nearest Neighbour
%     Dis_NN_frame = cell(1,corFrames);
%     val_Dis_NN = cell(1,corFrames);
%     indx_Dis_NN = cell(1,corFrames);
%     Corrected = cell(corFrames,1);
%     for k = 1:corFrames % Column 24 25 26 (uncorrected x y and angle) in totalPSFfit1 are 18 19 20 in totalPSFfits
%         temp = totalPSFfits(ismember(totalPSFfits(:,1),dataFrames(k)),18:20);
%         Dis_NN_frame{k} = sqrt((bsxfun(@minus,x_NHA,temp(:,1)')).^2 ...
%             + (bsxfun(@minus,y_NHA,temp(:,2)').^2));
%         [val_Dis_NN{k},indx_Dis_NN{k}] = min(Dis_NN_frame{k});
%
%         for SMs = 1:length(indx_Dis_NN{k})
%             Corrected{k}(SMs,1) = temp(SMs,1) - interp1(Corr_NHA{indx_Dis_NN{k}(SMs),1}(:,1),...
%                 Corr_NHA{indx_Dis_NN{k}(SMs),1}(:,2),temp(SMs,3),'spline');
%             Corrected{k}(SMs,2) = temp(SMs,2) - interp1(Corr_NHA{indx_Dis_NN{k}(SMs),1}(:,1),...
%                 Corr_NHA{indx_Dis_NN{k}(SMs),1}(:,3),temp(SMs,3),'spline');
%             Corrected{k}(SMs,3) = interp1(Corr_NHA{indx_Dis_NN{k}(SMs),1}(:,1),...
%                 Corr_NHA{indx_Dis_NN{k}(SMs),1}(:,4),temp(SMs,3),'spline');
%         end
%     end
%     totalPSFfits(:,28:30) = cell2mat(Corrected);
% end
% localCorrection = true; localThresh = 4000;
if spatialCorrection
    totalPSFfits(:,28:35) = NaN(size(totalPSFfits,1),8);
    dataFrames = unique(totalPSFfits(:,1));
    corFrames = length(dataFrames);
    x_NHA = x_NHA; y_NHA = y_NHA;
    % Calculate distance and find the Nearest Neighbour
    %     Dis_NN_frame = cell(1,corFrames);
    if corFrames < 4*length(x_NHA)
        val_Dis_NN = cell(1,corFrames);
        indx_Dis_NN = cell(1,corFrames);
        for k = 1:corFrames % n totalPSFfit1 are 18 19 20 in totalPSFfits
            temp = totalPSFfits(ismember(totalPSFfits(:,1),dataFrames(k)),18:19);
            Dis_NN_frame = sqrt((bsxfun(@minus,x_NHA,temp(:,1)')).^2 ...
                + (bsxfun(@minus,y_NHA,temp(:,2)').^2));
            [val_Dis_NN{k},indx_Dis_NN{k}] = min(Dis_NN_frame);
        end
        val_NN = cell2mat(val_Dis_NN);
        indxs = cell2mat(indx_Dis_NN);
    else
        val_NN = inf(size(totalPSFfits,1),1);
        indxs = ones(size(totalPSFfits,1),1);
        %         gap = 25;
        %         for k = 1:gap:length(x_NHA)
        %             %             Dis_NN = bsxfun(@hypot,totalPSFfits(:,18)-x_NHA(k),totalPSFfits(:,19)-y_NHA(k));
        %             %           Dis_NN = sqrt((totalPSFfits(:,18)-x_NHA(k)).^2 + (totalPSFfits(:,19)-y_NHA(k)).^2 );
        %             [Dis_NN, indx_Dis_NN] = min(bsxfun(@hypot,bsxfun(@minus,totalPSFfits(:,18),x_NHA(k:min(k+gap-1,length(x_NHA)))'),bsxfun(@minus,totalPSFfits(:,19),y_NHA(min(k+gap-1,length(x_NHA)))')),[],2);
        %             gNN = Dis_NN < val_NN;
        %             val_NN(gNN) = Dis_NN(gNN);
        %             indxs(gNN) = indx_Dis_NN(gNN)+k-1;
        %         end
        for k = 1:length(x_NHA)
            Dis_NN = bsxfun(@hypot,totalPSFfits(:,18)-x_NHA(k),totalPSFfits(:,19)-y_NHA(k));
            gNN = Dis_NN < val_NN;
            val_NN(gNN) = Dis_NN(gNN);
            indxs(gNN) = k;
        end
    end
    totalPSFfits(:,31:32) = [val_NN,indxs];
    %     if localCorrection
    %         for k = 1:corFrames
    %             indx_Dis_NN{k}(val_Dis_NN{k}>localThresh) = nan;
    %         end
    %     end
    
    indxU = unique(indxs(~isnan(indxs)));
    
    for SMs = 1:length(indxU)
        NHA = indxs == indxU(SMs);
        totalPSFfits(NHA,28) = totalPSFfits(NHA,18) - interp1(xy_shift_adj{indxU(SMs),1}(:,1),...
            xy_shift_adj{indxU(SMs),1}(:,2),totalPSFfits(NHA,20),'spline',NaN);
        totalPSFfits(NHA,29) = totalPSFfits(NHA,19) - interp1(xy_shift_adj{indxU(SMs),1}(:,1),...
            xy_shift_adj{indxU(SMs),1}(:,3),totalPSFfits(NHA,20),'spline',NaN);
        totalPSFfits(NHA,30) = interp1(Corr_Z{indxU(SMs),1}(:,1),...
            Corr_Z{indxU(SMs),1}(:,2),totalPSFfits(NHA,20),'spline',NaN);
        if exist('Corr_Z_L','var')
            totalPSFfits(NHA,33) = totalPSFfits(NHA,18) - interp1(xy_shift_adj_L{indxU(SMs),1}(:,1),...
                xy_shift_adj_L{indxU(SMs),1}(:,2),totalPSFfits(NHA,20),'spline',NaN);
            totalPSFfits(NHA,34) = totalPSFfits(NHA,19) - interp1(xy_shift_adj_L{indxU(SMs),1}(:,1),...
                xy_shift_adj_L{indxU(SMs),1}(:,3),totalPSFfits(NHA,20),'spline',NaN);
            totalPSFfits(NHA,35) = interp1(Corr_Z_L{indxU(SMs),1}(:,1),...
                Corr_Z_L{indxU(SMs),1}(:,2),totalPSFfits(NHA,20),'spline',NaN);
        end
    end
    %
    lL = nan(size(Corr_Z));
    uL = lL;lL_L = lL;uL_L = lL;
    for a = 1:size(Corr_Z,1)
        lL(a) = min(Corr_Z{a,1}(:,1));
        uL(a) = max(Corr_Z{a,1}(:,1));
        if exist('Corr_Z_L','var')
            lL_L(a) = min(Corr_Z_L{a,1}(:,1));
            uL_L(a) = max(Corr_Z_L{a,1}(:,1));
        end
    end
    totalPSFfits(totalPSFfits(:,20) < lL(totalPSFfits(:,32)) |...
        totalPSFfits(:,20) > uL(totalPSFfits(:,32)),28:30) = NaN;
    if exist('Corr_Z_L','var')
        totalPSFfits(totalPSFfits(:,20) < lL_L(totalPSFfits(:,32)) |...
            totalPSFfits(:,20) > uL_L(totalPSFfits(:,32)),33:35) = NaN;
    end
    %             totalPSFfits(totalPSFfits(NHA,20) < min(Corr_Z_L{indxU(SMs),1}(:,1)) | totalPSFfits(NHA,20) > max(Corr_Z_L{indxU(SMs),1}(:,1)),33:35) = nan;
    %         totalPSFfits(totalPSFfits(NHA,20) < min(Corr_Z{indxU(SMs),1}(:,1)) | totalPSFfits(NHA,20) > max(Corr_Z{indxU(SMs),1}(:,1)),28:30) = nan;
end
%% output data to external file
clear dataWindow;
save('-v7.3', [outputFilePrefix{stack} 'molecule fits.mat']);
delete([outputFilePrefix{stack} 'partial molecule fits.mat']);
%     else
%         error('No frames selected!')

% profile off
% profile viewer;
% p = profile('info');
% profsave(p,'profile_results_ffitSMs');

end
% end
