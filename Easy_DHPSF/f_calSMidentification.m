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
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [templateFrames, ROI, dataFile1, dataPath, darkFile,darkPath,...
    EMGain, templateLocs, outputFilePrefix,nhaData,tempThresh, iROI] = f_calSMidentification(calFile,calBeadIdx,...
    templateFile, boxRadius,channel,sigmaBounds,gaussianFilterSigma,minDistBetweenSMs, fittingMethod)
% f_calSMidentification is a module in easy_dhpsf that prepares the
% templates from f_calDHPSF and uses them to generate a series of template
% matches. These are then used to judge an appropriate threshold for
% f_fitSMs. This module also sets the file and other parameters used for
% f_fitSMs, as well as some parameters for f_trackFiducials.
% if this is true, select rect ROI for channel and poly for excluding fids,
% edge of FoV outside iris, any other impinging 'non-image' feature
usePolyROI = true;
% profile on

% initializing variables for ROI selection
FOVmask = [];

%Ask if user wants to use SSIM to estimate thresholds. Use this to
%automatically estimate thresholds. The ssimThresh vairable should be
%adjusted to get satisfactory thresholds. The output threshold values may 
%still need to be manually adjusted.
dlg_title = 'SSIM';
prompt = {'Do you want to use SSIM to estimate thresholds?'};
def = {num2str(1)};
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
useSSIM = str2num(inputdialog{1});
ssimThresh = 0.04;

% Instrument Specific Parameters

dlg_title = 'Set EM Gain';
prompt = {  'EM Gain (1 if no gain):' };
def = {'1'};

meanCFocPos = zeros(120);
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
if ~exist(calFile)
    [calFile, calPath] = uigetfile({'*.mat'},...
        'calFile not found! Please reselect path:');
    calFile = [calPath calFile];
end
load(calFile);

num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);

if isempty(inputdialog)
    error('User cancelled the program')
end

EMGain = str2double(inputdialog{1});
if EMGain < 1 || isnan(EMGain)
    warning('EMGain should be >= 1. Setting to 1...');
    EMGain = 1;
end

frameNum = 1;
scrsz = get(0,'ScreenSize');
% Options for lsqnonlin
options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
%    'FinDiffType','central','DerivativeCheck','on');

%% ask user for relevant datafiles
[dataFile, dataPath] = uigetfile({'*.dcimg'; '*.tif'},...
    'Open SMACM image stack(s) for data processing',...
    'MultiSelect', 'on');
if isequal(dataFile,0)
    error('User cancelled the program');
end

if isstr(dataFile)
    isTif = strsplit(dataFile, '.');
    strIsTif = isTif{2};
    isDcimg =strcmpi(strIsTif, 'dcimg');
else
    isTif = strsplit(dataFile{1}, '.');
    strIsTif = isTif{2};
    isDcimg =strcmpi(strIsTif, 'dcimg');
end

if (isDcimg)
    if ischar(dataFile)==1
        dataFile1 = [dataPath dataFile];
        dataFile1 = {dataFile1};
    end
else
    if ischar(dataFile)==1
        dataFile = {dataFile};
    end
    for i = 1:length(dataFile)
        dataFile{i}= [dataPath dataFile{i}];
    end
    dataFile1 = dataFile;
end

% allows user to limit the number of files in case there are many, many
% large files (i.e., very long acquisitions where the sample does not
% change very much)
if length(dataFile1) > 1
    dlg_title = 'Select Files';
    prompt = {  'Choose files for thresholding' };
    def = {num2str(1:length(dataFile))};
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(inputdialog)
        error('User cancelled the program')
    end
    selectedFiles = str2num(inputdialog{1});
else
    selectedFiles = 1:length(dataFile1);
end

% allows user to select subset of frames, at the beginning, for all files
dlg_title = 'Select Frames (default: use all frames)';
num_lines = 1;
def = {};
prompt = {};
% populates 'fileInfo' and 'numFrames' for all files, and generates the
% fields for the frame selection dlg
for i = 1:length(selectedFiles)
    %     fileInfoAll{i} = imfinfo([dataPath dataFile{selectedFiles(i)}]);
    %     numFramesAll(i) = length(fileInfoAll{i});
    if  (isDcimg)
        dcimgfile = fullfile(dataPath, dataFile);
        dcimgfile = strrep(dcimgfile, '\', '\\');
        [framedata,totalframes]= dcimgmatlab(0, dcimgfile);
        totalframes= double(totalframes);
        framedatatrans = transpose(framedata);
    else
        fileInfoAll{i} = imfinfo([dataFile{selectedFiles(i)}]);
        numFramesAll(i) = length(fileInfoAll{i});
        totalframes = numFramesAll(i);
    end
    def{i} = ['[1:' num2str(totalframes) ']'];
    prompt{i} = ['Choose frames for ' dataFile1{selectedFiles(i)}];
end

inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
interpVal = [];
for i = 1:length(selectedFiles)
    framesAll{i} = str2num(inputdialog{i});
end

% profile on;
for stack = selectedFiles
    
    fileIdx = find(selectedFiles == stack);
    %     fileInfo = fileInfoAll{fileIdx};
    frames = framesAll{fileIdx};
    frames = fliplr(frames);
    numFrames = frames;
    if isDcimg
        [imgHeight, imgWidth] = size(framedatatrans);
    else
        dataFileInfo = fileInfoAll{fileIdx};
        imgHeight = dataFileInfo.Height;
        imgWidth = dataFileInfo.Width;
    end
    if stack == selectedFiles(1)
        
        
        if ~exist(templateFile)
            templateFile2 = strsplit(templateFile,'\');
            templateFile = [calPath, templateFile2{end}];
            if ~exist(templateFile)
            [templateFile, templatePath] = uigetfile({'*.mat'},...
                'templateFile not found! Please reselect path:',calPath);
            templateFile = [templatePath templateFile];
            end
        end
        load(templateFile);
        
        templateSize = size(template,2);
        %         end
        clear templateFrames
        %load(calFile);
        
        if ~exist('templateFrames')
            goodFit_forward = logical(goodFit_f(calBeadIdx,:));
            % compute templates to use based upon fitted DG angles
        templateFrames = interp1(meanAngles(calBeadIdx,goodFit_forward),...
                1:length(meanAngles(calBeadIdx,goodFit_forward)),linspace(-80,80,7),'nearest'); %90:-30:-60
        end
        % if some angles are out of range, use the ends of the template
        % stack. first determine whether frames are increasing or
        % decreasing so the the 'ends' are chosen correctly.
        
        % use 2nd frame since 1st seems to give wrong angle (at least for
        % NHAs)
        if any(isnan(templateFrames))
            if nanmean(diff(templateFrames)) >= 0
                endFrames = [2 sum(goodFit_forward)];
            elseif nanmean(diff(templateFrames)) < 0
                endFrames = [sum(goodFit_forward) 2];
            else
                endFrames = nan;
            end
            if isnan(templateFrames(1))
                templateFrames(1) = endFrames(1);
            end
            if isnan(templateFrames(end))
                templateFrames(end) = endFrames(2);
            end
        end
        
        temp = inputdlg({['Input sets of frames corresponding to each template or the index of templates to use (ex. ' ...
            mat2str(templateFrames) ')']},...
            'Input template numbers',1, ...
            {mat2str(templateFrames)}); ...     
            templateFrames = str2num(temp{1});
        
    
        
        [darkFile, darkPath] = uigetfile({'*.dcimg'; '*.tif'},'Open image stack with dark counts (same parameters as SMACM data)',dataPath);
      
        
        %inputdlg(prompt,title,nl,def,options)
        temp = inputdlg({'Use (M)edian filtering, (W)avelet filtering or treat as (N)HA?',...
            'If median filtering, what sigma size in pix (0 for no smoothing)',...
            'If median filtering, what total window size in frames?',...
            'If median filtering, number of frames to "interpolate"',...
            'When fitting data, display raw (not background-subtracted) data?'},...
            'Input filtering options',1, ...
            {'M','15','101','10','1'});
        filterType = temp{1};
        if strcmp(filterType,'m')
            filterType = 'M';
        elseif strcmp(filterType,'w')
            filterType = 'W';
        elseif strcmp(filterType,'n')
            filterType = 'N';
        end
        
        medianBlurSigma = str2num(temp{2});
        windowSize = str2num(temp{3});
        interpVal = str2num(temp{4}); % number of frames to interpolate in median bg estimation. 1 uses each frame.
        dispRaw = logical(str2num(temp{5}));
        
        if strcmp(filterType,'N')
            nhaData = true;
            medianFilter = false;
        elseif strcmp(filterType,'M')
            medianFilter = true;
            nhaData = false;
        else
            medianFilter = false;
            nhaData = false;
        end
    end
   
    %% create output log filenames
    % saves in labeled directory if a channel is selected
    if channel == '0'
        if isDcimg
            outputFilePrefix{stack} = [dataFile1{1}(1:length(dataFile1{1})-6) filesep 'threshold ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        else
            outputFilePrefix{stack} = [dataFile{1}(1:length(dataFile{1})-11) filesep 'threshold ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        end
    else
        if isDcimg
            outputFilePrefix{stack} = [dataFile1{1}(1:length(dataFile1{1})-6) filesep channel(1) 'threshold ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        else
            outputFilePrefix{stack} = [dataFile{1}(1:length(dataFile{1})-11) filesep channel(1) 'threshold ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        end
    end
    iROI = [1,1];
    if channel == 'g'
        load('gainGreen.mat', 'gain');
    elseif channel == 'r'
        load('gainRed.mat', 'gain');
    else
        gain = ones(2048); %2048 is the size of the camera chip
    end

    gain = gain;
    if ~isequal(size(gain),[imgHeight imgWidth])
        warning('Gain and data image stack are not the same size. Please update initial coordinates.');
        
        dlg_title = 'Enter initial coordinates';
        prompt = { 'XO','YO'};
        def = { '1','1'};
        num_lines = 1;
        inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
        
        iROI(1) = str2double(inputdialog{1});
        iROI(2) = str2double(inputdialog{2});
        
        gain = gain(iROI(2):iROI(2)+imgHeight-1,iROI(1):iROI(1)+imgWidth-1);
    end
    
    
    mkdir(outputFilePrefix{stack});
    
    if stack== selectedFiles(1)
        %% Compute darkAvg counts
        
        if ~isequal(darkFile,0)
            if ~exist([darkPath darkFile])
                [darkFile, darkPath] = uigetfile({'*.Dcimg','*.tiff'},...
                    'darkFile not found! Please reselect path:');
            end
            
            if isDcimg
                darkFile1 = [darkPath darkFile];
            else
                darkFile1 = darkFile;
                darkFileInfo = imfinfo([darkPath darkFile]);
            end
            % Computes average of darkAvg frames for background subtraction
            if isDcimg
                dcimgfile_D = fullfile(darkPath, darkFile);
                
                dcimgfile_D = strrep(dcimgfile_D, '\', '\\');
                
                [framedata_D,totalframes_D]= dcimgmatlab(0, dcimgfile_D);
                framedatatrans_D = transpose (framedata_D);
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
                    [framedata_D,~]= dcimgmatlab(frame-1, dcimgfile_D);
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
        
        templateFrames = templateFrames';
        numTemplates = size(templateFrames,1);
        templateColors = jet(numTemplates);
        templateLocs = zeros(numTemplates,5);
        fitParam = zeros(1,8);
        [xIdx, yIdx] = meshgrid(1:templateSize,1:templateSize);
        
        
        hTemplate=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        for a=1:numTemplates
            
            % make minimum count level in template equal to 0
            template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                - min(min(template(templateFrames(a),:,:)));
            % normalize energy contained (sum of all counts) in the template
            template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                / sum(sum(template(templateFrames(a),:,:)));
            % finally, make mean of template equal to 0
            if ~nhaData
                template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                    - mean(mean(template(templateFrames(a),:,:)));
            end
            
            % find two largest peaks in template
            [tempY, tempX] = ind2sub([templateSize templateSize],find(imregionalmax(template(templateFrames(a),:,:))));
            temp = sortrows([tempX tempY template(sub2ind(size(template),templateFrames(a)*ones(length(tempX),1),tempY,tempX))],-3);
            
            % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
            fitParam(3) = temp(1,1);
            fitParam(4) = temp(1,2);
            fitParam(5) = temp(2,1);
            fitParam(6) = temp(2,2);
            fitParam(1) = temp(1,3);
            fitParam(2) = temp(2,3);
            fitParam(7) = mean(sigmaBounds);
            fitParam(8) = mean(sigmaBounds);
            
            
            lowerBound = [0 0 1 1 1 1 sigmaBounds(1) sigmaBounds(1)];
            upperBound = [max(max(template(templateFrames(a),:,:))) max(max(template(templateFrames(a),:,:))) ...
                templateSize templateSize templateSize templateSize ...
                sigmaBounds(2) sigmaBounds(2)];
            
            % Fit with lsqnonlin
            
            fitParam = lsqnonlin(@(x) ...
                f_doubleGaussianVector(x,squeeze(template(templateFrames(a),:,:)),0,xIdx,yIdx),...
                fitParam,lowerBound,upperBound,options);
            
            templateLocs(a,1:2) = fitParam(3:4);
            templateLocs(a,3:4) = fitParam(5:6);
            % calculate rough angle between peaks
            templateLocs(a,5) = 180/pi*atan2(templateLocs(a,2)-templateLocs(a,4), ...
                templateLocs(a,1)-templateLocs(a,3));
            
            subplot(1,numTemplates,a);
            imagesc(squeeze(template(templateFrames(a),:,:)));
            axis image;colormap hot;colorbar;
            hold on;
            plot(templateLocs(a,1),templateLocs(a,2),'.','MarkerEdgeColor', templateColors(a,:));
            plot(templateLocs(a,3),templateLocs(a,4),'.','MarkerEdgeColor', templateColors(a,:));
            title({['Frames ' mat2str(templateFrames(a,:))] ...
                ['Angle = ' num2str(templateLocs(a,5)) ' deg']});
        
        end
        imwrite(frame2im(getframe(hTemplate)),[outputFilePrefix{stack} 'templates.tif']);
        clear templateInfo tempX tempY temp xIdx yIdx;
        
        close(hTemplate); % closes template figure
        %% user picks ROI
        % pick region of interest by reading first frame and having user select
        % region
        
        % Compute average image
        avgImg = zeros(imgHeight,imgWidth);
        initialframestoavg = 200;
     
        avgImgFrames = min(initialframestoavg,length(frames));
        for a = 1:avgImgFrames
            if isDcimg
                dcimgfile = fullfile(dataPath, dataFile);
                dcimgfile = strrep(dcimgfile, '\', '\\');
                
                [framedata,totalframes]= dcimgmatlab(a-1, dcimgfile);
                totalframes = double(totalframes);
                framedatatrans = transpose (framedata);
                avgImg = avgImg + double(framedatatrans) - darkAvg;
            else
                dcimgfile = fullfile(dataPath, dataFile);
                frames1= fliplr(frames);
                avgImg = avgImg + double(imread([dataFile{stack}],frames1(a),'Info', fileInfo)) - darkAvg;
            end
        end
        avgImg = avgImg/avgImgFrames;
        
        hROI = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(avgImg,[max(0,min(avgImg(:))-5*std(avgImg(:))), max(0,min(avgImg(:)))+5*std(avgImg(:))]);axis image;colormap hot;
        if channel == 'g'
            ROI = imrect(gca,[1 1 270 270]);
        elseif channel == 'r'
            ROI = imrect(gca,[243 243 270 270]);
        else
            ROI = imrect(gca,[1 1 size(avgImg,1) size(avgImg,2)]);
        end
        
        title({'Shape box and double-click to choose region of interest for PSF extraction' ...
            ['[xmin ymin width height] = ' mat2str(ROI.getPosition)]...
            ['The displayed image is the average of the first ' num2str(avgImgFrames) ' frames']});
        addNewPositionCallback(ROI,@(p) title({'Shape box and double-click to choose region of interest for PSF extraction' ...
            ['[xmin ymin width height] = ' mat2str(p,3)]...
            'The displayed image is the average of the first 200 frames'}));
        % make sure rectangle stays within image bounds
        fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
        setPositionConstraintFcn(ROI,fcn);
        ROI = round(wait(ROI));
        % make sure ROI is an even number of pixels
        if mod(ROI(3),2)==1
            ROI(3) = ROI(3)-1;
        end
        if mod(ROI(4),2)==1
            ROI(4) = ROI(4)-1;
        end
        cropWidth = ROI(3);
        cropHeight = ROI(4);
        
        
        close(hROI) % closes ROI selection
        
        if usePolyROI
            hFOVmaskFig=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
            imagesc(avgImg(ROI(2):ROI(2)+ROI(4)-1, ...
                ROI(1):ROI(1)+ROI(3)-1),[0 max(0, min(avgImg(:)))+5*std(avgImg(:))]);
            axis image;colorbar
            colormap hot;
            title('Select ROIpoly of area to keep');
            [FOVmask, maskX, maskY] = roipoly;
            xCenter=(max(maskX)+min(maskX))/2;
            yCenter=(max(maskY)+min(maskY))/2;
            x1=(maskX-xCenter)*0.95+xCenter;
            y1=(maskY-yCenter)*0.95+yCenter;
            FOVmask1=roipoly(avgImg(ROI(2):ROI(2)+ROI(4)-1, ...
                ROI(1):ROI(1)+ROI(3)-1),x1,y1);
            close(hFOVmaskFig);
        end
        
 
        
        %% prepare template for template matching
        
        % pad template to same size as input
        templatePad = zeros(numTemplates,cropHeight,cropWidth);
        templateFT = zeros(numTemplates,cropHeight,cropWidth);
        for a=1:numTemplates
            
            
            templatePad(a,:,:) = padarray(squeeze(template(templateFrames(a),:,:)),...
                [(cropHeight-size(template,2))/2 ...
                (cropWidth-size(template,3))/2],min(min(template(templateFrames(a),:,:))));
            
            % multiplying by conjugate of template in FT domain is squivalent
            % to flipping the template in the real domain
            templateFT(a,:,:) = conj(fft2(squeeze(templatePad(a,:,:))));
        end
        clear templatePad temp;
        
        % apply Gaussian filter to phase correlation data to weight low frequencies
        % more heavily since SNR is higher there
        gaussianFilter = abs(fft2(fspecial('gaussian', [cropHeight cropWidth], gaussianFilterSigma)));
        
    end % end of the prep that is done only for first file.
    
    %% do template matching
    
    hMatchFig = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    totalPSFfits = zeros( 13000000, 6+15+3);
    numPSFfits = 0;
    startTime = tic;
    %     frameNum = 1;
    if size(frames,1) > 1 % make sure it will work in the for loop (need 1xn)
        frames = frames';
    end
    
    dataWindow = nan([2*floor(windowSize/2)+1,ROI(4),ROI(3)]);
    if medianBlurSigma~=0
        medBlurFilt = fspecial('gaussian',100,medianBlurSigma);
    else
        medBlurFilt=1;
    end
    lastFrame = nan;
    
    
    meanBG = nan(length(frames),1);
    meanSignal = nan(length(frames),1);
    
    PSFvals = cell(length(frames),1);
    meanBGvals = cell(length(frames),1);
    meanSignalvals = cell(length(frames),1);
    bkgndImgvals = cell(idivide(int32(length(frames)),int32(interpVal)),1);

    
  
    
    if medianFilter
        for c = frames(2:end)
            if length(frames)~= c
                currIdx = find(frames==c);               
                    if isnan(lastFrame)
                        if isDcimg
                            [bkgndImgMed, dataWindow] = f_medianFilter_DCImg(dcimgfile, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                        else
                            [bkgndImgMed, dataWindow] = f_medianFilter([dataFile{stack}], fileInfo, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                        end
                        lastFrame = c;
                       
                        bkgndImgvals{idivide(int32(c), int32(interpVal))+1}= bkgndImgMed;
              
                    elseif isempty(bkgndImgvals{idivide(int32(c), int32(interpVal))+1,1})
                        [bkgndImgMed, dataWindow] = f_medianFilter_DCImg(dcimgfile, darkAvg, ROI, frames, c, windowSize, dataWindow,lastFrame);
                        lastFrame =c;
                        bkgndImgvals{idivide(int32(c), int32(interpVal))+1}= bkgndImgMed;
                    end
            end
        end
    
       
        
        if isDcimg
            clear dataWindow mex
        else
            clear dataWindow
        end

    end

    
    
    dimMask = size(FOVmask);bufMask = 5;
    sMask = false(dimMask(1)+2*bufMask,dimMask(2)+2*bufMask,4);
    sMask(1+2*bufMask:dimMask(1)+2*bufMask,1+bufMask:dimMask(2)+bufMask,1) = FOVmask;
    sMask(1:dimMask(1),1+bufMask:dimMask(2)+bufMask,2) = FOVmask;
    sMask(1+bufMask:dimMask(1)+bufMask,1+2*bufMask:dimMask(2)+2*bufMask,3) = FOVmask;
    sMask(1+bufMask:dimMask(1)+bufMask,1:dimMask(2),4) = FOVmask;
    FTFOVmask = FOVmask & sMask(1+bufMask:dimMask(1)+bufMask,1+bufMask:dimMask(2)+bufMask,1) & sMask(1+bufMask:dimMask(1)+bufMask,1+bufMask:dimMask(2)+bufMask,2) & sMask(1+bufMask:dimMask(1)+bufMask,1+bufMask:dimMask(2)+bufMask,3) & sMask(1+bufMask:dimMask(1)+bufMask,1+bufMask:dimMask(2)+bufMask,4);
    
    forloopTime = 0;
    ROIinitial = ROI;
    %%
    gain = gain(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
    if stack == 1
        forloopST = tic;
        for c = 1:min(100,length(frames))
            
            if isDcimg
                [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dcimgfile,medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
            else
                [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dataFile{stack},medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
            end
            dataFT = fft2(data,cropHeight,cropWidth);
            maxPeakImg = zeros(cropHeight,cropWidth);
            % matrix PSFLocs stores information about double helices that were
            % found via template matching
            % rows are different matches
            % [xLocation yLocation matchingTemplateNumber matchConfidence];
            PSFLocs = zeros(100,4);
            numPSFLocs = 0;
            for b=1:numTemplates
                % try no prefiltering
                %H = 1;
                % try phase correlation
                %H = 1./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
                % try weighted phase correlation (emphasizing low frequency
                % components
                H = gaussianFilter./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
                if nhaData
                    H = gaussianFilter;
                end
                % normalize H so it doesn't add any energy to template match
                %H = H / sqrt(sum(abs(H(:)).^2));
                
                peakImg = ifftshift(ifft2(dataFT.*squeeze(templateFT(b,:,:)).*H));
                
                % normalize response of peakImg by dividing by number of pixels in
                % data
                %peakImg = peakImg / (cropHeight*cropWidth);
                maxPeakImg = max(maxPeakImg, peakImg);
                
                % only remember matches that are 5 standard deviations above the
                % mean
                peakThreshold = mean(peakImg(:))+5*std(peakImg(:));
                if nhaData
                    peakThreshold = mean(peakImg(:))+ std(peakImg(:));
                end
                peakImg(peakImg < peakThreshold) = peakThreshold;
                peakImg(~FTFOVmask) = peakThreshold;%This mask removes the excess brightness introduced by the fft/ifft of the FOV
                temp = find(imregionalmax(peakImg));
                % make sure threshold didn't eliminate all peaks and create
                % lots of matches
                if length(temp) < cropHeight*cropWidth/2;
                    [tempY, tempX] = ind2sub([cropHeight cropWidth],temp);
                    PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                        [tempX tempY b*ones(length(temp),1) peakImg(temp)];
                    numPSFLocs = numPSFLocs+length(temp);
                end
      
            end
            
            %% filter out extraneous matches due to very strong signals
            
            forloopTime = forloopTime + toc(forloopST);
            
            if numPSFLocs > 0
                % sort location matrix in decending order of confidence
                temp = sortrows(PSFLocs(1:numPSFLocs,:),-4);
                % copy most confident match to list of locations
                PSFLocs(1,:) = temp(1,:);
                numPSFLocs = 1;
                
                for b=2:size(temp,1)
                    % make sure that this candidate location is a minimum distance away
                    % from all other candidate locations
                    if sum((temp(b,1)-PSFLocs(1:numPSFLocs,1)).^2 + (temp(b,2)-PSFLocs(1:numPSFLocs,2)).^2 >= minDistBetweenSMs^2) == numPSFLocs
                        % add it to list of locations
                        numPSFLocs = numPSFLocs + 1;
                        PSFLocs(numPSFLocs,:) = temp(b,:);
                    end
                end
                if 2 > size(temp,1)
                    b = 1;
                end
            end
            
            %totalPSFfits(c, 1:numPSFLocs,1:6) = [b*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];
            
            PSFvals{c,1} = [frames(c)*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];
            %[totalPSFfits, numPSFfits]  = fixPSFfits(PSFLocs, numPSFLocs, totalPSFfits, b, numPSFfits); %#ok<PFTUS,PFTIN>
            
            %% output an example image for each threshold level
            %  so that user can pick appropriate threshold later
            %             matchInfo = zeros(numPSFLocs,3);
            for b=1:numPSFLocs
                moleThreshold = round(PSFLocs(b,4)*100000);
                if nhaData
                    moleThreshold=round(PSFLocs(b,4)); %makes number more reasonable
                end
                moleFileName = ['template ' num2str(PSFLocs(b,3)) ' threshold ' num2str(moleThreshold,'%g') '.jpeg']; %%f6.4 without scaling
                if isempty(dir([outputFilePrefix{stack} moleFileName]))
                    % create indices to isolate image of candidate molecule
                    [xIdx, yIdx] = meshgrid(PSFLocs(b,1)-9:PSFLocs(b,1)+9, ...
                        PSFLocs(b,2)-9:PSFLocs(b,2)+9);
                    % make sure indices are inside ROI
                    if min(xIdx(:)) < 1
                        xIdx = xIdx + (1-min(xIdx(:)));
                    end
                    if max(xIdx(:)) > cropWidth
                        xIdx = xIdx - (max(xIdx(:))-cropWidth);
                    end
                    if min(yIdx(:)) < 1
                        yIdx = yIdx + (1-min(yIdx(:)));
                    end
                    if max(yIdx(:)) > cropHeight
                        yIdx = yIdx - (max(yIdx(:))-cropHeight);
                    end
                    
                    %SSIM test
                    templateNum = PSFLocs(b,3);
                    templateImg = squeeze(template(templateFrames(templateNum),:,:));
                    templateCenterX = round((templateLocs(PSFLocs(b,3),1) + templateLocs(PSFLocs(b,3),3))/2);
                    if templateCenterX < 10
                        templateCenterX = 10;
                    elseif templateCenterX > 11
                        templateCenterX = 11;
                    end
                    templateCenterY = round((templateLocs(PSFLocs(b,3),2) + templateLocs(PSFLocs(b,3),4))/2);
                     if templateCenterY < 10
                        templateCenterY = 10;
                    elseif templateCenterY > 11
                        templateCenterY = 11;
                    end
                    end
                    templateImgResize = templateImg(templateCenterX-9:templateCenterX+9,templateCenterY-9:templateCenterY+9);
                    templateImgResize = 1+round(255*(templateImgResize-min(templateImgResize(:)))/(max(templateImgResize(:))-min(templateImgResize(:))));
                    
                    % output a picture of the
                    img = data(yIdx(:,1),xIdx(1,:));
                    img = 1+round(255*(img-min(img(:)))/(max(img(:))-min(img(:))));
                    img = double(img);
                    match = ssim(img,templateImgResize);
             
                    matchInfo{c}(b,1) = PSFLocs(b,3);
                    matchInfo{c}(b,2) = moleThreshold;
                    matchInfo{c}(b,3) = match;
                    if match > ssimThresh
                        moleFileName = ['template ' num2str(PSFLocs(b,3)) ' threshold ' num2str(moleThreshold,'%g') ' SSIM' num2str(match) '.jpeg']; %%f6.4 without scaling
                        imwrite(imresize(ind2rgb(img,hot(256)),3,'nearest'),[outputFilePrefix{stack} moleFileName]);
                    end
                    %end SSIM test
                   end
%             end
            
            
            numPSFfits = numPSFfits+numPSFLocs;
            
            
            %%  plot results of template matching and fitting
            set(0,'CurrentFigure',hMatchFig);
            subplot('Position',[0.025 0.025 .9/2 .95],'parent',hMatchFig);
            imagesc(maxPeakImg,[0 3*peakThreshold]);axis image;
            title({'Peaks correspond to likely template matches' ...
                [num2str(numPSFLocs) ' matches found']});
            
            subplot('Position',[0.525 0.025 .9/2 .95],'parent',hMatchFig);
            imagesc(data); axis image;colormap hot;
            hold on;
            for b=1:numPSFLocs
                plot(PSFLocs(b,1), PSFLocs(b,2), 'o', ...
                    'MarkerSize', 15*PSFLocs(b,4)/peakThreshold, ...
                    'MarkerEdgeColor', templateColors(PSFLocs(b,3),:));
            end
            hold off;
            title({['Frame ' num2str(frames(c)) ': raw data - bkgnd & dark offset'] ...
                ['ROI [xmin ymin width height] = ' mat2str(ROI)]});
            
            drawnow;
            
            meanSignalvals{c,1} = [frames(c), mean(data(FOVmask))];
            meanBGvals{c,1} = [frames(c), mean(bkgndImg(FOVmask))];
   
            
        end
        
        forloopTime = forloopTime
        
        parforloopTime = 0;
        parforLoopST = tic;
        fileInfo='';
        if (floor(length(frames)/100)*100) > 100
            for n = 101:100:length(frames)
                p = n + 99;
                p = min(p,length(frames));
                parfor c = n:p
                    if isDcimg
                        [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dcimgfile,medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
                    else
                        [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dataFile{stack},medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
                    end
                    dataFT = fft2(data,cropHeight,cropWidth);
                    maxPeakImg = zeros(cropHeight,cropWidth);
                    % matrix PSFLocs stores information about double helices that were
                    % found via template matching
                    % rows are different matches
                    % [xLocation yLocation matchingTemplateNumber matchConfidence];
                    PSFLocs = zeros(100,4);
                    numPSFLocs = 0;
                    for b=1:numTemplates
                        % try no prefiltering
                        %H = 1;
                        % try phase correlation
                        %H = 1./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
                        % try weighted phase correlation (emphasizing low frequency
                        % components
                        H = gaussianFilter./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
                        if nhaData
                            H = gaussianFilter;
                        end
                        % normalize H so it doesn't add any energy to template match
                        %H = H / sqrt(sum(abs(H(:)).^2));
                        
                        peakImg = ifftshift(ifft2(dataFT.*squeeze(templateFT(b,:,:)).*H));
                        
                        % normalize response of peakImg by dividing by number of pixels in
                        % data
                        %peakImg = peakImg / (cropHeight*cropWidth);
                        maxPeakImg = max(maxPeakImg, peakImg);
                        
                        % only remember matches that are 5 standard deviations above the
                        % mean
                        peakThreshold = mean(peakImg(:))+5*std(peakImg(:));
                        if nhaData
                            peakThreshold = mean(peakImg(:))+ std(peakImg(:));
                        end
                        peakImg(peakImg < peakThreshold) = peakThreshold;
                        temp = find(imregionalmax(peakImg));
                        % make sure threshold didn't eliminate all peaks and create
                        % lots of matches
                        if length(temp) < cropHeight*cropWidth/2;
                            [tempY, tempX] = ind2sub([cropHeight cropWidth],temp);
                            PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                                [tempX tempY b*ones(length(temp),1) peakImg(temp)];
                            numPSFLocs = numPSFLocs+length(temp);
                        end
                    end

                    %% filter out extraneous matches due to very strong signals
                    
                    
                    
                    if numPSFLocs > 0
                        % sort location matrix in decending order of confidence
                        temp = sortrows(PSFLocs(1:numPSFLocs,:),-4);
                        % copy most confident match to list of locations
                        PSFLocs(1,:) = temp(1,:);
                        numPSFLocs = 1;
                        for b=2:size(temp,1)
                            % make sure that this candidate location is a minimum distance away
                            % from all other candidate locations
                            if sum((temp(b,1)-PSFLocs(1:numPSFLocs,1)).^2 + (temp(b,2)-PSFLocs(1:numPSFLocs,2)).^2 >= minDistBetweenSMs^2) == numPSFLocs
                                % add it to list of locations
                                numPSFLocs = numPSFLocs + 1;
                                PSFLocs(numPSFLocs,:) = temp(b,:);
                            end
                        end
                        if 2 > size(temp,1)
                            b = 1;
                        end
                    end
                    
                    %totalPSFfits(c, 1:numPSFLocs,1:6) = [b*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];
                    PSFvals{c,1} = [frames(c)*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];
                    %[totalPSFfits, numPSFfits]  = fixPSFfits(PSFLocs, numPSFLocs, totalPSFfits, b, numPSFfits); %#ok<PFTUS,PFTIN>
                    
                    %% output an example image for each threshold level
                    %  so that user can pick appropriate threshold later
                    
                    for b=1:numPSFLocs
                        moleThreshold = round(PSFLocs(b,4)*100000);
                        if nhaData
                            moleThreshold=round(PSFLocs(b,4)); %makes number more reasonable
                        end
                        moleFileName = ['template ' num2str(PSFLocs(b,3)) ' threshold ' num2str(moleThreshold,'%g') '.jpeg']; %%f6.4 without scaling
                        if isempty(dir([outputFilePrefix{stack} moleFileName]))
                            % create indices to isolate image of candidate molecule
                            [xIdx, yIdx] = meshgrid(PSFLocs(b,1)-9:PSFLocs(b,1)+9, ...
                                PSFLocs(b,2)-9:PSFLocs(b,2)+9);
                            % make sure indices are inside ROI
                            if min(xIdx(:)) < 1
                                xIdx = xIdx + (1-min(xIdx(:)));
                            end
                            if max(xIdx(:)) > cropWidth
                                xIdx = xIdx - (max(xIdx(:))-cropWidth);
                            end
                            if min(yIdx(:)) < 1
                                yIdx = yIdx + (1-min(yIdx(:)));
                            end
                            if max(yIdx(:)) > cropHeight
                                yIdx = yIdx - (max(yIdx(:))-cropHeight);
                            end
                            
                            %SSIM test
                            templateNum = PSFLocs(b,3);
                            templateImg = squeeze(template(templateFrames(templateNum),:,:));
                            templateCenterX = round((templateLocs(PSFLocs(b,3),1) + templateLocs(PSFLocs(b,3),3))/2);
                            if templateCenterX < 10
                                templateCenterX = 10;
                            elseif templateCenterX > 11
                                templateCenterX = 11;
                            end
                            templateCenterY = round((templateLocs(PSFLocs(b,3),2) + templateLocs(PSFLocs(b,3),4))/2);
                            if templateCenterY < 10
                                templateCenterY = 10;
                            elseif templateCenterY > 11
                                templateCenterY = 11;
                            end
                                templateCenterY = round(size(templateImg,2)/2);
                        
                            
                            templateImgResize = templateImg(templateCenterX-9:templateCenterX+9,templateCenterY-9:templateCenterY+9);
                            templateImgResize = 1+round(255*(templateImgResize-min(templateImgResize(:)))/(max(templateImgResize(:))-min(templateImgResize(:))));
                            
                            
                            % output a picture of the example
                            img = data(yIdx(:,1),xIdx(1,:));
                            img = 1+round(255*(img-min(img(:)))/(max(img(:))-min(img(:))));
                            img = double(img);
                            match = ssim(img,templateImgResize);
                            
                            matchInfo{c}(b,1) = PSFLocs(b,3);
                            matchInfo{c}(b,2) = moleThreshold;
                            matchInfo{c}(b,3) = match;
                            if match > ssimThresh
                                moleFileName = ['template ' num2str(PSFLocs(b,3)) ' threshold ' num2str(moleThreshold,'%g') ' SSIM' num2str(match) '.jpeg']; %%f6.4 without scaling
                                imwrite(imresize(ind2rgb(img,hot(256)),3,'nearest'),[outputFilePrefix{stack} moleFileName]);
                            end
                        end
                  
                            %end SSIM test
                        end
                    end
                    
                    numPSFfits = numPSFfits+numPSFLocs;
                    meanSignalvals{c,1} = [frames(c), mean(data(FOVmask))];
                    meanBGvals{c,1} = [frames(c), mean(bkgndImg(FOVmask))];
     
                end
                if isDcimg
                    clear mex;
                end
                psays = ['Frames processed: ' num2str(p)];
                disp(psays);
        
        
        
    else
%         fileInfo=imfinfo(dataFile{stack});
        for c =1:length(frames)
            if isDcimg
                [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dcimgfile,medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
            else
                [data,bkgndImg] = FrameData(c,isDcimg,gain,ROI,frames,darkAvg,FOVmask,FOVmask1,dataFile{stack},medianFilter,nhaData,usePolyROI,bkgndImgvals,interpVal,medianBlurSigma,medBlurFilt);
            end
            dataFT = fft2(data,cropHeight,cropWidth);
            maxPeakImg = zeros(cropHeight,cropWidth);
            % matrix PSFLocs stores information about double helices that were
            % found via template matching
            % rows are different matches
            % [xLocation yLocation matchingTemplateNumber matchConfidence];
            PSFLocs = zeros(100,4);
            numPSFLocs = 0;
            for b=1:numTemplates
                % try no prefiltering
                %H = 1;
                % try phase correlation
                %H = 1./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
                % try weighted phase correlation (emphasizing low frequency
                % components
                H = gaussianFilter./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
                if nhaData
                    H = gaussianFilter;
                end
                % normalize H so it doesn't add any energy to template match
                %H = H / sqrt(sum(abs(H(:)).^2));
                
                peakImg = ifftshift(ifft2(dataFT.*squeeze(templateFT(b,:,:)).*H));
                
                % normalize response of peakImg by dividing by number of pixels in
                % data
                %peakImg = peakImg / (cropHeight*cropWidth);
                maxPeakImg = max(maxPeakImg, peakImg);
                
                % only remember matches that are 5 standard deviations above the
                % mean
                peakThreshold = mean(peakImg(:))+5*std(peakImg(:)); 
                if nhaData
                    peakThreshold = mean(peakImg(:))+ std(peakImg(:));
                end
                peakImg(peakImg < peakThreshold) = peakThreshold;
                temp = find(imregionalmax(peakImg));
                % make sure threshold didn't eliminate all peaks and create
                % lots of matches
                if length(temp) < cropHeight*cropWidth/2
                    [tempY, tempX] = ind2sub([cropHeight cropWidth],temp);
                    PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                        [tempX tempY b*ones(length(temp),1) peakImg(temp)];
                    numPSFLocs = numPSFLocs+length(temp);
                end
            end
            %% filter out extraneous matches due to very strong signals
            if numPSFLocs > 0
                % sort location matrix in decending order of confidence
                temp = sortrows(PSFLocs(1:numPSFLocs,:),-4);
                % copy most confident match to list of locations
                PSFLocs(1,:) = temp(1,:);
                numPSFLocs = 1;
                for b=2:size(temp,1)
                    % make sure that this candidate location is a minimum distance away
                    % from all other candidate locations
                    if sum((temp(b,1)-PSFLocs(1:numPSFLocs,1)).^2 + (temp(b,2)-PSFLocs(1:numPSFLocs,2)).^2 >= minDistBetweenSMs^2) == numPSFLocs
                        % add it to list of locations
                        numPSFLocs = numPSFLocs + 1;
                        PSFLocs(numPSFLocs,:) = temp(b,:);
                    end
                end
                if 2 > size(temp,1)
                    b = 1;
                end
            end
            PSFvals{c,1} = [frames(c)*ones(numPSFLocs,1) (1:numPSFLocs)' PSFLocs(1:numPSFLocs,:)];
            %% output an example image for each threshold level
            %  so that user can pick appropriate threshold later
            
            for b=1:numPSFLocs
                moleThreshold = round(PSFLocs(b,4)*100000);
                if nhaData
                    moleThreshold=round(PSFLocs(b,4)); %makes number more reasonable
                end
                moleFileName = ['template ' num2str(PSFLocs(b,3)) ' threshold ' num2str(moleThreshold,'%g') '.jpeg']; %%f6.4 without scaling
                if isempty(dir([outputFilePrefix{stack} moleFileName]))
                    % create indices to isolate image of candidate molecule
                    [xIdx, yIdx] = meshgrid(PSFLocs(b,1)-boxRadius:PSFLocs(b,1)+boxRadius, ...
                        PSFLocs(b,2)-boxRadius:PSFLocs(b,2)+boxRadius);
                    % make sure indices are inside ROI
                    if min(xIdx(:)) < 1
                        xIdx = xIdx + (1-min(xIdx(:)));
                    end
                    if max(xIdx(:)) > cropWidth
                        xIdx = xIdx - (max(xIdx(:))-cropWidth);
                    end
                    if min(yIdx(:)) < 1
                        yIdx = yIdx + (1-min(yIdx(:)));
                    end
                    if max(yIdx(:)) > cropHeight
                        yIdx = yIdx - (max(yIdx(:))-cropHeight);
                    end
                    
                    % output a picture of the
                    img = data(yIdx(:,1),xIdx(1,:));
                    img = 1+round(255*(img-min(img(:)))/(max(img(:))-min(img(:))));
                    imwrite(imresize(ind2rgb(img,hot(256)),3,'nearest'),[outputFilePrefix{stack} moleFileName]);
                end
            end
            
            
            numPSFfits = numPSFfits+numPSFLocs;
            
            
            
    
            
            if medianFilter
            meanSignalvals{c} = {currIdx, mean(data(FOVmask))};
            meanBGvals{c} = {currIdx, mean(bkgndImg(FOVmask))};
            else
                meanSignalvals{c} = {mean(data(FOVmask))};
                meanBGvals{c} = {mean(bkgndImg(FOVmask))};
            end
    
            
        end
        parforloopTime= parforloopTime;
        end
    end
    
    elapsedTime = toc(startTime);
    
    %SSIM
    if useSSIM == 1
        matchOutput = [];
        matchOutput = vertcat(matchOutput,matchInfo{:});
        matchOutput = matchOutput(matchOutput(:,3) > ssimThresh,:);
        
        threshold = 0.090;
        matchAvg = [];
        tempThresh1 = [];
        for a = 1:numTemplates
            temp = unique(matchOutput((matchOutput(:,1)==a),2));
            matchTemp = [];
            if ~isempty(temp)
                for b = 1:length(temp)
                    matchTemp(b,1) = a;
                    matchTemp(b,2) = temp(b);
                    matchTemp(b,3) = mean(matchOutput(matchOutput(:,1)==a & matchOutput(:,2)==temp(b),3));
                end
                matchAvg = vertcat(matchAvg,matchTemp);
            end
        end
        matchAvg = sortrows(matchAvg,[1 2]);
        tempThresh1 = zeros(numTemplates,2);
        for g = 1:numTemplates
            if ~isempty(matchAvg(matchAvg(:,1)==g,:))
                P = fit(matchAvg(matchAvg(:,1)==g,2),matchAvg(matchAvg(:,1)==g,3),'smoothingspline');
                fig1 = figure; plot(P,matchAvg(matchAvg(:,1)==g,2),matchAvg(matchAvg(:,1)==g,3));
                P1 = feval(P,matchAvg(matchAvg(:,1)==g,2));
                [minVal minIdx] = min(abs((P1(:,1) - threshold)));
                temp = matchAvg(matchAvg(:,1)==g,:);
                tempThresh1(g,1) = g;
                tempThresh1(g,2) = temp(minIdx,2);
                title({['Template ' num2str(g)],['Threshold: ' num2str(tempThresh1(g,2))]});
                filename = ['Template ' num2str(g) ' SSIM'];
                saveas(fig1, [outputFilePrefix{stack} filename ]);
            end
        end
        
        tempThresh = tempThresh1(:,2)';
    else
        tempThresh = ones(1,numTemplates);
    end

    
    clear data bkgnd residual fileInfo maxPeakImg reconstructImg xIdx yIdx temp;
    close(hMatchFig) % closes fitting figure to prevent messiness
    
    fps = length(frames)/elapsedTime
    moleculesPerSec = numPSFfits/elapsedTime
    
    %% output data to external files
    
    textHeader = {'frame number' 'molecule number' 'template x location in ROI (px)' 'template y location in ROI (px)' ...
        'matching template number' 'match confidence (au)' ...
        'amp 1 (counts)' 'amp 2 (counts)' 'x location 1 (px)' 'y location 1 (px)' ...
        'x location 2 (px)' 'y location 2 (px)' 'sigma 1 (px)' 'sigma 2 (px)' 'background mean (counts)' ...
        'total fit error (counts)' 'good fit flag' 'x center (nm)' 'y center (nm)' ...
        'angle (deg)' 'number of photons' 'aberration corrected x location (nm)' ...
        'aberration corrected y location (nm)' 'z location (nm)'};
    % save fit info to MATLAB mat file
    save('-v7.3', [outputFilePrefix{stack} 'threshold output.mat']);
    
    
end
end