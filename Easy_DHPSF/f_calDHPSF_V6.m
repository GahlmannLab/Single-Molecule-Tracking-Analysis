
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

function [outputFilePrefix, numBeads] = ...
    f_calDHPSF_V6(conversionGain,nmPerPixel,boxRadius,channel,sigmaBounds,lobeDistBounds, fittingMethod)
% f_calDHPSF is a module in easy_dhpsf that calibrates the z vs. angle
% response of the DH-PSF using one or more discrete fluorescent particles,
% usually fluorescent beads. To generate the data, an objective stepper
% should be used, and the order and size of these steps read out
% from a .dat log file. In addition to the calibration, f_calDHPSF
% generates a series of templates used for template matching.
%
% Instrument Specific Parameters
% profile on;
dlg_title = 'Set EM Gain, background subtraction type';
prompt = { 'EM Gain (1 if no gain):',...
    'Use wavelet subtraction?',...
    'Many NHA holes?',...
    'If so, are they at 45 degrees?',...
    'Is DHPSF horizontal when in focus?',...
    'Do you want to fine select bead positions?', ...
    'Do you want to reprocess a run?'};
def = { '1','1','0','0','0','1', '0'};
fileInfo = '';
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);

EMGain = str2double(inputdialog{1});
if EMGain < 1 || isnan(EMGain)
    warning('EMGain should be a number >= 1. Setting to 1...');
    EMGain = 1;
end

% whether to use wavelet background subtraction, or subtract mean of image
% Off by default to minimize any possible change to shape of
% templates.
useWaveSub = logical(str2num(inputdialog{2})); %#ok<*ST2NM>
fillNHA = logical(str2num(inputdialog{3}));
angled = logical(str2num(inputdialog{4}));
horizontal = logical(str2num(inputdialog{5}));
fineSelect = logical(str2num(inputdialog{6}));
Restart = logical(str2num(inputdialog{7}));
conversionFactor = conversionGain/EMGain;
ampRatioLimit = 0.5;
sigmaRatioLimit = 0.4;
blurSize = 0.5*160/nmPerPixel;
scrsz = get(0,'ScreenSize');

if ~Restart
    
    %% ask user for relevant datafiles
    [dataFile, dataPath] = uigetfile({'*.dcimg;*.tif;*.tiff;*.mat','Standard image types';'*.*','All Files'},'MultiSelect','on','Open image stack for data processing');
    if isequal(dataFile,0)
        error('User cancelled the program');
    end
    
    if isstr(dataFile)
        isTif = strsplit(dataFile, '.');
        strIsTif = isTif{2};
        isDcimg =strcmpi(strIsTif, 'dcimg');
        ismat =strcmpi(strIsTif, 'mat');
    else
        isTif = strsplit(dataFile{1}, '.');
        strIsTif = isTif{2};
        isDcimg =strcmpi(strIsTif, 'dcimg');
        ismat =strcmpi(strIsTif, 'mat');       
     
    end
    

    
    if (isDcimg)
        if ischar(dataFile)==1
            dataFile1 = [dataPath dataFile];
            dataFile1 = {dataFile1};
        end
    %This allows user to load mat files also
    elseif ismat
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
    % for i = 1:length(dataFile1)
    %     dataFile1{i}= [dataPath dataFile1{i}];
    % end
    [darkFile, darkPath] = uigetfile({'*.dcimg;*.tif;*.tiff','Standard image types';'*.*','All Files'},'Open image stack with dark counts (same parameters as calibration)',dataPath);
    [logFile, logPath] = uigetfile({'*.dat;*.txt','Text Files';'*.*','All Files'},'Open sequence log file for calibration',dataPath);
    logFile = [logPath logFile];
    if isequal(logFile,0)
        error('A sequence log file must be specified for the DHPSF calibration');
    end
    
    if isDcimg
        lastDir=dataPath;
        dcimgfile = fullfile(dataPath, dataFile);
        dcimgfile = strrep(dcimgfile, '\', '\\');
        [framedata,totalframes]= dcimgmatlab(0, dcimgfile);
        totalframes = double(totalframes);
        [imgHeight, imgWidth] = size(transpose(framedata));
        numFrames = double(totalframes);
        numFramesTotal = numFrames;
        numFiles =1;
    %Modified 07/27/17--Alma   
    elseif ismat
        fov=size((framedata_after{1,1}));
        frames_perloop=fov(3);
        loop_info=size(framedata_after);
        numloops=loop_info(2);
        totalframes=frames_perloop*numloops;
        imgHeight=fov(2);
        imgWidth=fov(1);
    else
        numFiles = length(dataFile);
        numFramesTotal = 0;
        for n=1:numFiles
            dataFileInfo = imfinfo(dataFile{n});
            % dataFileInfo = [dataFileInfo; dataFileInfo1];
            numFrames = length(dataFileInfo);
            numFramesTotal = numFramesTotal + numFrames;
            totalframes = double(numFramesTotal);
        end
        imgHeight = dataFileInfo.Height;
        imgWidth = dataFileInfo.Width;
    end
    % saves in labeled directory if a channel is selected
    if channel == '0'
        if isDcimg
            outputFilePrefix = [dataFile1{1}(1:length(dataFile1{1})-6) filesep 'calibration ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        else
            outputFilePrefix = [dataFile{1}(1:length(dataFile{1})-11) filesep 'calibration ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        end
    else
        if isDcimg
            outputFilePrefix = [dataFile1{1}(1:length(dataFile1{1})-6) filesep channel(1) ' calibration ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        else
            outputFilePrefix = [dataFile{1}(1:length(dataFile{1})-11) filesep channel(1) ' calibration ' ...
                datestr(now,'yyyymmdd HHMM') filesep];
        end
    end
    % make new directory in the current directory
    mkdir(outputFilePrefix);
    
    %% Analyze the C-focus scan series
    
    sifLogData =  importdata(logFile);
    sifLogData1 = sifLogData(~(sifLogData(:,1)==-1),:);   % transition entries "-1" removed
    sifLogData = [sifLogData(2:end,:); sifLogData(1,:)];  % move the first line with minus values to the end
    %  n=1;
    step_indices = find(sifLogData(:,1)==-1);
    numSteps = length(step_indices);
    if diff(diff(step_indices))==0
        frame_step=diff(step_indices);
    else
        disp ('Change the condition for finding step_indices')
    end
    meanCFocPos=zeros(1,numSteps);
    stdCFocPos = meanCFocPos;
    meanCXPos = meanCFocPos;
    stdCXPos = meanCFocPos;
    meanCYPos = meanCFocPos;
    stdCYPos = meanCFocPos;
    for i=1:length(step_indices)
        meanCFocPos(i) = mean(sifLogData(step_indices(i)-frame_step(1)+1:step_indices(i)-1,4))*1000; % in units of nm
        stdCFocPos(i) = std(sifLogData(step_indices(i)-frame_step(1)+1:step_indices(i)-1,4))*1000;   % in units of nm
        if size(sifLogData,2) > 4
            meanCXPos(i) = mean(sifLogData(step_indices(i)-frame_step(1)+1:step_indices(i)-1,5))*1000; % in units of nm
            stdCXPos(i) = std(sifLogData(step_indices(i)-frame_step(1)+1:step_indices(i)-1,5))*1000;   % in units
            meanCYPos(i) = mean(sifLogData(step_indices(i)-frame_step(1)+1:step_indices(i)-1,6))*1000; % in units of nm
            stdCYPos(i) = std(sifLogData(step_indices(i)-frame_step(1)+1:step_indices(i)-1,6))*1000;   % in units
        end
        startFrame(i)= step_indices(i)- frame_step(1)+2-i;
        endFrame(i) = step_indices(i)-i;
    end
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
    variance = ReadN_DI.^2; % changed by Ting Yan, 12/27/2016
    
    if ~isequal(size(gain),[imgHeight imgWidth])
        warning('Gain and data image stack are not the same size. Please update initial coordinates.');
        
        dlg_title = 'Enter initial coordinates';
        prompt = { 'XO','YO'};
        def = { '1','1'};
        num_lines = 1;
        inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
        
        iROI(1) = str2double(inputdialog{1});
        iROI(2) = str2double(inputdialog{2});
        
        %     if or(size(gain)< [iROI(1)+ imgWidth-1,iROI(2)+ imgHeight-1])
        %         error('Out of bounds!')
        %     end
        gain = gain(iROI(2):iROI(2)+imgHeight-1,iROI(1):iROI(1)+imgWidth-1);
        variance = variance(iROI(2):iROI(2)+imgHeight-1,iROI(1):iROI(1)+imgWidth-1);
    end
    
    %% Compute dark counts
    if ~isequal(darkFile,0)
        %darkFile = [darkPath darkFile];
        % Computes average of dark frames for background subtraction
        if isDcimg
            dcimgfile_D = fullfile(darkPath, darkFile);
            dcimgfile_D = strrep(dcimgfile_D, '\', '\\');
            [framedata_dark,numDarkFrames]= dcimgmatlab(0, dcimgfile_D);
            [darkHeight, darkWidth] = size(transpose(framedata_dark));
            darkAvg = zeros(darkHeight, darkWidth);
            for frame = 0:numDarkFrames-1
                [framedata_dark,numDarkFrames]= dcimgmatlab(frame, dcimgfile_D);
                framedatatrans = transpose (framedata_dark);
                darkAvg = darkAvg + double(framedatatrans);
            end
        else
            darkFile = [darkPath darkFile];
            % Computes average of dark frames for background subtraction
            darkFileInfo = imfinfo(darkFile);
            numDarkFrames = length(darkFileInfo);
            darkAvg = zeros(darkFileInfo(1).Height,darkFileInfo(1).Width);
            for frame = 1:numDarkFrames
                darkAvg = darkAvg + double(imread(darkFile,frame,'Info',darkFileInfo));
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
    clear darkFileInfo;
    
    %% Pick region of interest for analysis
    
    % Plots the avg image .tif
    dataAvg = zeros(imgHeight,imgWidth);
    if ~fillNHA
        for n=1:numFiles
            if isDcimg
                for frame = 0:min(190,totalframes)-1;
                    [framedata,totalframes]= dcimgmatlab(frame, dcimgfile);
                    totalframes = double(totalframes);
                    framedatatrans = transpose(framedata);
                    dataAvg = dataAvg + double(framedatatrans);
                end
            else
                for frame = 1:min(190,length(dataFileInfo));
                    dataAvg = dataAvg + double(imread(dataFile{n},frame,'Info',dataFileInfo));
                end
            end
        end
        dataAvg = dataAvg/190 - darkAvg;
    else
        for n=1:numFiles
            if isDcimg
                for frame = 1:10:totalframes
                    [framedata,~]= dcimgmatlab(frame-1, dcimgfile);
                    framedatatrans = transpose(framedata);
                    dataAvg = dataAvg + double(framedatatrans);
                end
            else
                totalframes = length(dataFileInfo);
                for frame = 1:10:totalframes;
                    dataAvg = dataAvg + double(imread(dataFile{n},frame,'Info',dataFileInfo));
                end
            end
        end
        dataAvg = dataAvg/(totalframes/10) - darkAvg;
        dataAvg(dataAvg< median(dataAvg(:))) = median(dataAvg(:));
    end
    
    if ~fillNHA
        % could change contrast here
        hROI=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        % imagesc(dataAvg, [0 max(max(dataAvg))/5]);colormap parula;
        imagesc(dataAvg, [min(min(dataAvg)) max(max(dataAvg))]); axis image;colormap jet;colorbar;
        ROI = imrect(gca,[imgWidth/4 imgHeight/4 imgWidth/2 imgHeight/2]);
        
        title({'Shape box and double-click to choose region of interest for PSF fitting' ...
            mat2str(ROI.getPosition)});
        addNewPositionCallback(ROI,@(p) title({'Double-click to choose region of interest for PSF fitting' ...
            ['[xmin ymin width height] = ' mat2str(p,3)]}));
        % make sure rectangle stays within image bounds
        fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
        setPositionConstraintFcn(ROI,fcn);
        ROI = round(wait(ROI));
        cropHeight = ROI(4);
        cropWidth = ROI(3);
        close(hROI);
    else
        ROI = [1,1,size(dataAvg)];
        cropHeight = ROI(4);
        cropWidth = ROI(3);
    end
    
    %% Ask user for bead location(s)
    % could change contrast here
    hLocs=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    % imagesc(dataAvg(ROI(2):ROI(2)+ROI(4)-1, ...
    %     ROI(1):ROI(1)+ROI(3)-1),[0 max(max(dataAvg))/10]); colormap parula;
    imagesc(dataAvg(ROI(2):ROI(2)+ROI(4)-1, ...
        ROI(1):ROI(1)+ROI(3)-1),[min(min(dataAvg)) max(max(dataAvg))]);axis image;colorbar;colormap jet;
    if ~fillNHA
        title('Use LMB to select fiducials. Hit enter to stop or use RMB to select the final fiducial.');
    end
    % if fillNHA && ~fourQuads
    %     title('Select four fiducials that form a large box. Then select fiducials to outline the edge of your "FOV".');
    % elseif fillNHA && fourQuads
    %     %     title('Select four fiducials to form a small box. Then pick edge of fiducials on the middle inside corner.');
    % end
    hold on;
    % User will click on beads in image and then text will be imposed on the
    % image corresponding to the bead number.
    % moleLocs is a n by 2 array -- x y values for each bead
    moleLocs = [];
    erase = ones(0,2);
    n = 0;
    % Loop, collecting bead locations and drawing text to mark them
    but = 1;
    if ~fillNHA
        while but == 1 % right click
            [xi,yi,but] = ginput(1);
            if isempty(xi)
                break
            end
            if but == 2; % left click
                erase = [erase; round([xi yi])];
                text(xi,yi,'*','color','red','fontsize',13,'fontweight','bold');
                but = 1;
                continue
            end
            n = n+1;
            text(xi,yi,num2str(n),'color','white','fontsize',13,'fontweight','bold');
            moleLocs(n,:) = round([xi yi]);
        end
        hold off;
        if fineSelect           % Select the positions for the beads with higher precision
            ROIsize = 20;
            hLocs_1=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
            for n = 1:size(moleLocs,1)
                imagesc(dataAvg(ROI(2)+moleLocs(n,2)-ROIsize-1:ROI(2)+moleLocs(n,2)+ROIsize-1,...
                    ROI(1)+moleLocs(n,1)-ROIsize-1:ROI(1)+moleLocs(n,1)+ROIsize-1),...
                    [0 max(max(dataAvg))/10]);axis image;colorbar;colormap default;
                title('Use LMB to select fiducials. Hit enter to stop or use RMB to select the final fiducial.');
                hold on
                [xi,yi,but] = ginput(1);
                text(xi,yi,'*','color','green','fontsize',20,'fontweight','bold');
                moleLocs(n,:) = round(moleLocs(n,:)-ROIsize+[xi,yi]-1);
                hold off;
            end
            close(hLocs_1)
        end
    else
        dataBlur = imfilter(dataAvg, fspecial('disk', boxRadius), 'same', 'conv');
        dataBlur(dataBlur < median(dataBlur(:))) = median(dataBlur(:));
        [tempY, tempX] = ind2sub(size(dataAvg),find(imregionalmax(dataBlur)));
        moleLocs = [tempX,tempY];
        %% Clean up duplicated locs
        dupes = abs(bsxfun(@minus,moleLocs(:,1),moleLocs(:,1)'))<boxRadius & abs(bsxfun(@minus,moleLocs(:,2),moleLocs(:,2)'))<boxRadius;
        dupes2 = unique(dupes(sum(dupes)>1,:),'rows');
        moleLocs2 = [];
        for n = 1:size(dupes2,1)
            moleLocsT = moleLocs(dupes2(n,:),:);
            moleLocsS = nan(size(moleLocsT,1),1);
            for i = 1:size(moleLocsT,1)
                moleLocsS(i,:) = dataBlur(moleLocsT(i,2),moleLocsT(i,1));
            end
            moleLocs2(n,:) = round(sum(bsxfun(@times, moleLocsT,moleLocsS))/sum(moleLocsS));
        end
        moleLocs(any(dupes2),:) = [];
        moleLocs = [moleLocs;moleLocs2];
        clear moleLocsS moleLocsT moleLocs2
        for n = 1:length(moleLocs)
            text(moleLocs(n,1),moleLocs(n,2),num2str(n),'color','white','fontsize',13,'fontweight','bold');
        end
    end
    saveas(hLocs,[outputFilePrefix 'bead map.png']);
    close(hLocs);
    clear dataBlur dataAvg tempY tempX mex;
    moleLocs(moleLocs(:,1)<2*boxRadius+2|moleLocs(:,2)<2*boxRadius+2|moleLocs(:,1)>ROI(3)-2*boxRadius-2|moleLocs(:,2)>ROI(4)-2*boxRadius-2,:) = [];
    %% Initialize data arrays
    % profile on
    % [frameNum moleNum amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2
    %  bkgndMean totalFitError goodFit xCenter yCenter angle numPhotons CFocusPosition]
    
    numBeads = size(moleLocs,1);
    PSFfits = nan(numFramesTotal*numBeads, 21);
    PSFfitsA = cell(numBeads,1);
    % phasemask_position = zeros(numFiles,1);
    startTime = tic;
    
    % startFrame = zeros(numFiles,200);
    % endFrame = zeros(numFiles,200);
    c=0;
    var = double((variance./gain)./gain);
    datavar = var(ROI(2):ROI(2)+ROI(4)-1, ...
        ROI(1):ROI(1)+ROI(3)-1);
    
    boxLength = 2*boxRadius+1;
    TxIdx = bsxfun(@plus,(-boxRadius:boxRadius)',moleLocs(:,1)');
    TyIdx = bsxfun(@plus,(-boxRadius:boxRadius)',moleLocs(:,2)');
    lB = [zeros(numBeads,2), min(TxIdx)', min(TyIdx)', min(TxIdx)', min(TyIdx)', repmat(sigmaBounds(1),numBeads,2), zeros(numBeads,1)];
    uB = [zeros(numBeads,2), max(TxIdx)', max(TyIdx)', max(TxIdx)', max(TyIdx)', repmat(sigmaBounds(2),numBeads,2), zeros(numBeads,1)];
    if numFiles > 1
        PSFfitsZ = cell(numFiles,1);
    end
    
    for n = 1:numFiles
        if isDcimg
            numFrames = double(totalframes);
        else
            dataFileInfo = imfinfo(dataFile{n});
            numFrames = double(length(dataFileInfo));
        end
        %% Fit chosen beads throughout entire image stack
        h = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        
        for a=1:numFrames%size(sifLogData,1)-2
            % for a=1:130
            %var = double((variance./gain)./gain);
            if isDcimg
                [framedata,totalframes]= dcimgmatlab(a-1, dcimgfile);
                framedatatrans = transpose (framedata);
                dataFrame = double(framedatatrans)-darkAvg;
                dataFrame = double(dataFrame./gain);
                data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
                    ROI(1):ROI(1)+ROI(3)-1);
                clear mex;
            else
                dataFrame = double(imread(dataFile{n},a,'Info',dataFileInfo))-darkAvg;
                dataFrame = double(dataFrame./gain);
                data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
                    ROI(1):ROI(1)+ROI(3)-1);
            end
            % subtract the background and continue
            if useWaveSub
                bkgndImg = f_waveletBackground(data);
            else
                % subtract median to get rough idea of #photons of fiducial
                bkgndImg = repmat(median(data(:)),size(data));
                %%Note: bkgnd handled differently by CMD
                %bkgndImg2= mean(mean(bkgndImg));
            end
            %%Note: bkgnd handled differently by CMD
            switch(fittingMethod)
                case 'LSQ with DG model'
                    data = data - bkgndImg;
                case 'MLE with DG model'
                    data = data + datavar;
            end
            % blur data for more robust peak finding
            dataBlur = imfilter(data, fspecial('gaussian', size(data), blurSize), 'same', 'conv');
            
            %% do fitting to extract exact locations of DH-PSFs
            % create reconstructed DH-PSF image from fitted data
            reconstructImg = zeros(cropHeight, cropWidth);
            PSFfitsT = nan(numBeads,21);
            PSFfitsS = cell(numBeads,1);
            if numBeads > 24
                parfor b=1:numBeads
                    if any(TxIdx(:,b)<1) || any(TyIdx(:,b)<1) || any(TxIdx(:,b)>ROI(3)) || any(TyIdx(:,b)>ROI(4))
                        PSFfitsS{b,1} = [a+c, b,nan(1,10), -1000];
                        continue;
                    end
                    PSFfitsS{b,1} = PSFfitting_calDHPSF(data,fittingMethod,cropHeight,cropWidth,TxIdx,TyIdx,sigmaBounds,moleLocs,dataBlur,datavar,bkgndImg,boxLength,a,c,b,uB,lB);
                end
            else
                for b=1:numBeads
                    if any(TxIdx(:,b)<1) || any(TyIdx(:,b)<1) || any(TxIdx(:,b)>ROI(3)) || any(TyIdx(:,b)>ROI(4))
                        PSFfitsS{b,1} = [a+c, b,nan(1,10), -1000];
                        continue;
                    end
                    PSFfitsS{b,1} = PSFfitting_calDHPSF(data,fittingMethod,cropHeight,cropWidth,TxIdx,TyIdx,sigmaBounds,moleLocs,dataBlur,datavar,bkgndImg,boxLength,a,c,b,uB,lB);
                end
            end
            PSFfitsT(:,1:13) = cell2mat(PSFfitsS);
            for b = 1:numBeads
                if min(TxIdx(:,b)) < 1 || max(TyIdx(:,b)) > cropHeight || max(TxIdx(:,b)) > cropWidth || min(TyIdx(:,b)) < 1
                    PSFfitsT(b,13) = -1000;
                    PSFfitsT(b,17) = NaN;
                    continue
                end
                %Below is a way to count the photons. It integrates the box and
                %subtracts the boxarea*offset from the fit. It is inherently flawed
                %if there happen to be bright pixels inside of the fitting region.
                %Units of photons for LSQ nonlin
                PSFfitsT(b,17) = sum(sum(data(TyIdx(:,b),TxIdx(:,b))));
            end
            % Calculate midpoint between two Gaussian spots
            % convert from pixels to nm
            PSFfitsT(:,14) = (PSFfitsT(:,5)+PSFfitsT(:,7))/2*nmPerPixel;
            PSFfitsT(:,15) = (PSFfitsT(:,6)+PSFfitsT(:,8))/2*nmPerPixel;
            
            % Below is the calculation of the angle of the two lobes.
            % Remember that two vertical lobes is focal plane because camera
            % outputs data that is rotated. Therefore, we want y2>y1 for all
            % angle calculations (so that -90<=angle<=90, and we use swap
            % the use of x and y for the atan2 calculation.
            if ~horizontal
                flop = PSFfitsT(:,8)>PSFfitsT(:,6);
                PSFfitsT(flop,16) = atan2(-(PSFfitsT(flop,7)-PSFfitsT(flop,5)),PSFfitsT(flop,8)-PSFfitsT(flop,6)) * 180/pi;
                PSFfitsT(~flop,16) = atan2(-(PSFfitsT(~flop,5)-PSFfitsT(~flop,7)),PSFfitsT(~flop,6)-PSFfitsT(~flop,8)) * 180/pi;
            elseif horizontal
                flop = PSFfitsT(:,7)>PSFfitsT(:,5);
                PSFfitsT(flop,16) = atan2(-(PSFfitsT(flop,8)-PSFfitsT(flop,6)),PSFfitsT(flop,7)-PSFfitsT(flop,5)) * 180/pi;
                PSFfitsT(~flop,16) = atan2(-(PSFfitsT(~flop,6)-PSFfitsT(~flop,8)),PSFfitsT(~flop,5)-PSFfitsT(~flop,7)) * 180/pi;
            end
            
            %The interlobe distance
            PSFfitsT(:,19) = sqrt((PSFfitsT(:,5)-PSFfitsT(:,7)).^2 + (PSFfitsT(:,6)-PSFfitsT(:,8)).^2);
            %Amplitude Ratio
            PSFfitsT(:,20) = abs(PSFfitsT(:,3) - PSFfitsT(:,4))./sum(PSFfitsT(:,3:4),2);
            % Gaussian width Ratio
            PSFfitsT(:,21) = abs(PSFfitsT(:,9) - PSFfitsT(:,10))./sum(PSFfitsT(:,9:10),2);
            
            
            %% Now evaluate the fits
            % Conditions for fits (play with these):
            % (1) Amplitude of both lobes > 0
            % (2) All locations x1,y1, x2,y2 lie inside area of small box
            % (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
            % (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
            % (5) Make sure amplitudes are within 100% of one another
            % (6) Make sure totalFitError/(total number of photons) < 1.05 (not
            %     enabled at the present time)
            
            
            % absolute amplitude > 0?
            PSFfitsT(PSFfitsT(:,3)<0,13) = -1001;
            PSFfitsT(PSFfitsT(:,4)<0,13) = -1001;
            % peaks inside box?
            % sigma ratio of lobes less than limit?
            PSFfitsT(PSFfitsT(:,21) > sigmaRatioLimit,13) = -1004;
            % interlobe distance within bounds?
            PSFfitsT(PSFfitsT(:,19) > lobeDistBounds(2),13) = -1005;
            PSFfitsT(PSFfitsT(:,19) < lobeDistBounds(1),13) = -1005;
            % amplitude ratio of lobes less than limit?
            PSFfitsT(PSFfitsT(:,20) > ampRatioLimit,13) = -1006;
            %             if PSFfits(rowIdx,12)*conversionFactor/PSFfits(rowIdx,17) > 1 || ...
            %                PSFfits(rowIdx,12)*conversionFactor/PSFfits(rowIdx,17) < 0
            %                 PSFfits(rowIdx,13) = -1007;
            %             end
            % does the fit jump a large amount (>400 nm) from previous fit?
            %         shifted = (PSFfitsT(:,14)-moleLocs(:,1).*nmPerPixel).^2 + (PSFfitsT(:,15)-moleLocs(:,2).*nmPerPixel).^2 > 400^2;
            %         PSFfitsT(shifted,13) = -1008;
            %         FIX this section ASAP hkhk
            % if fit was successful, use the computed center location as center
            % of box for next iteration
            if numBeads > 10
                PSFfitsT = RefitErrors_2(PSFfitsT,TxIdx,TyIdx,data,bkgndImg,datavar,fittingMethod,numBeads,nmPerPixel,horizontal,boxLength,a,c, sigmaBounds,lB,uB,sigmaRatioLimit,lobeDistBounds,ampRatioLimit, cropHeight, cropWidth);
            end
            %The actual position read by the C-focus encoder
            PSFfitsT(:,18) = repmat(sifLogData1(a+c,4),numBeads,1); %c adds the length of the previous file to a
            moleLocs(PSFfitsT(:,13) > 0,:) = round(PSFfitsT(PSFfitsT(:,13)>0,14:15)./nmPerPixel);
            
            TxIdx(:,PSFfitsT(:,13) > 0) = bsxfun(@plus,(-boxRadius:boxRadius)',moleLocs(PSFfitsT(:,13) > 0,1)');
            TyIdx(:,PSFfitsT(:,13) > 0) = bsxfun(@plus,(-boxRadius:boxRadius)',moleLocs(PSFfitsT(:,13) > 0,2)');
            lB = [zeros(numBeads,2), min(TxIdx)', min(TyIdx)', min(TxIdx)', min(TyIdx)', repmat(sigmaBounds(1),numBeads,2), zeros(numBeads,1)];
            uB = [zeros(numBeads,2), max(TxIdx)', max(TyIdx)', max(TxIdx)', max(TyIdx)', repmat(sigmaBounds(2),numBeads,2), zeros(numBeads,1)];
            PSFfitsA{a,1} = PSFfitsT;
            
            %% plot image reconstruction so that fits can be checked
            [xIdx, yIdx] = meshgrid(1:cropWidth,1:cropHeight);
            for b = 1:numBeads
                fitParam = PSFfitsT(b,3:10);
                reconstructImg = reconstructImg + ...
                    fitParam(1).*exp( -((xIdx-fitParam(3)).^2+(yIdx-fitParam(4)).^2.) / (2.*fitParam(7).^2)) ...
                    +fitParam(2).*exp( -((xIdx-fitParam(5)).^2+(yIdx-fitParam(6)).^2.) / (2.*fitParam(8).^2));
            end
            %  plot results of template matching and fitting
            set(0,'CurrentFigure',h);
            
            if ~isempty(PSFfitsT)
                subplot('Position',[0.025 0.025 .9/2 .95]);
%                 imagesc(data,[min(bkgndImg(:)),max([max(PSFfitsT(:,3)),max(PSFfitsT(:,4)),max(bkgndImg(:))])]);axis image;colormap parula;
                imagesc(data);axis image;colormap parula;
                title(['Frame ' num2str(a+c) ': raw data - dark counts']);
                
                subplot('Position',[0.525 0.025 .9/2 .95]);
%                 imagesc(reconstructImg+bkgndImg,[min(bkgndImg(:)),max([max(PSFfitsT(:,3)),max(PSFfitsT(:,4)),max(bkgndImg(:))])]);
                imagesc(reconstructImg+bkgndImg);
                axis image;
                title('Image reconstructed from fitted matches');
                
                drawnow;
            end
        end
        
        c=c+numFrames;
        elapsedTime = toc(startTime);
        clear data residual dataAvg reconstructImg xIdx yIdx temp;
        close(h);
        fps = numFrames/elapsedTime
        beadsPerSec = numFrames*numBeads/elapsedTime
        
        if numFiles == 1
            PSFfits = cell2mat(PSFfitsA);
        elseif n == numFiles
            PSFfitsZ{n,1} = cell2mat(PSFfitsA);
            PSFfits = cell2mat(PSFfitsZ);
        else
            PSFfitsZ{n,1} = cell2mat(PSFfitsA);
        end
    end
    
    absLocs = zeros(numBeads,2); % actual positions in nm
    for b = 1:numBeads
        absLocs(b,1:2) = mean(squeeze(PSFfits(PSFfits(:,2)==b,14:15)),1) + [ROI(1)-1+iROI(1)-1,ROI(2)-1+iROI(2)-1] * nmPerPixel;
    end
    % save fit info to MATLAB mat file
    save([outputFilePrefix 'raw fits.mat']);
else
    [rawFile, rawPath] = uigetfile({'*.mat'},'Please pick the Raw Fits file.');
    if isequal(rawFile,0)
        error('User cancelled the program');
    end
    load([rawPath, rawFile]);
    outputFilePrefix = rawPath;
end
%% Processing setup

meanAngles = nan(numBeads,numSteps);
stddevAngles = nan(numBeads,numSteps);
meanPhotons = nan(numBeads,numSteps);
stddevPhotons = nan(numBeads,numSteps);
meanX = nan(numBeads,numSteps);
stdX = nan(numBeads,numSteps);
meanY = nan(numBeads,numSteps);
stdY = nan(numBeads,numSteps);
numGoodFrames = nan(numBeads,numSteps);

meanInterlobeDistance = nan(numBeads,numSteps);
stdInterlobeDistance = nan(numBeads,numSteps);
meanAmp1 = nan(numBeads,numSteps);
stdAmp1 = nan(numBeads,numSteps);
meanAmp2 = nan(numBeads,numSteps);
stdAmp2 = nan(numBeads,numSteps);
meanAmpRatio = nan(numBeads,numSteps);
stdAmpRatio = nan(numBeads,numSteps);
meanSigma1 = nan(numBeads,numSteps);
stdSigma1 = nan(numBeads,numSteps);
meanSigma2 = nan(numBeads,numSteps);
stdSigma2 = nan(numBeads,numSteps);
meanSigmaRatio = nan(numBeads,numSteps);
stdSigmaRatio = nan(numBeads,numSteps);
goodFit_f = zeros(numBeads,numSteps);
goodFit_b = zeros(numBeads,numSteps);
zAngleZero = nan(numBeads,1);
%     figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');

textHeader = {'Start Frame' 'End Frame' 'Mean Angle (deg)' ...
    'x deviation (nm)' 'y deviation (nm)' 'z position (nm)' 'Mean number of Photons' ...
    'Angle Std Dev (deg)' 'x std dev (nm)' 'y std dev (nm)' 'z Std Dev (nm)' ...
    'Photons Std Dev' 'Num good frames'};

stepSize =  diff(meanCFocPos(:));
z0 = meanCFocPos(:) - meanCFocPos(1);
z1 = z0(2:length(z0)-1);
% flip definition of z if angle vs z slope is positive (should be
% negative)
[~, centralBeadIdx] = min((absLocs(:,1)/nmPerPixel-(size(gain,1)+iROI(1)-1)/2).^2+(absLocs(:,2)/nmPerPixel-(size(gain,2)+iROI(2)-1)/2).^2);
% [~, centralBeadIdx] = min((absLocs(:,1)/nmPerPixel-size(gain,1)/2).^2+(absLocs(:,2)/nmPerPixel-size(gain,2)/2).^2);

for bead = 1:numBeads
    % extract fitting parameters for this bead
    beadFitParam(:,:) = PSFfits(PSFfits(:,2) == bead, :);
    for step = 1:numSteps
        %     for step = 8:18
        % extract bead fitting parameters for this step
        stepFitParam = beadFitParam(startFrame(step):endFrame(step),:);
        
        x = stepFitParam(:,14);
        y = stepFitParam(:,15);
        angles = stepFitParam(:,16);
        
        % correct noise in angle measurement if oscillating between +90
        % and -90 degrees. 
        if (~isempty(angles(angles > 80)) && ~isempty(angles(angles < -80)))
            angles(angles < 0) = angles(angles < 0) + 180;
        end
        
        amp1 = stepFitParam(:,3);
        amp2 = stepFitParam(:,4);
        sigma1 = stepFitParam(:,9);
        sigma2 = stepFitParam(:,10);
        numPhotons = stepFitParam(:,17);
        lobeDist = stepFitParam(:,19);
        ampRatio = stepFitParam(:,20);
        sigmaRatio = stepFitParam(:,21);
        
        goodFit = stepFitParam(:,13)>0;
        
        meanX(bead,step) = mean(x(goodFit));
        stdX(bead,step) = std(x(goodFit));
        meanY(bead,step) = mean(y(goodFit));
        stdY(bead,step) = std(y(goodFit));
        meanAngles(bead,step) = mean(angles(goodFit));
        stddevAngles(bead,step) = std(angles(goodFit));
        meanPhotons(bead,step) = mean(numPhotons(goodFit));
        stddevPhotons(bead,step) = std(numPhotons(goodFit));
        numGoodFrames(bead,step) = length(angles(goodFit));
        
        meanInterlobeDistance(bead,step) = mean(lobeDist(goodFit));
        stdInterlobeDistance(bead,step) = std(lobeDist(goodFit));
        meanAmp1(bead,step) = mean(amp1(goodFit));
        stdAmp1(bead,step) = std(amp1(goodFit));
        meanAmp2(bead,step) = mean(amp2(goodFit));
        stdAmp2(bead,step) = std(amp2(goodFit));
        meanSigma1(bead,step) = mean(sigma1(goodFit));
        stdSigma1(bead,step) = std(sigma1(goodFit));
        meanSigma2(bead,step) = mean(sigma2(goodFit));
        stdSigma2(bead,step) = std(sigma2(goodFit));
        meanAmpRatio(bead,step) = mean(ampRatio(goodFit));
        stdAmpRatio(bead,step) = std(ampRatio(goodFit));
        meanSigmaRatio(bead,step) = mean(sigmaRatio(goodFit));
        stdSigmaRatio(bead,step) = std(sigmaRatio(goodFit));
        
    end
    
    %unwrap angles if we have greater than 180 degree range; need to convert
    %180 degrees to 2pi because unwrap function works on intervals of 2pi, not
    %pi
    meanAngles(bead,:) = unwrap(meanAngles(bead,:)*2*pi/180)*180/(2*pi);
    if isempty(meanAngles(bead,meanAngles(bead,:)>10))
        meanAngles(bead,:) = meanAngles(bead,:) + 180;
    end
    if isempty(meanAngles(bead,meanAngles(bead,:)<-10))
        meanAngles(bead,:) = meanAngles(bead,:) - 180;
    end
    % this computes the differences in angle between adjacent steps
    % these differences are sampled at half points between steps
    angleSlope = squeeze(squeeze(diff(meanAngles(bead,:))));
    % average adjacent differences to resample differences at integer points,
    % convert to units of angle/nm
    indices = 1:length(angleSlope)-1;
    anglePerNM = (angleSlope(indices+1)+angleSlope(indices))./(2*stepSize(indices)');  %check if correct
    stddevNM(bead,:) = stddevAngles(bead, 2:length(stddevAngles(bead,:))-1)./anglePerNM*sign(nanMedian(anglePerNM));
    
    % test to make sure angle vs z curve is monotonic -- remove any
    % outlying points
    %     slope = sign(meanAngles(bead,length(meanAngles(bead,:)))-meanAngles(bead,1));
    goodFit = ~isnan(meanAngles(bead,:));
    
    stepSlope = false(1,numSteps);
    indices_forward = find(sign(angleSlope)==-1);   % used to be 1 for 8a back
    stepSlope(indices_forward)=true;
    goodFit_forward = stepSlope & goodFit;
    stepSlope = false(1,numSteps);
    indices_backward = find(sign(angleSlope)==1);   % used to be -1 for 8a back
    stepSlope(indices_backward)=true;
    goodFit_backward = stepSlope & goodFit;
    
    % if calibration movie only contains one backward scan,
    % then replicate it as if it were moving forward
    % Try to use forward frames unless only a few are fit, or many more
    % backward are fit.
    if sum(goodFit_backward) > sum(goodFit_forward)*2
        gfTemp = goodFit_forward;
        goodFit_forward = goodFit_backward;
        goodFit_backward = gfTemp;
        clear gfTemp
    end
    
    %     goodFit = logical([0 ones(1,13) zeros(1,11)]) & ...
    %         squeeze(squeeze(~isnan(meanAngles(bead,:))))' & ...
    %         [true squeeze(sign(meanAngles(bead,2:size(meanAngles,3))-...
    %         meanAngles(bead,1:size(meanAngles,3)-1)))'==slope];
    
    %     goodFit = goodFit(2:length(goodFit)-1);
    %     z = 0:stepSize:stepSize*(numSteps-1);
    % if there aren't enough points for a good curve, skip this fiduciary
    if sum(goodFit_forward) < 5
        disp(['Fiducial number ' num2str(bead) ' had only ' num2str(sum(goodFit_forward)) ' usable steps and is excluded']);
        continue;
    end
    
    % compute xyz vector so that z = 0 corresponds to angle = 0
    xAngleZero = interp1(squeeze(meanAngles(bead,goodFit_forward)),...
        squeeze(meanX(bead,goodFit_forward)),0,'pchip');
    yAngleZero = interp1(squeeze(meanAngles(bead,goodFit_forward)),...
        squeeze(meanY(bead,goodFit_forward)),0,'pchip');
    zAngleZero(bead) = interp1(squeeze(meanAngles(bead,goodFit_forward)),...
        z0(goodFit_forward),0,'pchip');
    meanX_absolute(bead,goodFit) = meanX(bead,goodFit);
    meanY_absolute(bead,goodFit) = meanY(bead,goodFit);
    meanX(bead,goodFit) = meanX(bead,goodFit) - xAngleZero;
    meanY(bead,goodFit) = meanY(bead,goodFit) - yAngleZero;
    %% move 'local' calibration parameters to 'global' array
    goodFit_f(bead,:) = goodFit_forward;
    goodFit_b(bead,:) = goodFit_backward;
    
end
slope = sign(meanAngles(centralBeadIdx,round(numSteps*1/2)+3)-meanAngles(centralBeadIdx,round(numSteps*1/2)-3));
if slope>0
    z0 = -z0;
    z1 = -z1;
end
GoodInterp = zAngleZero < max(z0) & zAngleZero > min(z0);
maxZzero = max(zAngleZero(GoodInterp));
z0 = z0-maxZzero;
z1 = z1-maxZzero;
z = z0;
%% write calibration parameters to a file
save([outputFilePrefix 'calibration.mat'], ...
    'meanAngles', 'meanX', 'meanY', 'meanX_absolute', 'meanY_absolute', 'z','xAngleZero','yAngleZero',...
    'zAngleZero', 'goodFit_f', 'absLocs',...
    'goodFit_b', 'meanPhotons', 'stddevPhotons','stdX','stdY',...
    'stddevAngles','meanInterlobeDistance','stdInterlobeDistance',...
    'meanAmpRatio','stdAmpRatio','meanAmp1','meanAmp2','meanCFocPos');

if numBeads < 25
    %% plot calibration parameters
    h1 = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    for bead = 1:numBeads
        goodFit_backward = logical(goodFit_b(bead,:)); goodFit_forward= logical(goodFit_f(bead,:));
        if sum(goodFit_forward) < 5
            continue;
        end
        if logical(sum(goodFit_backward))
            scanLegend = {'forward scan','backward scan'};
            precLegend = {'x forward','x backward','y forward','y backward'};
        else
            scanLegend = {'calibration scan'};
            precLegend = {'x precision','y precision'};
        end
        
        subplot(2,3,1:2);
        errorbar(z0(goodFit_forward),meanAngles(bead,goodFit_forward),stddevAngles(bead,goodFit_forward),'-');
        if logical(sum(goodFit_backward))
            hold on
            errorbar(z0(goodFit_backward),meanAngles(bead,goodFit_backward),stddevAngles(bead,goodFit_backward),':');
            hold off
        end
        axis tight;
        title({dataFile1{1} ['Fiduciary ' num2str(bead)]},'interpreter','none');
        legend(scanLegend);
        xlabel('z position (nm)');
        ylabel('Angle (deg)');
        
        subplot(2,3,3);
        plot(z0(goodFit_forward),meanX(bead,goodFit_forward)','b');
        hold on;
        plot(z0(goodFit_forward),meanY(bead,goodFit_forward)','r');
        if logical(sum(goodFit_backward))
            plot(z0(goodFit_backward),meanX(bead,goodFit_backward)',':b');
            plot(z0(goodFit_backward),meanY(bead,goodFit_backward)',':r');
        end
        axis tight;
        %     legend('x forward','x backward','y forward','y backward');
        xlabel('z position (nm)');
        ylabel('xy position (nm)');
        hold off;
        
        subplot(2,3,4);
        plot(z0(goodFit_forward),squeeze(stdX(bead,goodFit_forward)),'b');
        hold on;
        plot(z0(goodFit_forward),squeeze(stdY(bead,goodFit_forward)),'r');
        if logical(sum(goodFit_backward))
            plot(z0(goodFit_backward),squeeze(stdX(bead,goodFit_backward)),':b');
            plot(z0(goodFit_backward),squeeze(stdY(bead,goodFit_backward)),':r');
        end
        axis tight;
        ylim([0 40]);
        legend(precLegend);
        xlabel('z position (nm)');
        ylabel('localization precision (nm)');
        hold off;
        
        subplot(2,3,5);
        plot(z1(goodFit_forward(2:length(goodFit_forward)-1)),stddevNM(bead,goodFit_forward(2:length(goodFit_forward)-1)),'b');
        if logical(sum(goodFit_backward))
            hold on
            plot(z1(goodFit_backward(2:length(goodFit_backward)-1)),stddevNM(bead,goodFit_backward(2:length(goodFit_backward)-1)),':b');
            hold off
        end
        axis tight;
        ylim([0 40]);
        xlabel('z position (nm)');
        ylabel('z localization precision (nm)');
        
        subplot(2,3,6);
        errorbar(z0(goodFit_forward),meanPhotons(bead,goodFit_forward),stddevPhotons(bead,goodFit_forward),'b');
        if logical(sum(goodFit_backward))
            hold on
            errorbar(z0(goodFit_backward),meanPhotons(bead,goodFit_backward),stddevPhotons(bead,goodFit_backward),':b');
            hold off
        end
        axis tight;
        xlabel('z position (nm)');
        ylabel('number of photons');
        drawnow
        print(h1,'-dpng',[outputFilePrefix 'bead ' num2str(bead) ' stats_1.png']);
        % close this window if there are a lot of beads that were fitted
    end
    if numBeads > 10
        close(h1);
    end
    %% plot additional pertaining to the DH-PSF
    h2=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    for bead = 1:numBeads
        goodFit_backward = logical(goodFit_b(bead,:)); goodFit_forward= logical(goodFit_f(bead,:));
        if sum(goodFit_forward) < 5
            continue;
        end
        subplot(2,3,1);
        errorbar(z0(goodFit_forward),meanAngles(bead,goodFit_forward),stddevAngles(bead,goodFit_forward),'-');
        hold on
        %     errorbar(z(goodFit_backward),meanAngles(bead,goodFit_backward),stddevAngles(bead,goodFit_backward),':');
        axis tight;
        title({dataFile1{1} ['Fiduciary ' num2str(bead)]},'interpreter','none');
        legend('forward scan'); %,'backward scan');
        xlabel('z position (nm)');
        ylabel('Angle (deg)');
        hold off
        
        subplot(2,3,4);
        errorbar(z0(goodFit_forward),meanInterlobeDistance(bead,goodFit_forward),stdInterlobeDistance(bead,goodFit_forward),'-');
        axis tight;
        xlabel('z position (nm)');
        ylabel('Interlobe Distance (pix)');
        hold off;
        
        subplot(2,3,2);
        errorbar(z0(goodFit_forward),meanAmp1(bead,goodFit_forward),stdAmp1(bead,goodFit_forward),'-');
        hold on;
        errorbar(z0(goodFit_forward),meanAmp2(bead,goodFit_forward),stdAmp2(bead,goodFit_forward),'-r');
        axis tight;
        legend('Lobe 1','Lobe 2');
        xlabel('z position (nm)');
        ylabel('Gaussian Amplitudes (counts)');
        hold off;
        
        subplot(2,3,5);
        errorbar(z0(goodFit_forward),meanAmpRatio(bead,goodFit_forward),stdAmpRatio(bead,goodFit_forward),'-');
        hold on;
        axis tight;
        xlabel('z position (nm)');
        ylabel('Gaussian Amplitude Ratio');
        hold off;
        
        subplot(2,3,3);
        errorbar(z0(goodFit_forward),meanSigma1(bead,goodFit_forward),stdSigma1(bead,goodFit_forward),'-');
        hold on;
        errorbar(z0(goodFit_forward),meanSigma2(bead,goodFit_forward),stdSigma2(bead,goodFit_forward),'-r');
        axis tight;
        legend('Lobe 1','Lobe 2');
        xlabel('z position (nm)');
        ylabel('Gaussian Sigma (pixel)');
        hold off;
        
        subplot(2,3,6);
        errorbar(z0(goodFit_forward),meanSigmaRatio(bead,goodFit_forward),stdSigmaRatio(bead,goodFit_forward),'-');
        hold on;
        axis tight;
        xlabel('z position (nm)');
        ylabel('Gaussian Sigma Ratio');
        hold off;
        
        print(h2,'-dpng',[outputFilePrefix 'bead ' num2str(bead) ' stats_2.png']);
    end
    if numBeads > 10
        close(h2);
    end
else
    %% plot calibration parameters
    h1 = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    parfor bead = 1:numBeads
        goodFit_backward = logical(goodFit_b(bead,:)); goodFit_forward= logical(goodFit_f(bead,:));
        if sum(goodFit_forward) < 5
            continue;
        end
        if logical(sum(goodFit_backward))
            scanLegend = {'forward scan','backward scan'};
            precLegend = {'x forward','x backward','y forward','y backward'};
        else
            scanLegend = {'calibration scan'};
            precLegend = {'x precision','y precision'};
        end
        
        subplot(2,3,1:2);
        errorbar(z0(goodFit_forward),meanAngles(bead,goodFit_forward),stddevAngles(bead,goodFit_forward),'-');
        if logical(sum(goodFit_backward))
            hold on
            errorbar(z0(goodFit_backward),meanAngles(bead,goodFit_backward),stddevAngles(bead,goodFit_backward),':');
            hold off
        end
        axis tight;
        title({dataFile1{1} ['Fiduciary ' num2str(bead)]},'interpreter','none');
        legend(scanLegend);
        xlabel('z position (nm)');
        ylabel('Angle (deg)');
        
        subplot(2,3,3);
        plot(z0(goodFit_forward),meanX(bead,goodFit_forward)','b');
        hold on;
        plot(z0(goodFit_forward),meanY(bead,goodFit_forward)','r');
        if logical(sum(goodFit_backward))
            plot(z0(goodFit_backward),meanX(bead,goodFit_backward)',':b');
            plot(z0(goodFit_backward),meanY(bead,goodFit_backward)',':r');
        end
        axis tight;
        %     legend('x forward','x backward','y forward','y backward');
        xlabel('z position (nm)');
        ylabel('xy position (nm)');
        hold off;
        
        subplot(2,3,4);
        plot(z0(goodFit_forward),squeeze(stdX(bead,goodFit_forward)),'b');
        hold on;
        plot(z0(goodFit_forward),squeeze(stdY(bead,goodFit_forward)),'r');
        if logical(sum(goodFit_backward))
            plot(z0(goodFit_backward),squeeze(stdX(bead,goodFit_backward)),':b');
            plot(z0(goodFit_backward),squeeze(stdY(bead,goodFit_backward)),':r');
        end
        axis tight;
        ylim([0 40]);
        legend(precLegend);
        xlabel('z position (nm)');
        ylabel('localization precision (nm)');
        hold off;
        
        subplot(2,3,5);
        plot(z1(goodFit_forward(2:length(goodFit_forward)-1)),stddevNM(bead,goodFit_forward(2:length(goodFit_forward)-1)),'b');
        if logical(sum(goodFit_backward))
            hold on
            plot(z1(goodFit_backward(2:length(goodFit_backward)-1)),stddevNM(bead,goodFit_backward(2:length(goodFit_backward)-1)),':b');
            hold off
        end
        axis tight;
        ylim([0 40]);
        xlabel('z position (nm)');
        ylabel('z localization precision (nm)');
        
        subplot(2,3,6);
        errorbar(z0(goodFit_forward),meanPhotons(bead,goodFit_forward),stddevPhotons(bead,goodFit_forward),'b');
        if logical(sum(goodFit_backward))
            hold on
            errorbar(z0(goodFit_backward),meanPhotons(bead,goodFit_backward),stddevPhotons(bead,goodFit_backward),':b');
            hold off
        end
        axis tight;
        xlabel('z position (nm)');
        ylabel('number of photons');
        drawnow
        print(h1,'-dpng',[outputFilePrefix 'bead ' num2str(bead) ' stats_1.png']);
        % close this window if there are a lot of beads that were fitted
    end
    close(h1);
    %% plot additional pertaining to the DH-PSF
    h2=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    parfor bead = 1:numBeads
        goodFit_backward = logical(goodFit_b(bead,:)); goodFit_forward= logical(goodFit_f(bead,:));
        if sum(goodFit_forward) < 5
            continue;
        end
        subplot(2,3,1);
        errorbar(z0(goodFit_forward),meanAngles(bead,goodFit_forward),stddevAngles(bead,goodFit_forward),'-');
        hold on
        %     errorbar(z(goodFit_backward),meanAngles(bead,goodFit_backward),stddevAngles(bead,goodFit_backward),':');
        axis tight;
        title({dataFile1{1} ['Fiduciary ' num2str(bead)]},'interpreter','none');
        legend('forward scan'); %,'backward scan');
        xlabel('z position (nm)');
        ylabel('Angle (deg)');
        hold off
        
        subplot(2,3,4);
        errorbar(z0(goodFit_forward),meanInterlobeDistance(bead,goodFit_forward),stdInterlobeDistance(bead,goodFit_forward),'-');
        axis tight;
        xlabel('z position (nm)');
        ylabel('Interlobe Distance (pix)');
        hold off;
        
        subplot(2,3,2);
        errorbar(z0(goodFit_forward),meanAmp1(bead,goodFit_forward),stdAmp1(bead,goodFit_forward),'-');
        hold on;
        errorbar(z0(goodFit_forward),meanAmp2(bead,goodFit_forward),stdAmp2(bead,goodFit_forward),'-r'); %#ok<*PFBNS>
        axis tight;
        legend('Lobe 1','Lobe 2');
        xlabel('z position (nm)');
        ylabel('Gaussian Amplitudes (counts)');
        hold off;
        
        subplot(2,3,5);
        errorbar(z0(goodFit_forward),meanAmpRatio(bead,goodFit_forward),stdAmpRatio(bead,goodFit_forward),'-');
        hold on;
        axis tight;
        xlabel('z position (nm)');
        ylabel('Gaussian Amplitude Ratio');
        hold off;
        
        subplot(2,3,3);
        errorbar(z0(goodFit_forward),meanSigma1(bead,goodFit_forward),stdSigma1(bead,goodFit_forward),'-');
        hold on;
        errorbar(z0(goodFit_forward),meanSigma2(bead,goodFit_forward),stdSigma2(bead,goodFit_forward),'-r');
        axis tight;
        legend('Lobe 1','Lobe 2');
        xlabel('z position (nm)');
        ylabel('Gaussian Sigma (pixel)');
        hold off;
        
        subplot(2,3,6);
        errorbar(z0(goodFit_forward),meanSigmaRatio(bead,goodFit_forward),stdSigmaRatio(bead,goodFit_forward),'-');
        hold on;
        axis tight;
        xlabel('z position (nm)');
        ylabel('Gaussian Sigma Ratio');
        hold off;
        
        print(h2,'-dpng',[outputFilePrefix 'bead ' num2str(bead) ' stats_2.png']);
    end
    close(h2);
end

if numBeads < 25
    %% Generate template Stack of a chosen bead
    templateSize = 20;      % 2*round(10*160/nmPerPixel);  % 26;
    numFrames_temp = double(totalframes);
    if ~isDcimg
        dataFileInfo = imfinfo(dataFile{1});
        numFrames_temp =  length(dataFileInfo);
    end
    for bead = 1:numBeads
        numTemplates=length(startFrame(logical((goodFit_f(bead, :)))));
        forwardStartFrames = startFrame(logical((goodFit_f(bead, :))));
        forwardEndFrames = endFrame(logical((goodFit_f(bead, :))));
        
        template = zeros(numTemplates,templateSize,templateSize);
        for a=1:numTemplates
            for b = round((forwardStartFrames(a)+forwardEndFrames(a))/2):forwardEndFrames(a)
                % load data frames
                %                     for n = 1:numFiles  % it seems we do not need a loop
                %                     here, cause the value of b will determine which frame
                %                     is required to import. And we also need to use the
                %                     value of b as a if condition, this might mass up the
                %                     for loop.
                if b <= numFrames_temp
                    
                    
                    if isDcimg
                        [framedata,totalframes]= dcimgmatlab(b-1, dcimgfile);
                        framedatatrans = transpose (framedata);
                        dataFrame = double(framedatatrans)-darkAvg;
                    else
                        dataFileInfo = imfinfo(dataFile{1});
                        dataFrame = double(imread([dataFile{1}],b,'Info',dataFileInfo))-darkAvg;
                    end
                    data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
                        ROI(1):ROI(1)+ROI(3)-1);
                elseif b > numFrames_temp & b <= 2*numFrames_temp
                    dataFileInfo = imfinfo(dataFile{2});
                    dataFrame = double(imread([dataFile{2}],b-numFrames_temp,'Info',dataFileInfo))-darkAvg;
                    data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
                        ROI(1):ROI(1)+ROI(3)-1);
                elseif b > 2*numFrames_temp & b <= numFramesTotal
                    dataFileInfo = imfinfo(dataFile{3});
                    dataFrame = double(imread([dataFile{3}],b-2*numFrames_temp,'Info',dataFileInfo))-darkAvg;
                    data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
                        ROI(1):ROI(1)+ROI(3)-1);
                end
                
                % subtract the background and continue
                if useWaveSub
                    bkgndImg = f_waveletBackground(data);
                    data = data - bkgndImg;
                end
                
                good = PSFfits(:,1)==b & PSFfits(:,2)==bead;
                x_Pos = round(PSFfits(good,14)/nmPerPixel);
                y_Pos = round(PSFfits(good,15)/nmPerPixel);
                
                % check to make sure x_Pos,y_Pos are valid values
                if y_Pos-floor(templateSize/2) < 1
                    y_Pos = floor(templateSize/2)+1;
                end
                if y_Pos+floor(templateSize/2) > size(data,1)
                    y_Pos = size(data,1)-floor(templateSize/2);
                end
                if x_Pos-floor(templateSize/2) < 1
                    x_Pos = floor(templateSize/2)+1;
                end
                if x_Pos+floor(templateSize/2) > size(data,2)
                    x_Pos = size(data,2)-floor(templateSize/2);
                end
                DHPSF_Image = data(y_Pos-floor(templateSize/2):y_Pos+floor(templateSize/2),...
                    x_Pos-floor(templateSize/2):x_Pos+floor(templateSize/2));
                DHPSF_Image = DHPSF_Image(1:templateSize,1:templateSize);   % resize to the correct dimension if necessary
                template(a,:,:) = squeeze(template(a,:,:)) + DHPSF_Image;
                
            end
            %     imagesc(squeeze(template(a,:,:)))
            %     axis square
            %     drawnow
            %     pause(1)
            
        end
        
        save([outputFilePrefix 'bead ' num2str(bead) ' templates.mat'], 'template')
        if isDcimg
            clear data template framedatatrans dataFrame mex;
        else
            clear data template framedatatrans dataFrame;
        end
    end
else
    %% Generate template Stack of a chosen bead
    templateSize = 20;      % 2*round(10*160/nmPerPixel);  % 26;
    
    numFrames_temp = double(totalframes);
    if ~isDcimg
        dataFileInfo = imfinfo(dataFile{1});
        numFrames_temp =  length(dataFileInfo);
    end
    numTemplates = zeros(numBeads,1); templateVals = cell(numBeads,1);
    forwardStartFrames = cell(numBeads,1); forwardEndFrames = cell(numBeads,1);
    for bead = 1:numBeads
        numTemplates(bead,1) =length(startFrame(logical((goodFit_f(bead, :)))));
        forwardStartFrames{bead,1} = startFrame(logical((goodFit_f(bead, :))));
        forwardEndFrames{bead,1} = endFrame(logical((goodFit_f(bead, :))));
        templateVals{bead,1} = zeros(numTemplates(bead,1),templateSize,templateSize);
    end
    maxsteps = max(numTemplates);
    
    StartFrames = nan(maxsteps,numBeads);
    EndFrames = nan(maxsteps,numBeads);
    for bead = 1:numBeads
        StartFrames(1:numTemplates(bead,1),bead) = forwardStartFrames{bead,1};
        EndFrames(1:numTemplates(bead,1),bead) = forwardEndFrames{bead,1};
    end
    framestoprocessStart = unique(StartFrames(~isnan(StartFrames)));
    framestoprocessEnd = unique(EndFrames(~isnan(EndFrames)));
    templateImgs = cell(length(framestoprocessStart),1);
    parfor process = 1:length(framestoprocessStart)
        [FrameIdx{process},BeadIdx{process}]=ind2sub(size(StartFrames),find(framestoprocessStart(process)==StartFrames)');
        templateImgs{process,1} = zeros(length(BeadIdx{process}),templateSize,templateSize);
        for b = round((framestoprocessStart(process)+framestoprocessEnd(process))/2):framestoprocessEnd(process)
            data = zeros(ROI(4),ROI(3)); bkgndImg = data;
            dataFrame = zeros(size(gain));
            if b <= numFrames_temp
                if isDcimg
                    [framedata,~]= dcimgmatlab(b-1, dcimgfile);
                    framedatatrans = transpose (framedata);
                    dataFrame = double(framedatatrans)-darkAvg;
                else
                    dataFileInfo = imfinfo(dataFile{1});
                    dataFrame = double(imread([dataFile{1}],b,'Info',dataFileInfo))-darkAvg;
                end
                data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
                    ROI(1):ROI(1)+ROI(3)-1);
            elseif b > numFrames_temp && b <= 2*numFrames_temp
                dataFileInfo = imfinfo(dataFile{2});
                dataFrame = double(imread([dataFile{2}],b-numFrames_temp,'Info',dataFileInfo))-darkAvg;
                data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
                    ROI(1):ROI(1)+ROI(3)-1);
            elseif b > 2*numFrames_temp && b <= numFramesTotal
                dataFileInfo = imfinfo(dataFile{3});
                dataFrame = double(imread([dataFile{3}],b-2*numFrames_temp,'Info',dataFileInfo))-darkAvg;
                data = dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
                    ROI(1):ROI(1)+ROI(3)-1);
            end
            if useWaveSub
                bkgnd = f_waveletBackground(dataFrame(ROI(2):ROI(2)+ROI(4)-1, ...
                    ROI(1):ROI(1)+ROI(3)-1));
                data = data - bkgnd;
            end
            for c = 1:length(BeadIdx{process})
                good = PSFfits(:,1)==b & PSFfits(:,2)==BeadIdx{process}(c);
                x_Pos = round(PSFfits(good,14)/nmPerPixel);
                y_Pos = round(PSFfits(good,15)/nmPerPixel);
                
                % check to make sure x_Pos,y_Pos are valid values
                if y_Pos-floor(templateSize/2) < 1
                    y_Pos = floor(templateSize/2)+1;
                end
                if y_Pos+floor(templateSize/2) > size(data,1)
                    y_Pos = size(data,1)-floor(templateSize/2);
                end
                if x_Pos-floor(templateSize/2) < 1
                    x_Pos = floor(templateSize/2)+1;
                end
                if x_Pos+floor(templateSize/2) > size(data,2)
                    x_Pos = size(data,2)-floor(templateSize/2);
                end
                DHPSF_Image = data(y_Pos-floor(templateSize/2):y_Pos+floor(templateSize/2),...
                    x_Pos-floor(templateSize/2):x_Pos+floor(templateSize/2));
                DHPSF_Image = DHPSF_Image(1:templateSize,1:templateSize);   % resize to the correct dimension if necessary
                templateImgs{process,1}(c,:,:) = squeeze(templateImgs{process,1}(c,:,:)) + DHPSF_Image;
            end
        end
    end
    if isDcimg
        clear mex;
    end
    for b = 1:length(BeadIdx)
        for c = 1:length(BeadIdx{b})
            templateVals{BeadIdx{b}(c),1}(FrameIdx{b}(c),:,:) = templateImgs{b,1}(c,:,:);
        end
    end
    for bead = 1:numBeads
        template = templateVals{bead,1};
        save([outputFilePrefix 'bead ' num2str(bead) ' templates.mat'], 'template')
    end
    
    clear template templateImgs templateVals
    %%  Correction section
    maxddeg = 8;
    f_steps = 30:89;
    steps = 1:length(f_steps);
    diffAngles = diff(meanAngles(:,f_steps)');
    diffAngles = [-1*ones(numBeads,1),diffAngles'];
    meanAngles_NHA = meanAngles(:,f_steps);
    meanAngles_NHA(abs(diffAngles)>maxddeg) = nan;
    meanAngles_NHA(abs(meanAngles_NHA)>90) = nan;
    meanX_NHA = meanX(:,f_steps);
    meanY_NHA = meanY(:,f_steps);
    diffX = diff(meanX_NHA')';
    diffY = diff(meanY_NHA')';
    dTXY = 3.*max([nanmean(nanstd(diffX)),nanmean(nanstd(diffY))]);
    badXY = abs([zeros(numBeads,1),diffY])>dTXY|abs([zeros(numBeads,1),diffX])>dTXY;
    meanAngles_NHA(badXY) = nan;
    meanX_NHA(badXY) = nan;
    meanY_NHA(badXY) = nan;
    meanAngles_NHA = LocalElimination2(absLocs,meanAngles_NHA,outputFilePrefix);
    
    meanB = nan(numBeads,length(f_steps));
    for a = 1:numBeads
        meanA = meanAngles_NHA(a,:);
        if sum(~isnan(meanA))> (length(f_steps)/4)
            meanB(a,:) = interp1(steps(~isnan(meanA)),meanA(~isnan(meanA)),steps,'pchip',-90);
        end
    end
    meanAngles_NHA = meanAngles(:,f_steps);
    meanAngles_NHA(abs(meanAngles_NHA-meanB)>maxddeg) = nan;
    meanAngles_NHA = LocalElimination2(absLocs,meanAngles_NHA,outputFilePrefix);
    meanAngles_NHA(isnan(meanX_NHA) | isnan(meanY_NHA)) = nan;
    
    x_zero = nan(numBeads,1);
    y_zero = nan(numBeads,1);
    for a = 1:numBeads
        if sum(~isnan(meanAngles_NHA(a,:)))> (length(f_steps)/4)
            x_zero(a) = interp1(meanAngles_NHA(a,~isnan(meanAngles_NHA(a,:))),meanX_NHA(a,~isnan(meanAngles_NHA(a,:))),0,'linear');
            y_zero(a) = interp1(meanAngles_NHA(a,~isnan(meanAngles_NHA(a,:))),meanY_NHA(a,~isnan(meanAngles_NHA(a,:))),0,'linear');
        end
    end
    
    badBeads = isnan(x_zero) | isnan(y_zero) | (sum(isnan(meanAngles_NHA'))>round(length(f_steps)/3))';
    %     meanAngles_NHA = nan(numBeads,length(f_steps));
    %     meanAngles_NHA(abs(diffAngles)<maxddeg) = meanAngles(abs(diffAngles)<maxddeg);
    %     meanAngles_NHA = meanAngles_NHA(:,f_steps);
    %     meanAngles_NHA(bsxfun(@times,diffAngles,sign(nanmean(diffAngles)')')<-1) = nan;
    
    %     meanAngles_NHA = meanAngles;
    %     mc = 1;
    %     while mc ~= 0
    %         diffAngles = diff(meanAngles_NHA');
    %         diffAngles = [ones(numBeads,1),diffAngles'];
    %         ma= nanmedian(diffAngles);
    %         ms= median(nanstd(diffAngles));
    %         upbe = ma+3*ms;
    %         lpbe = ma-3*ms;
    %         mc = sum(sum(bsxfun(@gt,diffAngles,upbe)|bsxfun(@lt,diffAngles,lpbe)));
    %         meanAngles_NHA(bsxfun(@gt,diffAngles,upbe)|bsxfun(@lt,diffAngles,lpbe)) = nan;
    %     end
    %     mc = 1;
    %     while mc ~= 0
    %         ma= nanmedian(meanAngles_NHA);
    %         ms= nanstd(meanAngles_NHA);
    %         upbe = ma+3*ms;
    %         lpbe = ma-3*ms;
    %         mc = sum(sum(bsxfun(@gt,meanAngles_NHA,upbe)|bsxfun(@lt,meanAngles_NHA,lpbe)));
    %         meanAngles_NHA(bsxfun(@gt,meanAngles_NHA,upbe)|bsxfun(@lt,meanAngles_NHA,lpbe)) = nan;
    %     end
    %     badBeads = sum(isnan(meanAngles_NHA'))>round(length(f_steps)/4);
    meanAngles_NHA = meanAngles_NHA(~badBeads,:);
    meanX_NHA = meanX_NHA(~badBeads,:);
    meanY_NHA = meanY_NHA(~badBeads,:);
    x_NHA = absLocs(~badBeads,1); y_NHA = absLocs(~badBeads,2);
    numBeads_NHA = sum(~badBeads);
    
    %     x_NHA = absLocs(:,1); y_NHA = absLocs(:,2);
    %     numBeads_NHA = numBeads;
    
    meanCFocPos_good = meanCFocPos-meanCFocPos(60);
    meanCFocPos_NHA_rel_max = meanCFocPos_good(f_steps);
    %
    %     passFilt = true(1,numBeads_NHA);
    %     for q = 1:numBeads_NHA
    %         angles_NHA_good = meanAngles_NHA(q,f_steps);
    %         goodpoints = ~isnan(angles_NHA_good);
    %         SampleAngles = min(angles_NHA_good(goodpoints)):max(angles_NHA_good(goodpoints));
    %         SampleZ = interp1(angles_NHA_good(goodpoints),meanCFocPos_NHA_rel_max(goodpoints),SampleAngles,'spline');
    %         if max(SampleZ) > 1.05*max(meanCFocPos_NHA_rel_max) || min(SampleZ) < 1.05*min(meanCFocPos_NHA_rel_max)
    %             passFilt(q) = false;
    %         end
    %     end
    %
    %         meanAngles_NHA = meanAngles_NHA(:,f_steps);
    %     meanX_NHA = meanX_NHA(:,f_steps);
    %     meanY_NHA = meanY_NHA(:,f_steps);
    %     x_NHA = x_NHA(passFilt); y_NHA = y_NHA(passFilt);
    %     numBeads_NHA = sum(passFilt);
    for a = 1:numBeads_NHA
        xAngleZero = interp1(meanAngles_NHA(a,~isnan(meanAngles_NHA(a,:))),meanX_NHA(a,~isnan(meanAngles_NHA(a,:))),0,'pchip');
        yAngleZero = interp1(meanAngles_NHA(a,~isnan(meanAngles_NHA(a,:))),meanY_NHA(a,~isnan(meanAngles_NHA(a,:))),0,'pchip');
        meanX_NHA(a,:) = meanX_NHA(a,:) - xAngleZero;
        meanY_NHA(a,:) = meanY_NHA(a,:) - yAngleZero;
    end
    
    
    astps = repmat(steps,numBeads_NHA,1);
    vstps = nan(size(astps));
    vstps(~isnan(meanAngles_NHA)) = astps(~isnan(meanAngles_NHA));
    L_steps = max(min(vstps')):min(max(vstps'));
    clear vstps astps;
    meanAngles_NHA_L = meanAngles_NHA(:,L_steps);
    meanX_NHA_L = meanX_NHA(:,L_steps);
    meanY_NHA_L = meanY_NHA(:,L_steps);
    
    %% B spline and generate fine points
    order = 20;
    numSteps_Corr = length(f_steps);
    xy_shift_adj = cell(numBeads_NHA,1);
    xy_shift_adj_L = cell(numBeads_NHA,1);
    if (numSteps_Corr < order)
        display([' !!! Error: Choose n >= order=',num2str(order),' !!!']);
        return;
    end
    T = linspace(0,1,numSteps_Corr-order+2);
    y = linspace(0,1,10000);
    
    %% Generate 10000 points for x/y shift
    parfor a = 1:numBeads_NHA
        xy_shift = [meanAngles_NHA(a,~isnan(meanAngles_NHA(a,:)))',meanX_NHA(a,~isnan(meanAngles_NHA(a,:)))',meanY_NHA(a,~isnan(meanAngles_NHA(a,:)))'];
        xy_x_shift = DEBOOR(T,xy_shift(:,[1,2]),y,order);
        xy_y_shift = DEBOOR(T,xy_shift(:,[1,3]),y,order);
        x_shift_zero = interp1(xy_x_shift(:,1),xy_x_shift(:,2),0,'spline');
        y_shift_zero = interp1(xy_y_shift(:,1),xy_y_shift(:,2),0,'spline');
        xy_shift_adj{a,1} = [xy_x_shift(:,1),xy_x_shift(:,2) - x_shift_zero,xy_y_shift(:,2) - y_shift_zero];
        xy_shift_adj_L{a,1} = xy_shift_adj{a,1}(xy_shift_adj{a,1}(:,1)>= min(meanAngles_NHA_L(a,:)) & xy_shift_adj{a,1}(:,1)<= max(meanAngles_NHA_L(a,:)),:);
    end
    
    Corr_Z = cell(numBeads_NHA,1);
    Corr_Z_L = cell(numBeads_NHA,1);
    parfor a = 1:numBeads_NHA
        iPts = ~isnan(meanAngles_NHA(a,:));
        %         temp = [meanAngles_NHA(a,f_steps)',meanX_NHA(a,f_steps)',meanY_NHA(a,f_steps)',meanCFocPos_NHA_rel_max'];
        temp = [meanAngles_NHA(a,:)',meanCFocPos_NHA_rel_max'];
        Corr_Z{a,1} = temp(iPts,:);
        Corr_Z_L{a,1} = Corr_Z{a,1}(Corr_Z{a,1}(:,1)>= min(meanAngles_NHA_L(a,:)) & Corr_Z{a,1}(:,1)<= max(meanAngles_NHA_L(a,:)),:);
    end
    
    %
    %
    %     order = 20;
    %     numSteps_Corr = length(L_steps);
    %     xy_shift_adj_L = cell(numBeads_NHA,1);
    %     if (numSteps_Corr < order)
    %         display([' !!! Error: Choose n >= order=',num2str(order),' !!!']);
    %         return;
    %     end
    %     T = linspace(0,1,numSteps_Corr-order+2);
    %     y = linspace(0,1,10000);
    %
    %     %% Generate 10000 points for x/y shift
    %     parfor a = 1:numBeads_NHA
    %         xy_shift = [meanAngles_NHA_L(a,~isnan(meanAngles_NHA_L(a,:)))',meanX_NHA_L(a,~isnan(meanAngles_NHA_L(a,:)))',meanY_NHA_L(a,~isnan(meanAngles_NHA_L(a,:)))'];
    %         xy_x_shift = DEBOOR(T,xy_shift(:,[1,2]),y,order);
    %         xy_y_shift = DEBOOR(T,xy_shift(:,[1,3]),y,order);
    %         x_shift_zero = interp1(xy_x_shift(:,1),xy_x_shift(:,2),0,'spline');
    %         y_shift_zero = interp1(xy_y_shift(:,1),xy_y_shift(:,2),0,'spline');
    %         xy_shift_adj_L{a,1} = [xy_x_shift(:,1),xy_x_shift(:,2) - x_shift_zero,xy_y_shift(:,2) - y_shift_zero];
    %     end
    %     Corr_Z_L = cell(numBeads_NHA,1);
    %     parfor a = 1:numBeads_NHA
    %         iPts = ~isnan(meanAngles_NHA_L(a,:));
    %         %         temp = [meanAngles_NHA(a,f_steps)',meanX_NHA(a,f_steps)',meanY_NHA(a,f_steps)',meanCFocPos_NHA_rel_max'];
    %         temp = [meanAngles_NHA_L(a,:)',meanCFocPos_NHA_rel_max(L_steps)'];
    %         Corr_Z_L{a,1} = temp(iPts,:);
    %     end
    %
    %     %% Triple Error-Checked for Maximum Coverage
    %     for a = 1:numBeads_NHA
    %         temp = diff(Corr_Z{a,1});
    %         temp2 = temp(:,1)./temp(:,2);
    %         while sum(abs(temp2)>(abs(mean(temp2))+2.5*std(temp2))) > 0
    %             Corr_Z{a,1}(abs(temp2)>abs(mean(temp2))+2.5*std(temp2),:) = [];
    %             temp = diff(Corr_Z{a,1});
    %             temp2 = temp(:,1)./temp(:,2);
    %         end
    %     end
    save([outputFilePrefix 'CorrectionCalibration.mat'],'x_NHA','y_NHA',...
        'Corr_Z','xy_shift_adj','Corr_Z_L','xy_shift_adj_L');
end


end
