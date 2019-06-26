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

function [outputFilePrefix] = f_trackFiducials(dataFile,dataPath,calFile,calBeadIdx,templateFile,templateFrames,peakThreshold,...
    darkFile,darkPath,boxRadius,channel, gaussianFilterSigma,minDistBetweenSMs,...
    lobeDistBounds,conversionGain,nmPerPixel,EMGain,templateLocs,sigmaBounds, fittingMethod, iROI)
% f_trackFiducials is a module in easy_dhpsf that extracts the position of
% one or more fiducial beads in an image stack. This position is then used
% to correct the fit localizations from that image stack.

dlg_title = 'background subtraction type';
prompt = { 'Use wavelet subtraction?' };
def = { '0' };
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);

% whether to use wavelet background subtraction, or subtract mean of image
% Off by default to minimize any possible change to shape of
% templates.
useWaveSub = logical(str2num(inputdialog{1}));

conversionFactor = conversionGain/EMGain;
ampRatioLimit = 0.7;            % sets maximum allowed amplitude ratio
sigmaRatioLimit = 0.3;          % sets maximum allowed sigma ratio

scrsz = get(0,'ScreenSize');
numSyncFrames = 25;
% Options for lsqnonlin
options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
dataFile_stringname = dataFile{1}((length(dataPath)+1):length(dataFile{1}));
outputFilePrefix = cell(1,length(dataFile_stringname));

% sets absolute frame number for relating data to sequence log
absFrameNum = 1;

for stack = 1:length(dataFile)
    dcimgfile = fullfile(dataPath, dataFile_stringname);
    dcimgfile = strrep(dcimgfile, '\', '\\');
    [framedata,totalframes]= dcimgmatlab(0, dcimgfile);
    framedatatrans = transpose (framedata);
    
dlg_title = 'Select Frame Range';
prompt = { 'Select Frame Range' };
def = { ['[1:' num2str(totalframes) ']']};
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
frames = str2num(inputdialog{1});

numFrames = length(frames);
    
    [imgHeight,imgWidth]=size(framedatatrans);
    
    if stack == 1
        
        
        if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
            
            templateInfo = imfinfo(templateFile);
            if templateInfo(1).Height ~= templateInfo(1).Width
                error('Template is not square');
            end
            templateSize = templateInfo(1).Height;
        else
            if ~exist(templateFile)
                [templateFile, templatePath] = uigetfile({'*.mat'},...
                    'templateFile not found! Please reselect path:');
                templateFile = [templatePath templateFile];
            end
            load(templateFile);
            templateSize = size(template,2);
        end
        
    end
    
    % saves in labeled directory if a channel is selected
    if channel == '0'
        outputFilePrefix{stack} = [dataFile{stack}(1:length(dataFile{stack})-6) filesep 'fiduciaries ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
    else
        outputFilePrefix{stack} = [dataFile{stack}(1:length(dataFile{stack})-6) filesep channel(1) ' fiduciaries ' ...
            datestr(now,'yyyymmdd HHMM') filesep];
    end
    mkdir(outputFilePrefix{stack});
    
    
    if stack == 1
        
        % Compute dark counts
        
        if ~isequal(darkFile,0)
            % Computes average of dark frames for background subtraction
            
            darkFile1 = [darkPath darkFile];
            dcimgfile_D = fullfile(darkPath, darkFile);
            dcimgfile_D = strrep(dcimgfile_D, '\', '\\');
            [framedata_D,totalframes_D]= dcimgmatlab(0, dcimgfile_D);
            framedatatrans_D = transpose (framedata_D);
            numDarkFrames = totalframes_D;
            [darkHeight, darkWidth] = size(framedatatrans_D);
            
            darkAvg = zeros(darkHeight,darkWidth);
            for dframe = 1:numDarkFrames
                [framedata_D,totalframes_D]= dcimgmatlab(dframe-1, dcimgfile_D);
                framedatatrans_D = transpose (framedata_D);
                darkAvg = darkAvg + double(framedatatrans_D);
            end
            darkAvg = darkAvg/double(numDarkFrames);
            if ~isequal(size(darkAvg),[imgHeight imgWidth])
                warning('Dark count image and data image stack are not the same size. Resizing dark count image...');
                darkAvg = imresize(darkAvg,[imgHeight imgWidth]);
            end
        else
            darkAvg = 0;
        end
        clear darkFileInfo;
        
        %% Compute average image
        dataAvg = zeros(imgHeight,imgWidth);
        for frame = 1:200
            % for frame = 1:(numFrames/10)
            dcimgfile = fullfile(dataPath, dataFile_stringname);
            dcimgfile = strrep(dcimgfile, '\', '\\');
            
            [framedata,totalframes]= dcimgmatlab(frame-1, dcimgfile);
            framedatatrans = transpose (framedata);
            dataAvg = dataAvg + double(framedatatrans) - darkAvg;
        end
        dataAvg = dataAvg/200;
        
        % definitions related to templates
        numTemplates = length(templateFrames);
        %% Pick region of interest for analysis
        
        hROI = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(dataAvg);axis image;colormap hot;
        %% adjust for red vs. green gain and var
        if channel == 'g' || channel == '0'
            load('gainGreen.mat');
            load('VarGreen.mat');
            gain =gain;
%             variance = ReadN_DI;
            variance = ReadN_DI.^2; 
            ROI = imrect(gca,[1 1 64 64]);
        elseif channel == 'r'
            load('gainRed.mat');
            load('VarRed.mat');
            gain =gain;
%             variance = ReadN_DI;
            variance = ReadN_DI.^2;
            ROI = imrect(gca,[243 243 64 64]);
        end
        if channel == '0'
            load('gainRed.mat');
            load('VarRed.mat');
            gain =gain;
            %             variance = ReadN_DI;
            variance = ReadN_DI.^2; % changed by Ting Yan, 12/27/2016
        end
        
        title({'Double-click to choose region of interest for PSF extraction' ...
            mat2str(ROI.getPosition)});
        addNewPositionCallback(ROI,@(p) title({'Double-click to choose region of interest for PSF extraction' ...
            ['[xmin ymin width height] = ' mat2str(p,3)]}));
        %       make sure rectangle stays within image bounds
        fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
        setPositionConstraintFcn(ROI,fcn);
        ROI = round(wait(ROI));
        
        % if odd, rounds down to nearest even number for width and height
        % (must be divisible by two for padding step)
        ROI(3) = 2*floor(ROI(3)/2);
        ROI(4) = 2*floor(ROI(4)/2);
        
        cropHeight = ROI(4);
        cropWidth = ROI(3);
        close(hROI)
        %% Ask user for molecule location(s)
        
        % Plots the avg image .tif
        hMLocs = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(dataAvg(ROI(2):ROI(2)+ROI(4)-1, ...
            ROI(1):ROI(1)+ROI(3)-1));axis image;colorbar;colormap hot;
        title('Use LMB to select fiducials. Hit enter to stop or use RMB to select the final fiducial.');
        hold on;
        % User will click on molecules in image and then text will be imposed on the
        % image corresponding to the molecule number.
        % pointstofit is a n by 2 array -- x y values for each molecule
        % hold on;
        moleLocs = [];
        n = 0;
        % Loop, picking up the points.
        but = 1;
        while but == 1
            [xi,yi,but] = ginput(1);
            if isempty(xi)
                break
            end
            n = n+1;
            text(xi,yi,num2str(n),'color','white','fontsize',13,'fontweight','bold');
            moleLocs(n,:) = round([xi yi]);
        end
        hold off;
        
        numMoles = size(moleLocs,1);
        
        
        % save the image with the numbered molecules for future reference
        saveas(hMLocs,[outputFilePrefix{stack} 'molecule map.png']);
        close(hMLocs)
        

        
        %% prepare template for template matching
        hMFit = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        
        % pad template to same size as input
        templatePad = zeros(numTemplates,cropHeight,cropWidth);
        templateFT = zeros(numTemplates,cropHeight,cropWidth);
        for a=1:numTemplates
            
            if strcmp(templateFile(length(templateFile)-2:length(templateFile)),'tif')
                templatePad(a,:,:) = padarray(squeeze(template(a,:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(a,:,:))));
            else
                templatePad(a,:,:) = padarray(squeeze(template(templateFrames(a),:,:)),...
                    [(cropHeight-size(template,2))/2 ...
                    (cropWidth-size(template,3))/2],min(min(template(templateFrames(a),:,:))));
            end
            
            % multiplying by conjugate of template in FT domain is equivalent
            % to flipping the template in the real domain
            templateFT(a,:,:) = conj(fft2(squeeze(templatePad(a,:,:))));
        end
        clear templatePad;
        
        % apply Gaussian filter to phase correlation data to weight low frequencies
        % more heavily since SNR is higher there
        gaussianFilter = abs(fft2(fspecial('gaussian', [cropHeight cropWidth], gaussianFilterSigma)));
        
        
    end
    
    %% Fit chosen molecules throughout entire image
    
    % [frameNum moleNum amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2
    %  bkgndMean totalFitError goodFit xCenter yCenter angle numPhotons ...
    %  xAberrCorrected yAberCorrected zLocation]
    PSFfits = zeros(numFrames*numMoles, 23);
    numPSFfits = 0;
    bkgndFits = zeros(numFrames,8);
    startTime = tic;

    selectedFrames=frames;
 
    for a=1:numFrames
        
        if ~logical(sum(frames(a)==selectedFrames))
            for b=1:numMoles
                rowIdx = (a-1)*numMoles+b;
                PSFfits(rowIdx,13) = -1100;     % this frame was not analyzed
                PSFfits(rowIdx,1:2) = [a b];
            end
            continue
        end
        %%divide by GAIN!! if MLE need to add var!
        
        dcimgfile = fullfile(dataPath, dataFile_stringname);
        dcimgfile = strrep(dcimgfile, '\', '\\');
        
        [framedata,~]= dcimgmatlab(frames(a)-1, dcimgfile);
        framedatatrans = transpose (framedata);
        
        data = double(framedatatrans)-darkAvg;
        % data = double(data./gain);
        data = double(data*conversionFactor);
         
         if ~isequal(size(gain),[imgHeight imgWidth])

             gain = gain(iROI(2):iROI(2)+imgHeight-1,iROI(1):iROI(1)+imgWidth-1);
             variance = variance(iROI(2):iROI(2)+imgHeight-1,iROI(1):iROI(1)+imgWidth-1);
         end
    
        if strcmp(fittingMethod, 'MLE with DG model')
            datavar = double((variance./gain)./gain);
            data = double(data + datavar);
            datavar = datavar(ROI(2):ROI(2)+ROI(4)-1, ...
                ROI(1):ROI(1)+ROI(3)-1);
        end
        data = data(ROI(2):ROI(2)+ROI(4)-1, ...
            ROI(1):ROI(1)+ROI(3)-1);        % crop image to ROI
        
        % subtract the background and continue
        if useWaveSub || ~exist('useWaveSub')
            bkgndImg = f_waveletBackground(data);
        else
            % subtract median to get rough idea of #photons of fiducial
            bkgndImg = repmat(median(data(:)),size(data));
        end
        if strcmp(fittingMethod, 'LSQ with DG model')
            data = data - bkgndImg;
        end
        
        
      
        
        %% do fitting to extract exact locations of DH-PSFs

        
        % create reconstructed DH-PSF image from fitted data

        reconstructImg = zeros(cropHeight, cropWidth);
        
        for b=1:numMoles
            
            rowIdx = (a-1)*numMoles+b;
            dataFT = fft2(data,cropHeight,cropWidth);
            maxPeakImg = zeros(cropHeight,cropWidth);
            % matrix PSFLocs stores information about double helices that were
            % found via template matching
            % rows are different matches
            % [xLocation yLocation matchingTemplateNumber matchConfidence];
            PSFLocs = zeros(100,4);
            numPSFLocs = 0;
            for c=1:numTemplates
                % try no prefiltering
                %H = 1;
                % try phase correlation
                %H = 1./(abs(dataFT).*abs(squeeze(templateFT(b,:,:))));
                % try weighted phase correlation (emphasizing low frequency
                % components
                H = gaussianFilter./(abs(dataFT).*abs(squeeze(templateFT(c,:,:))));
                
                % normalize H so it doesn't add any energy to template match
                %H = H / sqrt(sum(abs(H(:)).^2));
                
                peakImg = ifftshift(ifft2(dataFT.*squeeze(templateFT(c,:,:)).*H));
                
                % normalize response of peakImg by dividing by number of pixels in
                % data
                %peakImg = peakImg / (cropHeight*cropWidth);
                maxPeakImg = max(maxPeakImg, peakImg);
                
                %threshold = mean(peakImg(:))+peakThreshold*std(peakImg(:));
                % These thresholds are defined relative to the single
                % molecule thresholds.  This can become a problem, if the
                % fiducial starts to bleach during the experiment.
                peakImg(peakImg < 3*peakThreshold(stack,c)) = 3*peakThreshold(stack,c);
                
                if isreal(peakImg) && sum(sum(isnan(peakImg)))==0
                    temp = find(imregionalmax(peakImg));
                else
                    % Write log file
                    fileID = fopen('peakImg log.txt','a');
                    fprintf(fileID,[datestr(now) '\r\n' dataFile{stack}]);
                    fprintf(fileID,'\r\nROI: [%d %d %d %d]',ROI);
                    fprintf(fileID,'\r\nFrame: %d\r\npeakImg matrix:\r\n',a);
                    dlmwrite('peakImg log.txt',peakImg,'-append','delimiter','\t','newline','pc')
                    fprintf(fileID,'\r\n********************NEXT********************\r\n\r\n');
                    fclose('all');
                    
                    logFlag = logFlag + 1;
                    peakImg(isnan(peakImg)) = 0; %inserted to deal with NaNs -AC 6/22
                    temp = find(imregionalmax(real(peakImg)));
                end
                
                
                % make sure threshold didn't eliminate all peaks and create
                % lots of matches
                if length(temp) < cropHeight*cropWidth/2;
                    [tempY tempX] = ind2sub([cropHeight cropWidth],temp);
                    PSFLocs(numPSFLocs+(1:length(temp)),:) = ...
                        [tempX tempY c*ones(length(temp),1) peakImg(temp)];
                    numPSFLocs = numPSFLocs+length(temp);
                end
            end
            PSFLocs = PSFLocs(PSFLocs(:,1)~=0,:);
            clear H dataFT peakImg
            
            %% filter out extraneous matches due to very strong signals
            
            if numPSFLocs == 0
                % This error flag indicates that no good PSF candiates
                % were found by the template matching given the
                % thresholds
                PSFfits(rowIdx,13) = -1101;
                PSFfits(rowIdx,1:2) = [a b];
                continue;
                
            elseif numPSFLocs > 0
                % keep only the matches that are within the
                % coordinates a previous successfil bead fit
                temp = [];
                for c=1:size(PSFLocs,1)
                    if (moleLocs(b,1)-PSFLocs(c,1))^2 +(moleLocs(b,2)-PSFLocs(c,2))^2 <= (minDistBetweenSMs)^2
                        temp = [temp; PSFLocs(c,:)];
                    end
                end
                % sort location matrix in decending order of confidence
                if size(temp,2)>0
                    temp = sortrows(temp,-4);
                    % copy most confident match
                    PSFLocs = temp(1,:);
                    %                 PSFLocs(1,:) = temp(1,:);
                    numPSFLocs = 1;
                else
                    % This error flag indicates that the PSF candiates
                    % were found only too far away from the
                    % coordinates a previous successfil bead fit
                    PSFfits(rowIdx,13) = -1102;
                    PSFfits(rowIdx,1:2) = [a b];
                    continue
                end
   
            end
            
            if numPSFLocs > 0
                % create indices to use for fitting
                % To allow for sample drift increase the box size for
                % fitting by 2 pixels
                [xIdx yIdx] = meshgrid(moleLocs(b,1)-(boxRadius+2):moleLocs(b,1)+(boxRadius+2), ...
                    moleLocs(b,2)-(boxRadius+2):moleLocs(b,2)+(boxRadius+2));

                % compute initial parameters from the location of two spots in
                % the templates
                % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
                fitParam(3) = PSFLocs(1) + templateLocs(PSFLocs(3),1)-(templateSize/2+0.5);
                fitParam(4) = PSFLocs(2) + templateLocs(PSFLocs(3),2)-(templateSize/2+0.5);
                fitParam(5) = PSFLocs(1) + templateLocs(PSFLocs(3),3)-(templateSize/2+0.5);
                fitParam(6) = PSFLocs(2) + templateLocs(PSFLocs(3),4)-(templateSize/2+0.5);
               % make sure initial guess lies within the ROI: if not, move on
                if fitParam(3)<min(xIdx(:)) || fitParam(5)<min(xIdx(:)) ...
                        || fitParam(3)>max(xIdx(:)) || fitParam(5)>max(xIdx(:)) ...
                        || fitParam(4)<min(yIdx(:)) || fitParam(6)<min(yIdx(:)) ...
                        || fitParam(4)>max(yIdx(:)) || fitParam(6)>max(yIdx(:))
                    PSFfits(rowIdx,13) = -1103;
                    PSFfits(rowIdx,1:2) = [a b];
                    continue;
                end
                fitParam(1) = data(round(fitParam(4)),round(fitParam(3)));
                fitParam(2) = data(round(fitParam(6)),round(fitParam(5)));
                fitParam(7) = 1.8;
                fitParam(8) = 1.8;
                switch (fittingMethod)
                    case 'LSQ with DG model'
                        lowerBound = [0 0 min(xIdx(:)) min(yIdx(:)) min(xIdx(:)) min(yIdx(:)) ...
                            sigmaBounds(1) sigmaBounds(1)];
                        upperBound = [max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                            max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                            max(xIdx(:)) max(yIdx(:)) max(xIdx(:)) max(yIdx(:)) ...
                            sigmaBounds(2) sigmaBounds(2)];
                        
                        %% Fit with lsqnonlin
                        [fitParam,temp,residual,exitflag] = lsqnonlin(@(x) ...
                            f_doubleGaussianVector(x,data(yIdx(:,1),xIdx(1,:)),0,xIdx,yIdx),...
                            fitParam,lowerBound,upperBound,options);
                        PSFfits(rowIdx,1:13) = [a b fitParam 0 sum(abs(residual)) exitflag];
                        
                    case'MLE with DG model'
                        fitParam(9) = 0;
                        bkgndImg2 = median(median(bkgndImg));
                        %                     bkgndImg2 = bkgndImg;
                        
                        lowerBound = [0 0 min(xIdx(:)) min(yIdx(:)) min(xIdx(:)) min(yIdx(:)) ...
                            sigmaBounds(1) sigmaBounds(1) -bkgndImg2];
                        upperBound = [max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                            max(max(data(yIdx(:,1),xIdx(1,:)))) ...
                            max(xIdx(:)) max(yIdx(:)) max(xIdx(:)) max(yIdx(:)) ...
                            sigmaBounds(2) sigmaBounds(2) bkgndImg2];
                        arr = [0, 0, 0, 0, 0, 0, 0, 0, 0];
                        arr2 = 10;
                        Aeq = [];
                        beq =[];
                        nonlcon = [];
                        
                       
                        hMLE= @(fitParam)Likelihood( fitParam, data(yIdx(:,1),xIdx(1,:)), bkgndImg2,  datavar(yIdx(:,1),xIdx(1,:)), xIdx,yIdx);
                        
                        options = optimset('Display','off');
                        [fitParam, fval, exitflag] = fmincon(hMLE, fitParam, arr, arr2, Aeq, beq, lowerBound, upperBound,  nonlcon, options);
                        PSFfits(rowIdx,1:13) = [a b fitParam  NaN exitflag];
                        
                end
                % shift coordinates relative to entire dataset (not just ROI)
                PSFfits(rowIdx,5) = PSFfits(rowIdx,5) + ROI(1)-1;
                PSFfits(rowIdx,6) = PSFfits(rowIdx,6) + ROI(2)-1;
                PSFfits(rowIdx,7) = PSFfits(rowIdx,7) + ROI(1)-1;
                PSFfits(rowIdx,8) = PSFfits(rowIdx,8) + ROI(2)-1;
                % Calculate midpoint between two Gaussian spots
                % shift coordinates relative to entire dataset (not just ROI) and
                % convert from pixels to nm
                PSFfits(rowIdx,14) = ((fitParam(3)+fitParam(5))/2 + ROI(1)-1)*nmPerPixel;
                PSFfits(rowIdx,15) = ((fitParam(4)+fitParam(6))/2 + ROI(2)-1)*nmPerPixel;
               
                
                % Below is the calculation of the angle of the two lobes.
                % Remember that two vertical lobes is focal plane because camera
                % outputs data that is rotated. Therefore, we want y2>y1 for all
                % angle calculations (so that -90<=angle<=90, and we use swap
                % the use of x and y for the atan2 calculation.
                x1 = fitParam(3);
                x2 = fitParam(5);
                y1 = fitParam(4);
                y2 = fitParam(6);
                % swap if y1>y2
                if (y1 > y2)
                    tx = x1; ty = y1;
                    x1 = x2; y1 = y2;
                    x2 = tx; y2 = ty;
                    clear tx ty;
                end
                %Finds the angle
                PSFfits(rowIdx,16) = atan2(-(x2-x1),y2-y1) * 180/pi;
                clear x1 x2 y1 y2;
                
                %Below is a way to count the photons. It integrates the box and
                %subtracts the boxarea*offset from the fit. It is inherently flawed
                %if there happens to be bright pixels inside of the fitting region.
                totalCounts = sum(sum(data(yIdx(:,1),xIdx(1,:))));
                %Do not multiply by conversion factor: if LSQ we are
                %already in units of photons, if MLE we are already in
                %units of lambda
                %PSFfits(rowIdx,17) = totalCounts*conversionFactor;
                PSFfits(rowIdx,17) = totalCounts;
                
                %The interlobe distance
                lobeDist = sqrt((fitParam(3)-fitParam(5)).^2 + ...
                    (fitParam(4)-fitParam(6)).^2);
                PSFfits(rowIdx,18) = lobeDist;
                
                %Amplitude Ratio
                ampRatio = abs(fitParam(1) - fitParam(2))/sum(fitParam(1:2));
                PSFfits(b,19) = ampRatio;
                
                % Gaussian width Ratio
                simgaRatio = abs(fitParam(7) - fitParam(8))/sum(fitParam(7:8));
                PSFfits(b,20) = simgaRatio;
                
                %% Now evaluate the fits
                % Conditions for fits (play with these):
                % (1) Amplitude of both lobes > 0
                % (2) All locations x1,y1, x2,y2 lie inside area of small box
                % (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
                % (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
                % (5) Make sure amplitudes are within 100% of one another
                % (6) Make sure totalFitError/(total number of photons) < 1.05
                
                if exitflag > 0
                    if fitParam(1)<0 || fitParam(2)<0
                        PSFfits(rowIdx,13) = -1001;
                    end
                    if fitParam(3)<min(xIdx(:)) || fitParam(5)<min(xIdx(:)) ...
                            || fitParam(3)>max(xIdx(:)) || fitParam(5)>max(xIdx(:)) ...
                            || fitParam(4)<min(yIdx(:)) || fitParam(6)<min(yIdx(:)) ...
                            || fitParam(4)>max(yIdx(:)) || fitParam(6)>max(yIdx(:))
                        PSFfits(rowIdx,13) = -1002;
                    end
                    if fitParam(7)<=sigmaBounds(1) || fitParam(8)<=sigmaBounds(1) ...
                            || fitParam(7)>=sigmaBounds(2) || fitParam(8)>=sigmaBounds(2)
                        PSFfits(rowIdx,13) = -1003;
                    end
                    if simgaRatio > sigmaRatioLimit;
                        PSFfits(rowIdx,13) = -1004;
                    end
                    if lobeDist < lobeDistBounds(1) || lobeDist > lobeDistBounds(2)
                        PSFfits(rowIdx,13) = -1005;
                    end
                    if ampRatio > ampRatioLimit;
                        PSFfits(rowIdx,13) = -1006;
                    end

                    
                end
                
                % if fit was successful, use the computed center location as center
                % of box for next iteration
                if PSFfits(rowIdx,13) > 0
                    
                    
                    % find previously good fit
                    index = find(PSFfits(:,13)>0&PSFfits(:,2)==b);
                    % Alternatively, could average a couple of previous
                    % frames to make distance threshold adjustable on the
                    % fly
                    if ~isempty(index) && length(index)>1
                        
                        index = index(end-1);
                        
                        prevPos = [(PSFfits(index,5)+PSFfits(index,7))/2,...
                            (PSFfits(index,6)+PSFfits(index,8))/2];
                        currPos = [(PSFfits(rowIdx,5)+PSFfits(rowIdx,7))/2,...
                            (PSFfits(rowIdx,6)+PSFfits(rowIdx,8))/2];
                        euclidDist = sqrt(sum((prevPos - currPos).^2))*nmPerPixel;
                        
                        if euclidDist < 400  % maximum allowable XY drift per frame
                            
                            %update bead position
                            moleLocs(b,:) = [round((fitParam(3)+fitParam(5))/2), ...
                                round((fitParam(4)+fitParam(6))/2)];
                            
                            % plot image reconstruction so that fits can be checked
                            [xIdx yIdx] = meshgrid(1:cropWidth,1:cropHeight);
                            reconstructImg = reconstructImg + ...
                                fitParam(1).*exp( -((xIdx-fitParam(3)).^2+(yIdx-fitParam(4)).^2.) / (2.*fitParam(7).^2)) ...
                                +fitParam(2).*exp( -((xIdx-fitParam(5)).^2+(yIdx-fitParam(6)).^2.) / (2.*fitParam(8).^2));
                            
                        else
                            PSFfits(rowIdx,13) = -1104;     % bead drifted too much
                        end
                        
                    else
                        %update bead position witout backcheck
                        moleLocs(b,:) = [round((fitParam(3)+fitParam(5))/2), ...
                            round((fitParam(4)+fitParam(6))/2)];
                        
                        % plot image reconstruction so that fits can be checked
                        [xIdx yIdx] = meshgrid(1:cropWidth,1:cropHeight);
                        reconstructImg = reconstructImg + ...
                            fitParam(1).*exp( -((xIdx-fitParam(3)).^2+(yIdx-fitParam(4)).^2.) / (2.*fitParam(7).^2)) ...
                            +fitParam(2).*exp( -((xIdx-fitParam(5)).^2+(yIdx-fitParam(6)).^2.) / (2.*fitParam(8).^2));
                    end
                    
                end
                
                
            end
            PSFfits(rowIdx,1:2) = [a b];
        end
        
        %%  plot results of template matching and fitting
        if logical(sum(frames(a)==selectedFrames))
            %             hMFit = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
            set(0,'CurrentFigure',hMFit);
            
            subplot('Position',[0.025 0.025 .9/2 .95]);
            imagesc(data);axis image;colormap hot;
            title(['Frame ' num2str(a) ': raw data - dark counts']);
            
            
            subplot('Position',[0.525 0.025 .9/2 .95]);
            imagesc(reconstructImg,[min(data(:)) max(data(:))]);axis image;
            title('Image reconstructed from fitted matches');
            drawnow;
        end
        
        clear mex;
        
        
    end
    elapsedTime = toc(startTime);
    clear data dataFileInfo residual dataAvg reconstructImg xIdx yIdx temp;
    fps = numFrames/elapsedTime
    moleculesPerSec = numFrames*numMoles/elapsedTime
    %     close(hMFit)
    %% Translate angle into corrected x,y positions and z position
    
   if ~exist(calFile)
        [calFile, calPath] = uigetfile({'*.mat'},...
            'calFile not found! Please reselect path:');
        calFile = [calPath calFile];
    end
    load(calFile);%squeeze(meanAngles(1,calBeadIdx,goodFit_forward)
    goodFit_forward = logical(goodFit_f(calBeadIdx,:));
    PSFfits(:,21) = PSFfits(:,14) - interp1((meanAngles(calBeadIdx,goodFit_forward)),(meanX(calBeadIdx,goodFit_forward)),PSFfits(:,16),'spline');
    PSFfits(:,22) = PSFfits(:,15) - interp1((meanAngles(calBeadIdx,goodFit_forward)),(meanY(calBeadIdx,goodFit_forward)),PSFfits(:,16),'spline');
    PSFfits(:,23) = interp1((meanAngles(calBeadIdx,goodFit_forward)),(z(goodFit_forward)),PSFfits(:,16),'spline');
    %% Output raw fit data
    
    %     textHeader = {'frame number' 'molecule number' ...
    %         'amp 1 (counts)' 'amp 2 (counts)' 'x mean 1 (px)' 'y mean 1 (px)' ...
    %         'x mean 2 (px)' 'y mean 2 (px)' 'sigma 1 (px)' 'sigma 2 (px)' 'background mean (counts)' ...
    %         'total fit error (counts)' 'good fit flag' 'x center (nm)' 'y center (nm)' ...
    %         'angle (deg)' 'number of photons' 'interlobe distance' 'amplitude ratio' 'sigma ratio' 'aberration corrected x location (nm)' ...
    %         'aberration corrected y location (nm)' 'z location (nm)'};
    % save fit info to MATLAB mat file
    save([outputFilePrefix{stack} 'raw fits.mat']);
    
    
    %% compute movement of fiduciaries
    
    devX = zeros(numFrames,numMoles);
    devY = zeros(numFrames,numMoles);
    devZ = zeros(numFrames,numMoles);
    goodFitFlag = zeros(numFrames,numMoles);
    numPhotons = zeros(numFrames,numMoles);
    avgDevX = zeros(numFrames,1);
    avgDevY = zeros(numFrames,1);
    avgDevZ = zeros(numFrames,1);
    numValidFits = zeros(numFrames,1);
    
    syncFrames = find(PSFfits(:,13)>0);
    syncFrames = syncFrames(end-(numSyncFrames-1):end);
    syncFrames = syncFrames';
    
    
    for molecule = 1:numMoles
        %% extract fitting parameters for this molecule
        moleculeFitParam = PSFfits(PSFfits(:,2) == molecule, :);
        
        goodFitFlag(:,molecule) = moleculeFitParam(:,13);
        goodFit = goodFitFlag(:,molecule) > 0;
        goodFitX = goodFit;
        goodFitY = goodFit;
        goodFitZ = goodFit;
        
        % compute deviation with respect to bead location averaged over last
        % numSyncFrames frames of the movie
        devX(:,molecule) = moleculeFitParam(:,21) ...
            - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),21));
        devY(:,molecule) = moleculeFitParam(:,22) ...
            - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),22));
        devZ(:,molecule) = moleculeFitParam(:,23) ...
            - mean(moleculeFitParam(any(bsxfun(@eq,moleculeFitParam(:,1), syncFrames),2),23));
        numPhotons(:,molecule) = moleculeFitParam(:,17);
        

        
        hDevs = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        subplot(2,1,1);
        plot(moleculeFitParam((goodFit & goodFitX),1),devX((goodFit & goodFitX),molecule),'r');
        hold on;
        plot(moleculeFitParam((goodFit & goodFitY),1),devY((goodFit & goodFitY),molecule),'b');
        axis tight;
        legend('x','y');
        title(['Fiduciary ' num2str(molecule)]);
        xlabel('Frame #');
        ylabel('Position (nm)');
        
        subplot(2,1,2);
        plot(moleculeFitParam((goodFit & goodFitZ),1),devZ((goodFit & goodFitZ),molecule),'k');
        axis tight;
        legend('z');
        title(['Fiduciary ' num2str(molecule)]);
        xlabel('Frame #');
        ylabel('Position (nm)');
        print(hDevs,'-djpeg',[outputFilePrefix{stack} 'fiduciary ' num2str(molecule) ' stats.jpg']);
        
  
        % if particle was fit successfully, add its movement to the average
        avgDevX = avgDevX + (goodFit & goodFitX).*devX(:,molecule);
        avgDevY = avgDevY + (goodFit & goodFitY).*devY(:,molecule);
        avgDevZ = avgDevZ + (goodFit & goodFitZ).*devZ(:,molecule);
        numValidFits = numValidFits + (goodFit & goodFitX & goodFitY & goodFitZ);
    end
    
    %% compute average movement of all fiduciaries
    avgDevX = avgDevX./numValidFits;
    avgDevY = avgDevY./numValidFits;
    avgDevZ = avgDevZ./numValidFits;
    
  
    %% output fiduciary correction data
    save([outputFilePrefix{stack} 'fiduciary corrections.mat'],'avgDevX','avgDevY',...
        'avgDevZ','numValidFits','devX','devY','devZ','goodFitFlag','numPhotons');
    

    
    
end

end
