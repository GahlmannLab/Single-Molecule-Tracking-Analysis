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

% other approach to passing variables to this function
% function f_processFits(totalPSFfits,numFrames,ROI,conversionFactor,...
%     sigmaBounds,lobeDistBounds,ampRatioLimit,sigmaRatioLimit,nmPerPixel)
% current approach: load all relevant variables from f_fitSMs output
% directly
function f_processFits(catPSFfits,numFrames,fitFilePrefix, fidTrackX, fidTrackY, fidTrackZ, nmPerPixel,spatialCorr,useCurrent,currFidIdx)

currDir = pwd;
if ~exist(fitFilePrefix{1})
    fitFilePrefix{1}(1,1) = currDir(1);
end

plotSigNoise = 0;
plotPhotons = 0;
useTimeColors = 0;
plotAsTracks = 0;
plotOnTimes = 0;
plotClusters = 0;
numPhotonRange = [300 10000];
xyPrecRange = [0 150];
zPrecRange = [0 200];

numFramesAll = sum(numFrames);
if ~exist(fitFilePrefix{1},'dir')
    [xPath] = uigetdir(...
        'File not found! Please reselect path:');
    fitFilePrefix{1} = [xPath];
end
load([fitFilePrefix{1} filesep 'molecule fits.mat']);

threshVals = peakThreshold(1,:);

scatterSize = 30; %str2double(inputdialog{2});
wlShiftX = 133 * nmPerPixel; %str2double(inputdialog{3});
wlShiftY = 393 * nmPerPixel; %str2double(inputdialog{4});
wlnmPerPixel = 110;
if sum(isnan(catPSFfits(:,30))) == size(catPSFfits,1)
    warning(['No fiducial correction applied!']);
end
matchdist = 1000;
nSample = 1.33;         % index of refraction of sample
nOil = 1.518;           % index of immersion oil

scrsz = get(0,'ScreenSize');

%% open datafiles

[whiteLightFile whiteLightPath] = uigetfile({'*.tiff;*.tif;*.avi'},'Open image stack with white light image',dataPath);
if whiteLightFile ~= 0
    whiteLightFile = [whiteLightPath whiteLightFile];
    [whiteLightPath, whiteLightFile, whiteLightExt] = fileparts(whiteLightFile);
    whiteLightFile = [whiteLightFile whiteLightExt];
    if strcmp(whiteLightExt,'.avi')
        path(path,whiteLightPath);
        whiteLightInfo = VideoReader([whiteLightFile]);
        whiteLight = zeros(whiteLightInfo(1).Height, whiteLightInfo(1).Width);
        while hasFrame(whiteLightInfo)
            whiteLight = whiteLight + double(readFrame(whiteLightInfo));
        end
    elseif strcmp(whiteLightExt,'.tif') || strcmp(whiteLightExt,'.tiff')
        whiteLightFile = [whiteLightPath '\' whiteLightFile];
        whiteLightInfo = imfinfo(whiteLightFile);
        whiteLight = zeros(whiteLightInfo(1).Height, whiteLightInfo(1).Width);
        % average white light images together to get better SNR
        for a = 1:length(whiteLightInfo)
            whiteLight = whiteLight + double(imread(whiteLightFile, ...
                'Info', whiteLightInfo));
        end
    end
    whiteLight = imresize(whiteLight,wlnmPerPixel/nmPerPixel);
    WLHeight = size(whiteLight,1);
    WLWidth = size(whiteLight,2);
    whiteLight = rot90(whiteLight); %Phase detector is rotated 90 degrees counterclockwise
    whiteLight = (whiteLight-min(whiteLight(:)))/(max(whiteLight(:))-min(whiteLight(:)));
    if channel == 'r'
        wlShiftX = 133 * nmPerPixel; 
        wlShiftY = 393 * nmPerPixel; 
    elseif channel == 'g'
        wlShiftX = 150 * nmPerPixel; 
        wlShiftY = 179 * nmPerPixel; 
        whiteLight = flipud(whiteLight);
    end
end

ROI_initial = ROI;


%% Laser Intensity Profile
%Ask if user wants to use plot laser profile
dlg_title = 'Laser Profile';
prompt = {'Do you want to plot the laser profile?'};
def = {num2str(0)};
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
findLaserInt = str2num(inputdialog{1});
if findLaserInt == 1;
    
           temp = inputdlg({'What was the laser power at the objective? (in mW)'},...
                'Input laser power',...
                1,...
                {'0'});
            powerAtObjective = str2double(temp{1})/1000;
            
    
      [darkFileLaser, darkPathLaser] = uigetfile({'*.dcimg'; '*.tif'},...
                'Open Images with Laser Background Dark Counts',...
                'MultiSelect', 'on');

            if isDcimg
                dcimgfile_D_Laser = fullfile(darkPathLaser, darkFileLaser);
                
                dcimgfile_D_Laser = strrep(dcimgfile_D_Laser, '\', '\\');
                
                [framedata_D_Laser,totalframes_D_Laser]= dcimgmatlab(0, dcimgfile_D_Laser);
                framedatatrans_D_Laser = transpose(framedata_D_Laser);
                numDarkFrames_Laser = totalframes_D_Laser;
                
                [darkHeightLaser, darkWidthLaser] = size(framedatatrans_D_Laser);
            else
                numDarkFrames_Laser = length(darkFileLaserInfo);
                darkHeightLaser = darkFileLaserInfo(1).Height;
                darkWidthLaser = darkFileLaserInfo(1).Width;
            end
            darkAvgLaser = zeros(darkHeightLaser, darkWidthLaser);
            for frame = 1:numDarkFrames_Laser
                if isDcimg
                    [framedata_D_Laser,totalframes_D_Laser]= dcimgmatlab(frame-1, dcimgfile_D_Laser);
                    framedatatrans_D_Laser = transpose (framedata_D_Laser);
                    darkAvgLaser = darkAvgLaser + double(framedatatrans_D_Laser);
                else
                    darkAvgLaser = darkAvgLaser + double(imread([darkPathLaser darkFileLaser],frame,'Info', darkFileLaserInfo));
                end
            end
            darkAvgLaser = darkAvgLaser/double(numDarkFrames_Laser);
    
      [dataFileLaser, dataPathLaser] = uigetfile({'*.dcimg'; '*.tif'},...
                'Open Images with Laser Background',...
                'MultiSelect', 'on');
            
        dcimgfile = fullfile(dataPathLaser, dataFileLaser);
        dcimgfile = strrep(dcimgfile, '\', '\\');
        [framedataLaser,totalframesLaser]= dcimgmatlab(0, dcimgfile);
        totalframesLaser= double(totalframesLaser);
        framedatatransLaser = transpose(framedataLaser);
        
        [imgHeightLaser, imgWidthLaser] = size(framedatatransLaser);

            
            avgImgLaser = zeros(imgHeightLaser,imgWidthLaser);
            avgImgFramesLaser = min(200,length(frames));
            for a = 1:avgImgFramesLaser
                if isDcimg
                    dcimgfile = fullfile(dataPathLaser, dataFileLaser);
                    dcimgfile = strrep(dcimgfile, '\', '\\');
                    
                    [framedataLaser,totalframesLaser]= dcimgmatlab(a-1, dcimgfile);
                    totalframesLaser = double(totalframesLaser);
                    framedatatransLaser = transpose (framedataLaser);
                    laserBkgnd =  f_waveletBackground(framedatatransLaser);
                   

                    avgImgLaser = avgImgLaser + laserBkgnd - darkAvgLaser;  
                else
                    dcimgfile = fullfile(dataPath, dataFile);
                    avgImgLaser = avgImgLaser + double(imread([dataFile{stack}],frames(a),'Info', fileInfo)) - darkAvgLaser;
                end
            end
            avgImgLaser = avgImgLaser/avgImgFramesLaser;
                       
            if ~isequal(size(gain),[length(avgImgLaser(1,:)), length(avgImgLaser(:,1))])
                if channel == 'g'
                    load('gainGreen.mat', 'gain');
                elseif channel == 'r'
                    load('gainRed.mat', 'gain');
                else
                    load('gainRed.mat', 'gain');
                end
            end
            
            %This is measured in counts
              avgImgLaser = avgImgLaser./gain;
           

            
            %3D laser profile
            laserProfile=avgImgLaser*powerAtObjective/sum(sum(avgImgLaser))...
                /(nmPerPixel^2)*10^14;
            
%             %2D laser profile                           
%             imagesc(laserProfile);
%             axis square
%             colorbar
%             title('2D Laser Intensity Profile');
            
%             smooth 3D laser profile
            laserFilter = fspecial('disk',21);
            laserProfileSmooth = imfilter(laserProfile,laserFilter);
 
            figure
            surf(laserProfileSmooth)
            shading interp
            title('3D Laser Intensity Profile');
            zlabel('Intensity, W/cm^2');
            
end

%% Plot reconstructions

pass = 1;
anotherpass = true;

numFramesAll = sum(numFrames);
zRange = [-2000 2000];
frameRange = [1 numFramesAll];
fitErrorRange = [0 4];
%ask for exposure time for on time histogram later
dlg_title = 'Exposure time?';
         prompt = 'Exposure time in ms?';
         num_lines =1;
         def = {num2str(50)};
         inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
         exposureTime = str2double(inputdialog{1});

while anotherpass == true
    
    
    %% Plot the filter parameters
    
    if pass == 1
        
        fitErrorCol = 16;
        goodFitFlagCol = 17;
        numPhotonCol = 21;
        lobeDistCol = 22;
        ampRatioCol = 23;
        sigmaRatioCol = 24;
        templateNumCol = 5;
        templateStrCol = 6;
        
        initGoodFits = catPSFfits(:,goodFitFlagCol) > 0;
        
        figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        subplot(2,2,1)
        [n,xout] = hist(catPSFfits(:,lobeDistCol), 4:0.1:12);
        bar(xout,n)
        hold on
        [n,xout] = hist(catPSFfits(initGoodFits,lobeDistCol), 4:0.1:12);
        bar(xout,n, 'green')
        %         [n,xout] = hist(catPSFfits(~initGoodFits,lobeDistCol), 4:0.1:12);
        %         bar(xout,n, 'red')
        hold off
        xlabel('pixel'); ylabel('Frequency');
        title('Lobe Distance');
        legend('unfiltered','initially included')
        xlim([4 14]);
        
        subplot(2,2,2)
        unfilteredFitError = catPSFfits(:,fitErrorCol)*conversionFactor./catPSFfits(:,numPhotonCol);
        fitError = catPSFfits(initGoodFits,fitErrorCol)*conversionFactor./catPSFfits(initGoodFits,numPhotonCol);
        badFitError = catPSFfits(~initGoodFits,fitErrorCol)*conversionFactor./catPSFfits(~initGoodFits,numPhotonCol);
        [n,xout] = hist(unfilteredFitError(unfilteredFitError > 0 & unfilteredFitError < 20), 0:0.1:20);
        bar(xout,n)
        hold on
        [n,xout] = hist(fitError(fitError > 0 & fitError < 30), 0:0.1:30);
        bar(xout,n, 'green')
        hold off
        xlabel('Fit Error'); ylabel('Frequency');
        title('Fit Error');
        legend('unfiltered','initially included')
        xlim([0 30]);
        
        subplot(2,2,3)
        hist(catPSFfits(initGoodFits,ampRatioCol), 100)
        [n,xout] = hist(catPSFfits(:,ampRatioCol), 0:0.01:1);
        bar(xout,n)
        hold on
        [n,xout] = hist(catPSFfits(initGoodFits,ampRatioCol), 0:0.01:1);
        bar(xout,n, 'green')
        hold off
        xlabel('Amplitude Ratio'); ylabel('Frequency');
        title('Amplitude Ratio');
        xlim([-0.1 1]);
        legend('unfiltered','initially included')
        
        subplot(2,2,4)
        hist(catPSFfits(initGoodFits,sigmaRatioCol), 100)
        xlabel('Sigma Ratio'); ylabel('Frequency');
        title('Sigma Ratio');
        xlim([-0.1 1]);
        
    end
    %% Chose a desired parameter set for reconstruction
    
    dlg_title = 'Please Input Parameters';
    prompt = {  'Size of points in reconstruction',...
        'Temporal Color Coding',...
        'White light shift X (in nm)',...
        'White light shift Y (in nm)',...
        'Z range lower bound(in nm)',...
        'Z range upper bound(in nm)',...
        'First frame',...
        'Last frame'...
        'Plot as Tracks'...
        'Tracks: max match threshold',...
        'Plot S/N Ratio?'...
        'Plot numPhotons?'...
        'Plot Clusters?'...
        'Plot On Times?'
        };
    def = { ...
        num2str(scatterSize), ...
        num2str(useTimeColors), ...
        num2str(wlShiftX), ...
        num2str(wlShiftY), ...
        num2str(zRange(1)), ...
        num2str(zRange(2)), ...
        num2str(frameRange(1)), ...
        num2str(frameRange(2)), ...
        num2str(plotAsTracks), ...
        num2str(matchdist), ...
        num2str(plotSigNoise)...
        num2str(plotPhotons)...
        num2str(plotClusters)...
        num2str(plotOnTimes)
        };
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    
    scatterSize = str2double(inputdialog{1});
    useTimeColors = str2double(inputdialog{2});
    wlShiftX = str2double(inputdialog{3});
    wlShiftY = str2double(inputdialog{4});
    zRange = [str2double(inputdialog{5}) str2double(inputdialog{6})];
    frameRange = [str2double(inputdialog{7}) str2double(inputdialog{8})];
    plotAsTracks = str2double(inputdialog{9});
    matchdist = str2double(inputdialog{10});
    plotSigNoise = str2double(inputdialog{11});
    plotPhotons = str2double(inputdialog{12});
    plotClusters = str2double(inputdialog{13});
    plotOnTimes = str2double(inputdialog{14});


    
    dlg_title = 'Please Input Parameters';
    prompt = {  ...
        'Lobe distance lower bound (in pixel)',...
        'Lobe distance upper bound (in pixel)',...
        'Amplitude ratio limit',...
        'Sigma ratio limit',...
        'Photon weighted fit error lower bound',...
        'Photon weighted fit error upper bound',...
        'Number of photons lower bound',...
        'Number of photons upper bound',...
        'xyPrec lower bound (nm)',...
        'xyPrec upper bound (nm)',...
        'zPrec lower bound (nm)',...
        'zPrec upper bound (nm)',...
        'Minimum template match values',...
        };
    def = { ...
        num2str(lobeDistBounds(1)), ...
        num2str(lobeDistBounds(2)), ...
        num2str(ampRatioLimit), ...
        num2str(sigmaRatioLimit), ...
        num2str(fitErrorRange(1)), ...
        num2str(fitErrorRange(2)), ...
        num2str(numPhotonRange(1)), ...
        num2str(numPhotonRange(2)), ...
        num2str(xyPrecRange(1)),...
        num2str(xyPrecRange(2)),...
        num2str(zPrecRange(1)),...
        num2str(zPrecRange(2)),...
        ['[' num2str(threshVals) ']'],...
        };
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    
    lobeDistBounds = [str2double(inputdialog{1}) str2double(inputdialog{2})];  
    ampRatioLimit = str2double(inputdialog{3});
    sigmaRatioLimit = str2double(inputdialog{4});
    fitErrorRange = [str2double(inputdialog{5}) str2double(inputdialog{6})];
    numPhotonRange = [str2double(inputdialog{7}) str2double(inputdialog{8})];
    xyPrecRange = [str2double(inputdialog{9}) str2double(inputdialog{10})];
    zPrecRange = [str2double(inputdialog{11}) str2double(inputdialog{12})];
    threshVals = str2num(inputdialog{13});
    
    if plotClusters == 1
         dlg_title = 'Cluster Parameters';
         prompt = {'Epsilon (nm)','Minimum # of Points'};
         num_lines =1;
         def = {num2str(400),num2str(8)};
         inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
         epsilon = str2double(inputdialog{1});
         MinPts = str2double(inputdialog{2});
    end
    %% Plot the white light image if specified
    if whiteLightFile ~= 0
        [xWL, yWL] = meshgrid(0:WLHeight,0:WLWidth);
    end
    
    
    %% Now re-evaluate the goodness of the fits
    
    % Conditions for fits (play with these):
    % (1) Amplitude of both lobes > 0
    % (2) All locations x1,y1, x2,y2 lie inside area of small box
    % (3) All sigmas need to be > sigmaBound(1) and < sigmaBound(2)
    % (4) Distance between lobes needs to be > lobeDist(1) pixels and < lobeDist(2) pixels
    % (5) Make sure amplitudes are within 100% of one another
    % (6) Make sure totalFitError/(total number of photons) is within the fitErrorRange
    
    fitErrorCol = 16;
    goodFitFlagCol = 17;
    numPhotonCol = 21;
    bkgndCol = 15;
    lobeDistCol = 22;
    ampRatioCol = 23;
    sigmaRatioCol = 24;
    
    % Empirically determined amplitudes for fitting function based on
    % localization precisision calibration collected on 20120518 on 8a back setup.
    amplitude =  [  361035.867260138,22.2956414971275;...   %   [A1x  A2x]
        348907.934759022,28.3183226442783;...   %   [A1y  A2y]
        840446.405407229,23.3314294806927];      %   [A1z  A2z]
    
    numPhotons = catPSFfits(:,numPhotonCol);
    meanBkgnd = catPSFfits(:,bkgndCol);
    
    if any(meanBkgnd<0)
        meanBkgnd(meanBkgnd<0) = 0;
    end
    
    % Equation 4 of Stallinga and Rieger, ISBI, Barcelona conference proveedings
    sigmaX = sqrt(amplitude(1,1) .* (1./numPhotons) + amplitude(1,1)*4*amplitude(1,2) .* meanBkgnd./(numPhotons).^2 + amplitude(1,1) .* (1./numPhotons) .* sqrt((2*amplitude(1,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(1,2)*(meanBkgnd./numPhotons)))));
    sigmaY = sqrt(amplitude(2,1) .* (1./numPhotons) + amplitude(2,1)*4*amplitude(2,2) .* meanBkgnd./(numPhotons).^2 + amplitude(2,1) .* (1./numPhotons) .* sqrt((2*amplitude(2,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(2,2)*(meanBkgnd./numPhotons)))));
    sigmaZ = sqrt(amplitude(3,1) .* (1./numPhotons) + amplitude(3,1)*4*amplitude(3,2) .* meanBkgnd./(numPhotons).^2 + amplitude(3,1) .* (1./numPhotons) .* sqrt((2*amplitude(3,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(3,2)*(meanBkgnd./numPhotons)))));
    
    %Need to resent goodFitFlagCol so it does not automatically filter out
    %localizations with an error flag from f_fitSMs.  If not reset you will
    %only be able to filter out more localizations, not add more back in.
    %Does not interfere with flags -1001,-1002,or -1003 since they are only
    %calculated in f_fitSMs
    catPSFfits(catPSFfits(:,goodFitFlagCol) < -1003, 17) = 0;
    catPSFfits(catPSFfits(:,goodFitFlagCol) > -1001, 17) = 0;

    
    catPSFfits(catPSFfits(:,sigmaRatioCol) > sigmaRatioLimit,goodFitFlagCol) = -1004;
    catPSFfits(catPSFfits(:,lobeDistCol) < lobeDistBounds(1),goodFitFlagCol) = -1005;
    catPSFfits(catPSFfits(:,lobeDistCol) > lobeDistBounds(2),goodFitFlagCol) = -1005;
    catPSFfits(catPSFfits(:,ampRatioCol) > ampRatioLimit,goodFitFlagCol) = -1006;
    catPSFfits(catPSFfits(:,fitErrorCol)*conversionFactor./numPhotons > fitErrorRange(2),goodFitFlagCol) = -1007;
    catPSFfits(catPSFfits(:,fitErrorCol)*conversionFactor./numPhotons < fitErrorRange(1),goodFitFlagCol) = -1007;
    catPSFfits(sigmaX < xyPrecRange(1) | sigmaX > xyPrecRange(2) |...
        sigmaY < xyPrecRange(1) | sigmaY > xyPrecRange(2) |...
        sigmaZ < zPrecRange(1) | sigmaZ > zPrecRange(2),goodFitFlagCol) = -1008;
    catPSFfits(catPSFfits(:,templateStrCol) < threshVals(catPSFfits(:,templateNumCol))',goodFitFlagCol) = -1009;
    catPSFfits(catPSFfits(:,goodFitFlagCol)>= 0,goodFitFlagCol) = 3;


    %% load valid xyz locations
        
    goodFits = catPSFfits(:,17) >= 0; 
    badFits = catPSFfits(:,17) < 0;
    goodFits = goodFits & catPSFfits(:,1) >= frameRange(1) & catPSFfits(:,1) <= frameRange(2);
    goodFitsNoPhotFilt = sum(goodFits);
    goodFits = goodFits & catPSFfits(:,numPhotonCol) >= numPhotonRange(1) & catPSFfits(:,numPhotonCol) <= numPhotonRange(2);
    if goodFitsNoPhotFilt - sum(goodFits) > goodFitsNoPhotFilt / 2
        warning(['More than half of the fits were thrown out due to '...
            'restrictions on the # photons. Double check this limit.']);
    end
    if sum(goodFits) < 5
        warning('Very few (<5) fits passed the filters. Double-check limits.');
    end
    % corrects zRange for index mismatch (see below for the inverse
    % transformation to the z position)
    corrzRange = zRange * nOil/nSample;
    
    % extract data; check if fiducial correction is in play (i.e., col 30 is not all nans)
    % notate fid-corrected localizations separately
    if sum(isnan(catPSFfits(:,30))) ~= size(catPSFfits,1)
        goodFits = goodFits & catPSFfits(:,30) >= corrzRange(1) & catPSFfits(:,30) <= corrzRange(2); % try to fid-corrected values for z range
        xLocPix = catPSFfits(goodFits,18)/nmPerPixel;
        yLocPix = catPSFfits(goodFits,19)/nmPerPixel;
        goodPSFfits = catPSFfits(goodFits,:);
        xLoc = catPSFfits(goodFits,28);
        yLoc = catPSFfits(goodFits,29);
        zLoc = catPSFfits(goodFits,30);
        xLoc_bad = catPSFfits(badFits,28);
        yLoc_bad = catPSFfits(badFits,29);
        % still generate xLoc etc. as below, but call them 'raw' if fids
        % were used. these can then be used for registration, where both
        % focal shift correction and fiducial correction should be
        % performed AFTER registering the two channels.
        xLocPixRaw = catPSFfits(goodFits,18)/nmPerPixel;
        yLocPixRaw = catPSFfits(goodFits,19)/nmPerPixel;
        xLocRaw = catPSFfits(goodFits,25);
        yLocRaw = catPSFfits(goodFits,26);
        zLocRaw = catPSFfits(goodFits,27);
        goodPSFfits = catPSFfits(goodFits,:);
        xLoc_badRaw = catPSFfits(badFits,25);
        yLoc_badRaw = catPSFfits(badFits,26);
        
    else % if no fid-corrected traces were generated (column 30 is nan)
        goodFits = goodFits & catPSFfits(:,27) >= corrzRange(1) & catPSFfits(:,27) <= corrzRange(2);
        xLocPix = catPSFfits(goodFits,18)/nmPerPixel;
        yLocPix = catPSFfits(goodFits,19)/nmPerPixel;
        xLoc = catPSFfits(goodFits,25);
        yLoc = catPSFfits(goodFits,26);
        zLoc = catPSFfits(goodFits,27);
        goodPSFfits = catPSFfits(goodFits,:);
        xLoc_bad = catPSFfits(badFits,25);
        yLoc_bad = catPSFfits(badFits,26);
    end
    
    clear corrzRange
    zLoc_IndexCorrected = zLoc * nSample/nOil;
    if exist('xLocRaw','var')
        zLoc_IndexCorrectedRaw = zLocRaw * nSample/nOil;
    end
    
    numPhotons = catPSFfits(goodFits,21);
    lobeDist = catPSFfits(goodFits,22);
    ampRatio = catPSFfits(goodFits,23);
    sigmaRatio = catPSFfits(goodFits,24);
    
    meanBkgnd = catPSFfits(goodFits,15);  % output from template match is already in units of photons.
    frameNum = catPSFfits(goodFits,1);
    PSFfits_bad = catPSFfits(badFits,:);
    
    %% ask user what region to plot in superresolution image
    h2Dfig=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    if whiteLightFile ~= 0
        xRange = xWL(1,:)*nmPerPixel-wlShiftX;
        yRange = yWL(:,1)*nmPerPixel-wlShiftY;
        % pick region that contains background
        imagesc(xRange,yRange,whiteLight,[0 1]);colormap gray;hold on;
    end
    %     plot(xLoc,yLoc,'.','MarkerSize',1);
    scatter((xLoc(:)), (yLoc(:)), '.');
    xlim([min(xLoc(:))-500 max(xLoc(:))+500]);
    ylim([min(yLoc(:))-500 max(yLoc(:))+500]);
    xlabel('x (nm)');ylabel('y (nm)');
    axis ij;
    axis square;
    
    
    if pass == 1
        ROI = imrect(gca,[min(xLoc(:)) min(yLoc(:)) max(xLoc(:))-min(xLoc(:)) max(yLoc(:))-min(yLoc(:))]);
    else
        ROI = imrect(gca,[ROI(1) ROI(2) ROI(3) ROI(4)]);
    end
    
    title({'Double-click to choose region that will be plotted in 3D scatterplot' ...
        mat2str(ROI.getPosition)});
    addNewPositionCallback(ROI,@(p) title({'Double-click to choose region that will be plotted in 3D scatterplot' ...
        ['[xmin ymin width height] = ' mat2str(p,3)]}));
    % make sure rectangle stays within image bounds
    fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(ROI,fcn);
    ROI = wait(ROI);
    clear avgImg fcn
    
    %% filter out localizations outside of ROI
    
    validPoints = xLoc>=ROI(1) & xLoc<=ROI(1)+ROI(3) & yLoc>ROI(2) & yLoc<=ROI(2)+ROI(4) & numPhotons>0;
    invalidPoints = xLoc_bad>=ROI(1) & xLoc_bad<=ROI(1)+ROI(3) & yLoc_bad>ROI(2) & yLoc_bad<=ROI(2)+ROI(4) ;
    
    if ~any(validPoints)
        disp('You chose an area without any points');
        continue
    end
    
    xLocPix = xLocPix(validPoints);
    yLocPix = yLocPix(validPoints);
    xLoc = xLoc(validPoints);
    yLoc = yLoc(validPoints);
    zLoc = zLoc(validPoints);
    goodPSFfits = goodPSFfits(validPoints,:);
    zLoc_IndexCorrected = zLoc_IndexCorrected(validPoints);
    
    if exist('xLocRaw','var')
        xLocPixRaw = xLocPixRaw(validPoints);
        yLocPixRaw = yLocPixRaw(validPoints);
        xLocRaw = xLocRaw(validPoints);
        yLocRaw = yLocRaw(validPoints);
        zLocRaw = zLocRaw(validPoints);
        zLoc_IndexCorrectedRaw = zLoc_IndexCorrectedRaw(validPoints);
    end
    
    numPhotons = numPhotons(validPoints);
    meanBkgnd = meanBkgnd(validPoints);
    frameNum = frameNum(validPoints);
    lobeDist = lobeDist(validPoints);
    ampRatio = ampRatio(validPoints);
    sigmaRatio = sigmaRatio(validPoints);
    
    PSFfits_bad = PSFfits_bad(invalidPoints,:);
    hRejections = figure;
    subplot(8,1,1)
    x = -1008:1:-1001;
    hist(PSFfits_bad(PSFfits_bad(:,17)<-10,17),x)
    xlabel('Error Flag');ylabel('Frequency');
    title({[num2str(size(PSFfits_bad,1)) ' bad localizations']});
    hold on;
    for t = 1:numTemplates
        subplot(8,1,t+1);
        x = -1008:1:-1001;
        hist(PSFfits_bad(PSFfits_bad(:,17)<-10 & PSFfits_bad(:,5) == t,17),x)
        xlabel('Error Flag');ylabel('Frequency');
        title({['Template ' num2str(t) ': ' num2str(size(PSFfits_bad(PSFfits_bad(:,5) == t),1)) ' bad localizations (' num2str(length(goodPSFfits(goodPSFfits(:,5) == t))) ' total)']});
        hold on;
    end

    
    hTemplateFit = figure;
    subplot(7,1,1)
    hold on;
    for t = 1:numTemplates
        subplot(7,1,t);
        xbin = min(threshVals):(min(threshVals)*5)/500:(min(threshVals)*5);
        histogram(goodPSFfits((goodPSFfits(:,5) == t),6),xbin)
        xlabel('Fit Value');ylabel('Frequency');
        title({['Template ' num2str(t) ' Fit Values']});
        hold on;
    end
    
    
    %% compute localization precision as a function of the number of photons
    
    amplitude =  [  361035.867260138,22.2956414971275;...   %   [A1x  A2x]
        348907.934759022,28.3183226442783;...   %   [A1y  A2y]
        840446.405407229,23.3314294806927];      %   [A1z  A2z]
    
    % Equation 4 of Stallinga and Rieger, ISBI, Barcelona conference proveedings
    if any(meanBkgnd<0);
        correctedBG = find(meanBkgnd<0);
        warning([num2str(length(correctedBG))...
            'values of meanBkgnd (indices below) were negative! Changing to 0.']);
        meanBkgnd(correctedBG)=0;
    end
    
    sigmaX = sqrt(amplitude(1,1) .* (1./numPhotons) + amplitude(1,1)*4*amplitude(1,2) .* meanBkgnd./(numPhotons).^2 + amplitude(1,1) .* (1./numPhotons) .* sqrt((2*amplitude(1,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(1,2)*(meanBkgnd./numPhotons)))));
    sigmaY = sqrt(amplitude(2,1) .* (1./numPhotons) + amplitude(2,1)*4*amplitude(2,2) .* meanBkgnd./(numPhotons).^2 + amplitude(2,1) .* (1./numPhotons) .* sqrt((2*amplitude(2,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(2,2)*(meanBkgnd./numPhotons)))));
    sigmaZ = sqrt(amplitude(3,1) .* (1./numPhotons) + amplitude(3,1)*4*amplitude(3,2) .* meanBkgnd./(numPhotons).^2 + amplitude(3,1) .* (1./numPhotons) .* sqrt((2*amplitude(3,2)*(meanBkgnd./numPhotons))./(1+(4*amplitude(3,2)*(meanBkgnd./numPhotons)))));
    
    
    
    meanNumPhotons = mean(numPhotons);
   
    
    %% plot statistics on these localizations
    
    hStatsFig=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
    subplot(1,4,1);
    hist(numPhotons,round(length(xLoc)/20));
    xlabel('Number of photons per localization');
    subplot(1,4,2);
    hist(sigmaX,round(length(xLoc)/20));
    xlabel('\sigma_x (nm)');
    subplot(1,4,3);
    hist(sigmaY,round(length(xLoc)/20));
    xlabel('\sigma_y (nm)');
    subplot(1,4,4);
    hist(sigmaZ,round(length(xLoc)/20));
    xlabel('\sigma_z (nm)');
        
    
    if plotSigNoise == 1
        sigNoise =numPhotons./meanBkgnd;
        
        %plot signal to noise for each good localization, rounded to nearest
        %integer. Limits range due to some rare, extremely large values.
        sigNoiseRound = round(sigNoise);
        xLocSigNoise = xLoc(sigNoiseRound < 3000 & sigNoiseRound > -3000,1);
        yLocSigNoise = yLoc(sigNoiseRound < 3000 & sigNoiseRound > -3000,1);
        sigNoiseRound = sigNoiseRound(sigNoiseRound < 3000 & sigNoiseRound > - 3000);
        
        %histogram for signal to noise
        figure;
        histogram(sigNoiseRound);
        xlabel('Signal/Noise');
        ylabel('Counts');
        title('Signal to Noise for Each Localization');
        
        %plot of S/N over FOV 
        sigNoiseRoundColors = zeros(length(sigNoiseRound),3);
        figure;
        cLength(:,1) = unique(sigNoiseRound);
        cLength(:,2) = 1:length(cLength);
        cColors = jet(length(cLength));
        for k = 1:length(sigNoiseRound)
            sigNoiseRoundColors(k,:) = cColors(cLength(sigNoiseRound(k,1) == cLength(:,1),2),:);
        end
        c = sigNoiseRoundColors;
        a = 36;
        scatter(xLocSigNoise,yLocSigNoise,a,c);
        xlabel('x');
        ylabel('y');
        title('Signal/Noise');
        minLabel = min(cLength(:,1));
        maxLabel = max(cLength(:,1));
        tickLabels = linspace(minLabel,maxLabel,5);
        colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',tickLabels);
        colormap jet; 
        
        clear cLength;
    end
    
    if plotPhotons == 1 
        %plot numPhotons over FOV 
        numPhotonsColors = zeros(length(numPhotons),3);
        figure;
        cLength(:,1) = unique(numPhotons);
        cLength(:,2) = 1:length(cLength);
        cColors = jet(length(cLength));
        for k = 1:length(numPhotons)
            numPhotonsColors(k,:) = cColors(cLength(numPhotons(k,1) == cLength(:,1),2),:);
        end
        c = numPhotonsColors;
        a = 36;
        scatter(xLoc,yLoc,a,c);
        xlabel('x');
        ylabel('y');
        title('Number of Photons');
        minLabel = min(cLength(:,1));
        maxLabel = max(cLength(:,1));
        tickLabels = linspace(minLabel,maxLabel,5);
        colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',tickLabels);
        colormap jet;
        
        clear cLength;
    end
    
     figure;
     histogram(goodPSFfits(:,1),length(frames));
     xlabel('Frame #');
     ylabel('Frequency');
     title('Number of Localizations Per Frame');
     
     if plotOnTimes == 1
         tracking = cell(length(PSFvals),1);
         parfor a = 1:length(PSFvals)-3
             if ~isempty(PSFvals{a,1}{1,2}) && ~isempty(PSFvals{a+1,1}{1,2})
                 for b = 1:length(PSFvals{a,1}{1,2})
                     if any(goodPSFfits(:,2)==b) && any(goodPSFfits(:,1)==a)
                         trackMol = find(goodPSFfits(:,1)==a & goodPSFfits(:,2)==b);
                         trackMolx = xLoc(trackMol);
                         trackMoly = yLoc(trackMol);
                         trackMolz = zLoc_IndexCorrected(trackMol);
                         checkMolmatch = 0;
                         for c = 1:length(PSFvals{a+1}{1,2})
                             checkMol = find(goodPSFfits(:,1)==a+1 & goodPSFfits(:,2)==c);
                             checkMolx = xLoc(checkMol);
                             checkMoly = yLoc(checkMol);
                             checkMolz = zLoc_IndexCorrected(checkMol);
                             distMol = sqrt((trackMolx-checkMolx).^2+(trackMoly-checkMoly).^2+(trackMolz-checkMolz).^2);
                             if distMol < matchdist
                                 trackingVals = [a b trackMolx trackMoly trackMolz a+1 c checkMolx checkMoly checkMolz];
                                 tracking{a,1}{b,c}(1,:) = trackingVals;
                                 checkMolmatch = checkMolmatch+1;
                             end
                         end
%                         
                     end
                 end
             end
         end
           
         trackingFits = [];
         for a = 1:length(tracking)-1
             if ~isempty(tracking{a,1})
                 for b = 1:size(tracking{a,1},1)
                     for c = 1:size(tracking{a,1},2)
                         if ~isempty(tracking{a,1}{b,c})
                             trackingFits = cat(1,trackingFits,tracking{a,1}{b,c}(1,:));
                         end
                     end
                 end
             end
         end
         
         %start actual tracking
         numTrack = 0;
         trackingFits(:,11) = zeros(length(trackingFits(:,1)),1);
         for a=1:length(trackingFits(:,1))-1
             if trackingFits(a,11) == 0
                 numTrack = numTrack + 1;
                 trackingFits(a,11) = numTrack;
             end
             for b=1:length(trackingFits(:,1))-a
                 if trackingFits(a+b,1:2) == trackingFits(a,6:7)
                     trackingFits(a+b,11) = trackingFits(a,11);
                     
                 end
             end
         end
         
        markerColors = jet(length(unique(trackingFits(:,11)))-1);
        figure; hold on;
        for a = 1:length(unique(trackingFits(:,11)))-1
            plot3([trackingFits(trackingFits(:,11)==a,3), trackingFits(trackingFits(:,11)==a,8)],[trackingFits(trackingFits(:,11)==a,4), trackingFits(trackingFits(:,11)==a,9)],[trackingFits(trackingFits(:,11)==a,5), trackingFits(trackingFits(:,11)==a,10)],'LineWidth',scatterSize/12,...
                'Color',markerColors(a,:));
        end
        axis vis3d equal;
         %end tracking
         
        
         trackingFitsMean(:,1) = (trackingFits(:,3) + trackingFits(:,8))/2;
         trackingFitsMean(:,2) = (trackingFits(:,4) + trackingFits(:,9))/2;
         trackingFitsMean(:,3) = (trackingFits(:,5) + trackingFits(:,10))/2;
         trackingFitsMean(:,4) = trackingFits(:,1);
         molNum = 1;
         trackingFitsMean(:,5) = 0;
         for aa = 1:length(trackingFitsMean)-3
             if trackingFitsMean(aa,5) ==0;
                 trackingFitsMean(aa,5) = molNum;
                 molNum = molNum +1;
             end
                 bb = 1;
                 while (trackingFitsMean(aa+bb,4)-trackingFitsMean(aa,4) <= 1) && aa+bb <= length(trackingFitsMean)
                     distance = sqrt((trackingFitsMean(aa,1)-trackingFitsMean(aa+bb,1)).^2+(trackingFitsMean(aa,2)-trackingFitsMean(aa+bb,2)).^2+(trackingFitsMean(aa,3)-trackingFitsMean(aa+bb,3)).^2);
                    if distance < matchdist && (trackingFitsMean(aa+bb,4)-trackingFitsMean(aa,4) == 1)
                     trackingFitsMean(aa+bb,5) = trackingFitsMean(aa,5);
                    end
                     bb = bb +1;
                 end
         end
         
         totalNum = length(xLoc);
         totalCount = zeros(1,max(trackingFitsMean(:,5)));
         
         uniqueTrackingFitsMean = zeros(max(trackingFitsMean(:,5)),8);
         for dd = 1:max(trackingFitsMean(:,5))
             totalCount(1,dd) = sum(trackingFitsMean(:,5) == dd) + 1;
             uniqueTrackingFitsMean(dd,1:5) = mean(trackingFitsMean(trackingFitsMean(:,5) == dd,:),1);
             uniqueTrackingFitsMean(dd,6) = sum(trackingFitsMean(:,5) == dd)+1;
         end
         
         nmPerPixel = 108;
         uniqueTrackingFitsMean(:,7:8) = uniqueTrackingFitsMean(:,1:2);
         uniqueTrackingFitsMean(:,1:2) = round(uniqueTrackingFitsMean(:,1:2)/nmPerPixel);
         uniqueTrackingFitsImage = zeros(imgWidth,imgHeight);
         for track = 1:length(uniqueTrackingFitsMean(:,1))
             if uniqueTrackingFitsImage(uniqueTrackingFitsMean(track,1),uniqueTrackingFitsMean(track,2)) == 0;
                 uniqueTrackingFitsImage(uniqueTrackingFitsMean(track,1),uniqueTrackingFitsMean(track,2)) = 1;
             else
                 uniqueTrackingFitsImage(uniqueTrackingFitsMean(track,1),uniqueTrackingFitsMean(track,2)) = uniqueTrackingFitsImage(uniqueTrackingFitsMean(track,1),uniqueTrackingFitsMean(track,2)) + 1;
             end
         end
         
         
         maxNum = max(totalCount(1,:));
         onHist = [];
         onHist(1,1) = totalNum - length(trackingFitsMean);
         onHist2(1:onHist(1,1),1) = 1;
         for ee = 2:maxNum
             onHist(ee,1) = sum(totalCount == ee);
             onHistLength = length(onHist2(:,1));
             onHist2(onHistLength+1:onHistLength+onHist(ee,1),1) = ee;
         end
         
         onHist2 = onHist2.*exposureTime;
         
         figure;
         edges = onHist2(1)/2:onHist2(1):onHist2(end)+onHist2(1)/2;
         histogram(onHist2,edges);
         title(['On Times: ' num2str(exposureTime) ' ms frames']);
         xlabel('Time On (ms)');
         ylabel('# of Occurances');
         
         
         clear trackingFitsMean trackingFits onHist onHist2;
     end
    
    %% plot 3D scatterplot of localizations with white light
    
    h3Dfig = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'renderer','opengl', 'Toolbar', 'figure');
    if whiteLightFile~=0
        [x,y,z] = meshgrid(xRange,yRange,[min(zLoc_IndexCorrected)-100 max(zLoc_IndexCorrected)]);
        %MAKE SURE x, y, z, and b ARE ALL THE EXACT SAME SIZE
        xslice = []; yslice = []; zslice = [min(zLoc_IndexCorrected)-100];
        b = repmat(whiteLight,[1 1 2]);
        b(length(b(:,1,:))+1,:,:) = 0;
        b(:,length(b(1,:,:))+1,:) = 0;

        h=slice(x,y,z,b,xslice,yslice,zslice);
        set(h,'EdgeColor','none','FaceAlpha',0.75);
        colormap(h3Dfig,gray); hold on;
    end
    
    if plotAsTracks == 1

        tracking = cell(length(PSFvals),1);
        parfor a = 1:length(PSFvals)-1
            if ~isempty(PSFvals{a,1}{1,2}) && ~isempty(PSFvals{a+1,1}{1,2})
                for b = 1:length(PSFvals{a,1}{1,2})
                    if any(goodPSFfits(:,2)==b) && any(goodPSFfits(:,1)==a)
                        trackMol = find(goodPSFfits(:,1)==a & goodPSFfits(:,2)==b);
                        trackMolx = xLoc(trackMol);
                        trackMoly = yLoc(trackMol);
                        trackMolz = zLoc_IndexCorrected(trackMol);
                        checkMolmatch = 0;
                        for c = 1:length(PSFvals{a+1}{1,2})
                            checkMol = find(goodPSFfits(:,1)==a+1 & goodPSFfits(:,2)==c);
                            checkMolx = xLoc(checkMol);
                            checkMoly = yLoc(checkMol);
                            checkMolz = zLoc_IndexCorrected(checkMol);
                            distMol = sqrt((trackMolx-checkMolx).^2+(trackMoly-checkMoly).^2+(trackMolz-checkMolz).^2);
                            if distMol < matchdist
                                trackingVals = [a b trackMolx trackMoly trackMolz a+1 c checkMolx checkMoly checkMolz];
                                tracking{a,1}{b,c}(1,:) = trackingVals;
                                checkMolmatch = checkMolmatch+1;
                            end
                        end
                        if checkMolmatch == 0 && ~isempty(PSFvals{a+2,1}{1,2}) && length(PSFvals) >= (a+2)
                            for c = 1:length(PSFvals{a+2}{1,2})
                                checkMol = find(goodPSFfits(:,1)==a+2 & goodPSFfits(:,2)==c);
                                checkMolx = xLoc(checkMol);
                                checkMoly = yLoc(checkMol);
                                checkMolz = zLoc_IndexCorrected(checkMol);
                                distMol = sqrt((trackMolx-checkMolx).^2+(trackMoly-checkMoly).^2+(trackMolz-checkMolz).^2);
                                if distMol < matchdist
                                    trackingVals = [a b trackMolx trackMoly trackMolz a+2 c checkMolx checkMoly checkMolz];
                                    tracking{a,1}{b,c}(1,:) = trackingVals;
                                    checkMolmatch = checkMolmatch+1;
                                end
                            end
                        end
                        if checkMolmatch == 0 && ~isempty(PSFvals{a+3,1}{1,2}) && length(PSFvals) >= (a+3)
                            for c = 1:length(PSFvals{a+2}{1,2})
                                checkMol = find(goodPSFfits(:,1)==a+3 & goodPSFfits(:,2)==c);
                                checkMolx = xLoc(checkMol);
                                checkMoly = yLoc(checkMol);
                                checkMolz = zLoc_IndexCorrected(checkMol);
                                distMol = sqrt((trackMolx-checkMolx).^2+(trackMoly-checkMoly).^2+(trackMolz-checkMolz).^2);
                                if distMol < matchdist
                                    trackingVals = [a b trackMolx trackMoly trackMolz a+3 c checkMolx checkMoly checkMolz];
                                    tracking{a,1}{b,c}(1,:) = trackingVals;
                                    checkMolmatch = checkMolmatch+1;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        
        trackingFits = [];
        for a = 1:length(tracking)-1
            if ~isempty(tracking{a,1})
                for b = 1:size(tracking{a,1},1)
                    for c = 1:size(tracking{a,1},2)
                        if ~isempty(tracking{a,1}{b,c})
                            trackingFits = cat(1,trackingFits,tracking{a,1}{b,c}(1,:));
                        end
                    end
                end
            end
        end    
        markerColors = jet(length(PSFvals));
        for a = 1:length(trackingFits)
            plot3([trackingFits(a,3), trackingFits(a,8)],[trackingFits(a,4), trackingFits(a,9)],[trackingFits(a,5), trackingFits(a,10)],'LineWidth',scatterSize/12,...
                'Color',markerColors(trackingFits(a,1),:));
            set(gca,'FontSize',15)
            hold on;
        end
        
     elseif plotClusters == 1
        %uses DBSCAN to find clusters in 3D.
         X = [xLoc,yLoc,zLoc_IndexCorrected];    
         [IDX, isnoise]=DBSCAN(X,epsilon,MinPts);
         
         %plot clusters
%          clusterColor = jet(max(IDX));
%          figure; 
         hold on;
         for c = 1:max(IDX)
             tTrack = frameNum(IDX(:) ==c);
             xLocClust = xLoc(IDX(:) == c);
             yLocClust = yLoc(IDX(:) == c);
             zLocClust = zLoc_IndexCorrected(IDX(:) == c);
             plot3((xLocClust(:)),(yLocClust(:)),zLocClust(:),'.','MarkerSize',scatterSize/3,...
            'Color',[1 0 0]); 

         end
             axis vis3d equal;
           
    elseif useTimeColors == 0
        plot3((xLoc(:)),(yLoc(:)),zLoc_IndexCorrected(:),'.','MarkerSize',scatterSize/3,...
            'Color',[1 0 0]);
        %     scatter3(xLoc,yLoc,zLoc,scatterSize,[1 1 0],'filled');
    elseif useTimeColors == 1
        %         scatter3(xLoc(a),yLoc(a),zLoc(a),scatterSize,frameNum(1):frameNum(length(frameNum)),'filled')
        markerColors = jet(frameNum(length(frameNum))-frameNum(1)+1);
        %     for a = 1:length(frameNum)
        %         scatter3(xLoc(a)-min(xLoc(:)),yLoc(a)-min(yLoc(:)),zLoc(a),scatterSize,'filled',...
        %             'MarkerFaceColor', markerColors(frameNum(a)-frameNum(1)+1,:),...
        %             'MarkerEdgeColor', markerColors(frameNum(a)-frameNum(1)+1,:));
        %     end
        hold on
        for a = 1:length(frameNum)
            scatter3(xLoc(a),yLoc(a),zLoc_IndexCorrected(a),scatterSize,'filled',...
                'MarkerFaceColor', markerColors(frameNum(a)-frameNum(1)+1,:),...
                'MarkerEdgeColor', markerColors(frameNum(a)-frameNum(1)+1,:));
            ticks = 0:zRange(2)/10:zRange(2);
            c = colorbar('southoutside');
%             colormap jet;
        end
    end
    axis vis3d equal;
    xlim([ROI(1) ROI(1)+ROI(3)]);
    ylim([ROI(2) ROI(2)+ROI(4)]);
    xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
    
    ROICenterX = ROI(1)+ROI(3)/2;
    ROICenterY = ROI(2)+ROI(4)/2;
    

    if exist('laserProfile', 'var')
        
        meanIntensityInROI = sum(sum(laserProfile(iROI(1,1):(iROI(1,1)+length(FOVmask(:,1))-1),iROI(1,2):(iROI(1,2)+length(FOVmask(1,:))-1))))/(length(FOVmask(:,1))*length(FOVmask(1,:)));
        title({[num2str(length(xLoc)) ' localizations'];...
            ['Mean Number of Signal Photons = ' num2str(meanNumPhotons) ' per frame'];...
            ['Mean Number of Background Photons = ' num2str(mean(meanBkgnd)) ' per pixel per frame'];...
            ['Localization Precision \sigma_x = ' num2str(mean(sigmaX)) ' nm'];...
            ['Localization Precision \sigma_y = ' num2str(mean(sigmaY)) ' nm'];...
            ['Localization Precision \sigma_z = ' num2str(mean(sigmaZ)) ' nm'];...
            ['Laser Intensity = ' num2str(meanIntensityInROI) ' W/m^2']},...
            'color','k');
    else
        
        title({[num2str(length(xLoc)) ' localizations'];...
            ['Mean Number of Signal Photons = ' num2str(meanNumPhotons) ' per frame'];...
            ['Mean Number of Background Photons = ' num2str(mean(meanBkgnd)) ' per pixel per frame'];...
            ['Localization Precision \sigma_x = ' num2str(mean(sigmaX)) ' nm'];...
            ['Localization Precision \sigma_y = ' num2str(mean(sigmaY)) ' nm'];...
            ['Localization Precision \sigma_z = ' num2str(mean(sigmaZ)) ' nm']},...
            'color','k');
    end
    
    set(gca,'color','white');
    set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
    

    
    %% Construct a questdlg with three options
    
    h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);

    
    dlg_title = 'Replot';
    prompt = {'Would you like to replot with a different parameter set?'};
    def =       { 'Yes'  };
    questiondialog = questdlg(prompt,dlg_title, def);
    % Handle response
    switch questiondialog
        case 'Yes'
            pass = pass + 1;
            close
            close
            close
            close
            close
            close
            close
            close
        case 'No'
            anotherpass = false;
        case 'Cancel'
            error('User cancelled the program');
    end
    

end
%% prompt to save data
[saveFile, savePath] = uiputfile({'*.*'},'Enter a directory title for this ROI. Otherwise, click cancel.');
savePath = [savePath saveFile '/'];
mkdir(savePath);
if ~isequal(saveFile,0)
    save([savePath 'Output'],'xLocPix','yLocPix','xLoc','yLoc','zLoc','zLoc_IndexCorrected','numPhotons','meanBkgnd','sigmaX','sigmaY','sigmaZ','frameNum',...
        'zRange','frameRange','sigmaBounds','lobeDistBounds','ampRatioLimit','sigmaRatioLimit','fitErrorRange','numPhotonRange',...
        'lobeDist','ampRatio','sigmaRatio','wlShiftX', 'wlShiftY','goodFits','fidTrackX', 'fidTrackY', 'fidTrackZ', 'nmPerPixel','whiteLightFile','threshVals',...
        'spatialCorr','useCurrent','currFidIdx');
    if exist('xLocRaw');
        save([savePath 'Output'],'xLocRaw','yLocRaw','zLocRaw','zLoc_IndexCorrectedRaw','-append');
    end
    if exist('whiteLight');
        save([savePath 'Output'],'whiteLight','xWL','yWL','-append');
    end
end
%%

% Save figures
close(hTemplateFit);

saveas(h3Dfig,[savePath '3D.fig']);
close(h3Dfig)
saveas(hStatsFig,[savePath 'stats.fig']);
close(hStatsFig)
saveas(hRejections,[savePath 'rejections.fig']);
close(hRejections)
saveas(h2Dfig,[savePath '2D.fig'])
close(h2Dfig)


end
