tic;

isNoise = 1;

%create FOV mask
load('FOVmask_SIM.mat');

% usePolyROI = 1;
% [testFile, testPath] = uigetfile({'*.mat';'*.*'},'Open test file to create FOVmask');
% load ([testPath,testFile]);
% 
% if usePolyROI
% %     hFOVmaskFig=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
%     hFOVmaskFig=figure;
% %     imagesc(totalImg{1},[0 max(0, min(avgImg(:)))+5*std(avgImg(:))]);
%     totalImg{1} = totalImg{1}(1:96,1:96);
%     imagesc(totalImg{1});
%     axis image;colorbar
%     colormap hot;
%     title('Select ROIpoly of area to keep');
%     [FOVmask, maskX, maskY] = roipoly;
%     xCenter=(max(maskX)+min(maskX))/2;
%     yCenter=(max(maskY)+min(maskY))/2;
%     x1=(maskX-xCenter)*0.95+xCenter;
%     y1=(maskY-yCenter)*0.95+yCenter;
%     FOVmask1=roipoly(totalImg{1},x1,y1);
%     close(hFOVmaskFig);
% end

[zFile, zPath] = uigetfile({'*.mat';'*.*'},'Open zcalfile');
load([zPath zFile]);
s.calFilePrefix = zPath;

[beadFile, beadPath] = uigetfile({'*.mat';'*.*'},'Open bead 1 templates file');
load([beadPath beadFile]);

[dataFiles, dataPathOrig] = uigetfile({'*.mat';'*.*'},'MultiSelect', 'on','Open data files');

if iscell(dataFiles)
    loopNum = length(dataFiles);
else
    loopNum = 1;
end
for j = 1:loopNum

% file where data is saved
projFile = '';
% all saveable data is originally organized in a structure 's'
% if a channel is selected, then the 's' data is duplicated to that
% channel.
% list of TIFs with raw SM data
s.smacmRawFile = {};
s.smacmRawPath = {};
% list of corresponding files with dark counts
s.smacmDarkFile = {};
s.smacmDarkPath = {};
% (optional) sequence file with shutter state, z position for each frame
s.smacmSifFile = {};
s.smacmSifPath = {};
% processed SMACM localizations (corresponding to each raw TIF)
s.smacmLocFile = {};
% concatenated SMACM localizations (fid corrected if available)
s.smacmFullLocFile = '';
% DHPSF calibration data, also points to where templates are saved
% s.calFilePrefix = '';
% DHPSF fiducial tracking data location
s.fidFilePrefix = {};
% DHPSF SM fit data location
s.fitFilePrefix = {};
% DHPSF threshold data save location
s.threshFilePrefix = {};


% status of fitting project: 
% cal, fid, thresh, template match, saved output, saved project
s.projStatus = false(1,5);
% calibration bead images
s.templateImgs = [];
% number of calibration beads to choose from
s.numCalBeads = 1;
% selected calibration bead
s.calBeadIdx = 1;
% use fiducials?
s.useFids = false;
% list of template index numbers, corresponds to first dimension of
% templateImgs[]
s.templateIdxs = [];
% list of template thresholds (columns) for each raw smacm file (rows)
s.templateThreshs = [];
% selected file for specifying template matching thresholds
s.threshFileSelIdx = 1;
% selected template
s.templateSelIdx = 1;
% locations of peaks within templates
s.templateLocs = zeros(6,5);
% region of interest inside raw TIF to process
s.smacmRawROI = [0 0 0 0];
% EM gain setting used when acquiring SMACM data
s.smacmEMGain = 1;
% photons/count, camera setting, global to all modules
s.conversionGain = 0.5; %8A % 24.7; % 8B
% s.conversionGain = 26.93; %8A % 24.7; % 8B
% imaging system property, global to all modules
s.nmPerPixel = 108;  % old value is 125.78; 8B back = 160
% s.nmPerPixel = 12;  % old value is 125.78; 8B back = 160

% channel identifier
s.channel = '0';
%use SSIM if you want to automatically fill in the threshold values after
%f_calSMIdentification
s.tempThresh = [];
%Universal switch between MLE fitting method and LSQ non linear fitting
%method
%MLE fitting method uses a double gaussian point spread function and uses
%units of Lambda (see nature methods Bewersdorf paper for details)
%LSQ nonlinear fitting method uses a double gaussian point spread function
%and uses units of Photons as its counts
s.fittingMethod = 'LSQ with DG model';
% [minWidth maxWidth] of the two spots of the DHPSF, units of pixels
% relative to the original value of 160 nm / pix
% ***This has been set to a constant value based upon Moerner lab DH
% microscopes; different implementations may vary***
s.sigmaBounds = [1.0 1.5];
% [minSpacing maxSpacing] between the two spots of the DHPSF, units of
% pixels relative to the original value of 160 nm / pix
s.lobeDistBounds = [3.5 10]*160/s.nmPerPixel;   %[3.5 10]*160/s.nmPerPixel;
% half-width of box to extract when fitting DHPSF images, units of integer pixels
% it's independent of the magnification of the optical setup, so just
% hard-code a radius for now
s.boxRadius = 13;        % round(7*160/s.nmPerPixel);
% smoothing filter width for identifying DHPSF SMs, units of pixels
s.gaussianFilterSigma = 1.5*160/s.nmPerPixel;
% minimum lateral distance between identified SMs, units of pixels
s.minDistBetweenSMs = 7.5*160/s.nmPerPixel;
%
% channelChoices = ['0';'r';'g'];

% s.nmPerPixel = 108;
% update all dependent parameters
s.sigmaBounds = [1.0 1.5]*160/s.nmPerPixel;
s.lobeDistBounds = [3.5 10]*160/s.nmPerPixel;
s.boxRadius = 13;        % round(7*160/s.nmPerPixel);
s.gaussianFilterSigma = 1.5*160/s.nmPerPixel;
s.minDistBetweenSMs = 7.5*160/s.nmPerPixel;
% s.projStatus(5) = false;


%% SMID outputs

% [zFile, zPath] = uigetfile({'*.mat';'*.*'},'Open zcalfile');
% load([zPath zFile]);
% s.calFilePrefix = zPath;
% 
% [beadFile, beadPath] = uigetfile({'*.mat';'*.*'},'Open bead 1 templates file');
% load([beadPath beadFile]);
% 
% [dataFile, dataPath] = uigetfile({'*.mat';'*.*'},'Open data file');
if iscell(dataFiles)
dataFile = dataFiles{j};
else
dataFile = dataFiles;
end
% load ([dataPathOrig,dataFile]);
s.smacmRawFile = cell(1);
s.smacmRawFile{1} =[dataPathOrig dataFile];
s.smacmRawPath = dataPathOrig;
s.smacmRawROI = [1,1,96,96];


%compute template frames
      goodFit_forward = logical(goodFit_f(s.calBeadIdx,:));
            templateFrames = interp1(meanAngles(s.calBeadIdx,goodFit_forward),...
                1:length(meanAngles(s.calBeadIdx,goodFit_forward)),linspace(-80,80,7),'nearest');       
s.templateIdxs = templateFrames;
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
%end compute template frames

s.smacmDarkFile = dataFile;
s.smacmDarkPath = dataPathOrig;
s.smacmSifFile = [];
s.smacmSifPath = [];
s.smacmEMGain = 1;

%compute templateLocs
templateSize = size(template,2);
templateFrames = templateFrames';
        numTemplates = size(templateFrames,1);
        templateColors = jet(numTemplates);
        templateLocs = zeros(numTemplates,5);
        fitParam = zeros(1,8);
        [xIdx, yIdx] = meshgrid(1:templateSize,1:templateSize);
        
%         hTemplate=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        for a=1:numTemplates
            
            % make minimum count level in template equal to 0
            template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                - min(min(template(templateFrames(a),:,:)));
            % normalize energy contained (sum of all counts) in the template
            template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                / sum(sum(template(templateFrames(a),:,:)));
            % finally, make mean of template equal to 0
            
                template(templateFrames(a),:,:) = template(templateFrames(a),:,:)...
                    - mean(mean(template(templateFrames(a),:,:)));
            

            % find two largest peaks in template
            [tempY, tempX] = ind2sub([templateSize templateSize],find(imregionalmax(template(templateFrames(a),:,:))));
            temp = sortrows([tempX tempY template(sub2ind(size(template),...
                templateFrames(a)*ones(length(tempX),1),tempY,tempX))],-3);
            
            % [amp1 amp2 xMean1 yMean1 xMean2 yMean2 sigma1 sigma2]
            fitParam(3) = temp(1,1);
            fitParam(4) = temp(1,2);
            fitParam(5) = temp(2,1);
            fitParam(6) = temp(2,2);
            fitParam(1) = temp(1,3);
            fitParam(2) = temp(2,3);
            fitParam(7) = mean(s.sigmaBounds);
            fitParam(8) = mean(s.sigmaBounds);
                  
            lowerBound = [0 0 1 1 1 1 s.sigmaBounds(1) s.sigmaBounds(1)];
            upperBound = [max(max(template(templateFrames(a),:,:))) max(max(template(templateFrames(a),:,:))) ...
                templateSize templateSize templateSize templateSize ...
                s.sigmaBounds(2) s.sigmaBounds(2)];

            % Fit with lsqnonlin
            options = optimset('FunValCheck','on','Diagnostics','off','Jacobian','on', 'Display', 'off');
            fitParam = lsqnonlin(@(x) ...
                f_doubleGaussianVector(x,squeeze(template(templateFrames(a),:,:)),0,xIdx,yIdx),...
                fitParam,lowerBound,upperBound,options);
            
            templateLocs(a,1:2) = fitParam(3:4);
            templateLocs(a,3:4) = fitParam(5:6);
            % calculate rough angle between peaks
            templateLocs(a,5) = 180/pi*atan2(templateLocs(a,2)-templateLocs(a,4), ...
                templateLocs(a,1)-templateLocs(a,3));
        end

s.templateLocs = templateLocs;
%end compute templateLocs

%threshold ouput variables
FOVmask = FOVmask;
FOVmask1 = FOVmask1;
usePolyROI = 1;
template = template;
    %median filter variables
    medianFilter = 1;
    windowSize = 101;
    medianBlurSigma = 15;
    dispRaw = true;
    interpVal = 10;
    if medianBlurSigma~=0
        medBlurFilt = fspecial('gaussian',100,medianBlurSigma);
    else
        medBlurFilt=1;
    end
    %end median filter variables

x0 = 900;
y0 = 900;

save('threshold output.mat','usePolyROI','FOVmask','FOVmask1','template',...
    'medianFilter','windowSize','medianBlurSigma','dispRaw','interpVal','medBlurFilt','x0','y0');
s.threshFilePrefix = cell(1);
s.threshFilePrefix{1} = [pwd '\'];
%end threshold output variables

s.nhaData = false;

% s.tempThresh = [645,645,645,645,645,645,645];
% s.tempThresh = [1500,1500,1500,1500,1500,1500,1500];
% s.tempThresh = [1000,1000,1000,1000,1000,1000,1000];
% s.tempThresh = 750*ones(1,7); %10000
s.tempThresh = 500*ones(1,7);




if isempty(s.tempThresh)
    s.templateThreshs = zeros(length(s.smacmRawFile),length(s.templateIdxs));
else
    s.templateThreshs = s.tempThresh;
end

s.templateThreshs=repmat(s.templateThreshs(s.threshFileSelIdx,:),...
    length(s.smacmRawFile),1);


                                    
%% Fit SMs
if isNoise
[s.fitFilePrefix,s.smacmSifFile,s.smacmSifPath] = f_fitSMs_V7_JR_SIM_Noise(s.smacmRawFile, s.smacmRawPath, ...
            [s.calFilePrefix 'calibration.mat'],s.calBeadIdx,...
            s.templateIdxs,s.templateThreshs/100000, s.smacmDarkFile, s.smacmDarkPath, s.smacmSifFile, s.smacmSifPath, s.boxRadius, ...
            s.channel,s.sigmaBounds, s.gaussianFilterSigma,s.minDistBetweenSMs,...
            s.lobeDistBounds,s.conversionGain,s.nmPerPixel,s.smacmEMGain,s.templateLocs,...
            [s.threshFilePrefix{1} 'threshold output.mat'],s.smacmRawROI,s.nhaData, s.fittingMethod); 
else
    [s.fitFilePrefix,s.smacmSifFile,s.smacmSifPath] = f_fitSMs_V7_JR_SIM_noNoise(s.smacmRawFile, s.smacmRawPath, ...
            [s.calFilePrefix 'calibration.mat'],s.calBeadIdx,...
            s.templateIdxs,s.templateThreshs/100000, s.smacmDarkFile, s.smacmDarkPath, s.smacmSifFile, s.smacmSifPath, s.boxRadius, ...
            s.channel,s.sigmaBounds, s.gaussianFilterSigma,s.minDistBetweenSMs,...
            s.lobeDistBounds,s.conversionGain,s.nmPerPixel,s.smacmEMGain,s.templateLocs,...
            [s.threshFilePrefix{1} 'threshold output.mat'],s.smacmRawROI,s.nhaData, s.fittingMethod);
end
        
        
%% Filter Output

s.useFids = false;

[totalPSFfits, numFrames, fidTrackX, fidTrackY, fidTrackZ,spatialCorr,useCurrent]...
    = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix,s.smacmSifFile, s.smacmSifPath, s.channel,[s.calFilePrefix 'calibration.mat'],s.calBeadIdx);

%4-6-2018
filename = strrep(dataFile,'.mat', '');

if isNoise == 1
    [savePath] = f_processFits_vCR_SIM(totalPSFfits,numFrames,s.fitFilePrefix,fidTrackX, fidTrackY, fidTrackZ, s.nmPerPixel,spatialCorr,useCurrent,s.calBeadIdx,dataPathOrig,filename);
else
    [savePath] = f_processFits_vCR_SIM_noNoise(totalPSFfits,numFrames,s.fitFilePrefix,fidTrackX, fidTrackY, fidTrackZ, s.nmPerPixel,spatialCorr,useCurrent,s.calBeadIdx,dataPathOrig,filename);
end
%% Post Processing
clear fiberData;
dataFileFiber = create_fiberData_auto(savePath);
% Tracking_MSD_auto;

end
toc;