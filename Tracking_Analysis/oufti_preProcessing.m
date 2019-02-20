%oufti_preProcessing
%This function inputs the localizations obtained by the Easy-DHPSF software
%and the corresponding phase contrast images and crops the phase contrast
%image to the approximate region where localizations are present.


%Select channel corresponding to the correct fluorescence pathway
%(red/green)
dlg_title = 'Input info';
prompt = {'Channel (r/g)'};
def = {'g'};
num_lines = 1;
inputInfo = inputdlg(prompt,dlg_title,num_lines,def);
channel = inputInfo{1};

%Multiple field'so-of-view (FOV) can be input and process simultaneously
dlg_title = 'Input';
prompt = {'Number of files?'};
def = {'1'};
num_lines = 1;
inputInfo = inputdlg(prompt,dlg_title,num_lines,def);
numFiles = str2num(inputInfo{1});

%Input files output from Easy-DHPSF for each FOV
dataFile1 = cell(numFiles,1);
dataPath1 = cell(numFiles,1);
for a = 1:numFiles
%load localizations file and phase image
[dataFile1{a}, dataPath1{a}] = uigetfile({'*.mat';'*.*'},'Open files with localizations');
end

%Input the corresponding phase contrast images for the FOVs above
[phaseFile1, phasePath1] = uigetfile({'*.tiff';'*.*'},'Open phase contrast images','MultiSelect', 'on',dataPath1{1});

%Loop over the FOVs
for a = 1:length(dataFile1)
    %Load the localizations
load ([dataPath1{a},dataFile1{a}]);

%Load the phase contrast image
if iscell(phaseFile1)
    phaseInfo = imfinfo([phasePath1 phaseFile1{a}]);
    phaseImg = zeros(phaseInfo(1).Height, phaseInfo(1).Width);
    phaseImg = phaseImg + double(imread([phasePath1 phaseFile1{a}],'Info', phaseInfo));
else
    phaseInfo = imfinfo([phasePath1 phaseFile1]);
    phaseImg = zeros(phaseInfo(1).Height, phaseInfo(1).Width);
    phaseImg = phaseImg + double(imread([phasePath1 phaseFile1],'Info', phaseInfo));
end


%set pixel size (nm) values for the phase contrast image and the
%fluorescence image
phasenmPerPixel = 110;
nmPerPixel = 108;


%resize and rotate phase image to match the localizations
phaseImg = imresize(phaseImg,phasenmPerPixel/nmPerPixel);
phaseHeight = size(phaseImg,1);
phaseWidth = size(phaseImg,2);
phaseImg = rot90(phaseImg); 

%flip the image to get it to match the orientation of the
%localizations
if channel == 'g'
        phaseImg = flipud(phaseImg);
end

%crop phase image only where localizations are present, and add buffer
%region so that no cells with localizations are accidentally cropped out of
%the phase image
xLim = round([(min(xLoc(:))-5000) (max(xLoc(:))+5000)]);
xLim1 = round(xLim(1,1)) ;
xLim2 = round(xLim(1,2));
yLim = round([(min(yLoc(:))-5000) (max(yLoc(:))+5000)]);
yLim1 = round(yLim(1,1));
yLim2 = round(yLim(1,2));
phaseMask = zeros(phaseWidth,phaseHeight);
[xPhase, yPhase] = meshgrid(0:phaseHeight,0:phaseWidth);
xRange = xPhase(1,:)*nmPerPixel-wlShiftX;
yRange = yPhase(:,1)*nmPerPixel-wlShiftY;
phaseMask(yRange > yLim1 & yRange < yLim2,xRange > xLim1 & xRange < xLim2 ) = 1;
[temp phaseXLim1] = find(xRange == min(xRange(xRange > xLim1)));
[temp phaseXLim2] = find(xRange == max(xRange(xRange < xLim2)));
[phaseYLim1 temp] = find(yRange == min(yRange(yRange > yLim1)));
[phaseYLim2 temp] = find(yRange == max(yRange(yRange < yLim2)));
phaseCrop = phaseImg(phaseYLim1:phaseYLim2,phaseXLim1:phaseXLim2);
phaseCrop = uint16(phaseCrop); 

%Plot the new phase image
figure;
imagesc(phaseCrop);
colormap gray;
axis square;

%Save the new phase image
if iscell(phaseFile1)
    filename = strsplit(phaseFile1{a},'_');
    filename = filename{1};
    save([dataPath1{a},'oufti_preProcess_' filename '.mat']);
    imwrite(phaseCrop,[dataPath1{a},'Cropped_',phaseFile1{a}]);
else
    filename = strsplit(phaseFile1,'_');
    filename = filename{1};
    save([dataPath1{1},'oufti_preProcess_' filename '1.mat']);
    imwrite(phaseCrop,[dataPath1{1},'Cropped1_',phaseFile1]);
end
end
clear all;