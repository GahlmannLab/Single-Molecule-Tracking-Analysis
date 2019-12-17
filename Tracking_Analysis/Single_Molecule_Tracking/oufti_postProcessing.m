% Julian Rocha

%Oufti_postProcessing
% This function inputs the information output from oufti_preProcessing, the
% cell outlines output from Oufti, and the localizations obtained from
% Easy_DHPSF and precisely overlays the cell outlines with the localizations

%Input the relevent files ouput from outfti_preProcessing, Oufti, and
%Easy_DHPSF
[dataFile2, dataPath2] = uigetfile({'*.mat';'*.*'},'Open oufti pre process info');
load ([dataPath2,dataFile2]);

[dataFile, dataPath] = uigetfile({'*.mat';'*.*'},'Open file with oufti cell outlines');
load ([dataPath,dataFile]);

[dataFile1, dataPath1] = uigetfile({'*.mat';'*.*'},'Open file with localizations');
load ([dataPath1,dataFile1]);


%crop and shift white light file from pre process file
cropShiftX = phaseXLim1*nmPerPixel;
cropShiftY = phaseYLim1*nmPerPixel;


meshXmin = 1*10^10;
meshXmax = -1*10^10;
meshYmin = 1*10^10;
meshYmax = -1*10^10;
meshPts = [];

%Shift the positions of the outlines until they are closely (but not
%perfectly) overlaid on the localizations. The code will ask the user if
%the data is well aligned after shifting. If it is not, the process is
%repeated until the user is satisfied with the overlay.
aligned = 'No';
inputval = [-25000; -1600];
while strcmp(aligned,'No')
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    %Reset index
    idx = 0;
    shiftval = {'Enter Y shift value (nm)', 'Enter X shift value (nm)'};
    titleval = 'Input';
    dimval = [1 35];
    definputval = {num2str(inputval(1,1)), num2str(inputval(2,1))};
    inputval = inputdlg(shiftval,titleval,dimval,definputval);
    inputval = str2double(inputval);
    
    outlineShiftY = inputval(1,1);
    outlineShiftX = inputval(2,1);
    
    for a = 1:cellListN
        mesh = cellList.meshData{1,1}{1,a}.mesh;
        if mesh ~= 0
            mesh = mesh*phasenmPerPixel;
            mesh(:,1) = mesh(:,1) + cropShiftX + outlineShiftX;
            mesh(:,2) = mesh(:,2) + cropShiftY + outlineShiftY;
            mesh(:,3) = mesh(:,3) + cropShiftX + outlineShiftX;
            mesh(:,4) = mesh(:,4) + cropShiftY + outlineShiftY;
            
            Pts = double(unique(vertcat(mesh(:,1:2),flipud(mesh(:,3:4))),'rows','stable'));
            
            hold on;
            plot(mesh(:,1),mesh(:,2),'g',mesh(:,3),mesh(:,4),'g');
            if any(~isfinite(Pts))
                continue;
            end
            idx = idx + 1;
            meshPts{idx,1} = Pts;
            meshXmin = min(min(meshPts{idx,1}(:,1)), meshXmin);
            meshXmax = max(max(meshPts{idx,1}(:,1)), meshXmax);
            meshYmin = min(min(meshPts{idx,1}(:,2)), meshYmin);
            meshYmax = max(max(meshPts{idx,1}(:,2)), meshYmax);
            
            
        end
    end
    scatter(xLoc(:),yLoc(:),10,'.','r');
    xlim([min(min(xLoc(:)),meshXmin)-500 max(max(xLoc(:)),meshXmax)+500]);
    ylim([min(min(yLoc(:)),meshYmin)-500 max(max(yLoc(:)),meshYmax)+500]);
    xlabel('x (nm)');ylabel('y (nm)');
    axis ij;
    axis square;
    
    %Initiate iterations. If the data is not well aligned return the the
    %screen to select new shift values.
    quest = {'Is the data well aligned?'};
    title1 = 'Input';
    btn1 = 'Yes';
    btn2 = 'No';
    defbtn = btn2;
    aligned = questdlg(quest,title1,btn1,btn2,defbtn);

end    
meshPts_trans = []; 


%This section will ask the user to pick 5 control points in both the cell
%outline and the localization data sets. The user should pick the same
%relative 5 points in both data sets (i.e. The cell poles of the same 5
%cells should be selected in both data sets. If the initial try does not
%produce a well aligned transformation, repeat the procedure and choose
%different points. The 5 points should be equally sampled throughout the
%field of view (4 corners and the center).
alignedFine = 'No';
while strcmp(alignedFine,'No');
    figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    for a = 1:cellListN
        mesh = cellList.meshData{1,1}{1,a}.mesh;
        if mesh ~= 0
            mesh = mesh*phasenmPerPixel;
            mesh(:,1) = mesh(:,1) + cropShiftX + outlineShiftX;
            mesh(:,2) = mesh(:,2) + cropShiftY + outlineShiftY;
            mesh(:,3) = mesh(:,3) + cropShiftX + outlineShiftX;
            mesh(:,4) = mesh(:,4) + cropShiftY + outlineShiftY;
            
            Pts = double(unique(vertcat(mesh(:,1:2),flipud(mesh(:,3:4))),'rows','stable'));
            
            %plot the cell outlines
            hold on;
            plot(mesh(:,1),mesh(:,2),'g',mesh(:,3),mesh(:,4),'g');
            if any(~isfinite(Pts))
                continue;
            end
            idx = idx + 1;
            meshPts{idx,1} = Pts;
            meshXmin = min(min(meshPts{idx,1}(:,1)), meshXmin);
            meshXmax = max(max(meshPts{idx,1}(:,1)), meshXmax);
            meshYmin = min(min(meshPts{idx,1}(:,2)), meshYmin);
            meshYmax = max(max(meshPts{idx,1}(:,2)), meshYmax);
            
            
        end
    end
    %plot the localizations
    scatter(xLoc(:),yLoc(:),10,'.','r');
    xlim([min(min(xLoc(:)),meshXmin)-500 max(max(xLoc(:)),meshXmax)+500]);
    ylim([min(min(yLoc(:)),meshYmin)-500 max(max(yLoc(:)),meshYmax)+500]);
    xlabel('x (nm)');ylabel('y (nm)');
    axis ij;
    axis square;
    
    cp_channel1 = [];
    n = 0;
    but = 1;
    title([{'Select 5 points for rough transformation at cell poles (4 corners and center).';...
        'Select outline points first (purple), then hit enter.'}]);
    while but == 1
        [xi,yi,but] = ginput(1);
        if isempty(xi)
            break
        end
        n = n+1;
        text(xi,yi,num2str(n),'color','m','fontsize',13,'fontweight','bold');
        cp_channel1(n,:) = round([xi yi]);
    end

    hold on;
    cp_channel2 = [];
    n = 0;
    but = 1;
    title({'Now select 5 corresponding localization points (blue), then hit enter.'});
    while but == 1
        [xi,yi,but] = ginput(1);
        if isempty(xi)
            break
        end
        n = n+1;
        text(xi,yi,num2str(n),'color','blue','fontsize',13,'fontweight','bold');
        cp_channel2(n,:) = round([xi yi]);
    end

    %calculate the rough tform from the 5 points and plot again. This will
    %produce a rough overlay of the cell outlines and localizations
    [tform, FRE, TRE, FRE_full, TRE_full] = matlab_transformation(...
        cp_channel1, cp_channel2 , 'affine');

        figure
        set(gcf, 'Position', get(0, 'Screensize'));
        title('Select points inside of cells to delete, then hit enter');
        scatter(xLoc(:),yLoc(:),10,'.','r');
        hold on;
        meshPts_trans = [];
        for a = 1:length(meshPts)
            meshPts_trans{a,1} = tformfwd(meshPts{a,1},tform);
            if ~isempty(meshPts_trans{a,1})
            scatter(meshPts_trans{a,1}(:,1),meshPts_trans{a,1}(:,2),'.','g');
            end
        end
        meshPts_trans = meshPts_trans(~cellfun('isempty',meshPts_trans));
        hold on;

        xlabel('x (nm)');ylabel('y (nm)');
        axis ij;
        axis square; 
        
        
    %Initiate iterations. If the data is not well aligned return the the
    %screen to select new shift values. In this step the alignment must be
    %much more precise than in the previous rough overlay.
    quest = {'Is the data well aligned?'};
    title1 = 'Input';
    btn1 = 'Yes';
    btn2 = 'No';
    defbtn = btn2;
    alignedFine = questdlg(quest,title1,btn1,btn2,defbtn);
   
end

%Plot the overlay and delete unwanted cells. Cells at the edge of the
%field-of-view (FOV) and cells with unusually low/high numbers of
%localizations should be deleted as well.
figure; hold on;
set(gcf, 'Position', get(0, 'Screensize'));
title('Click points inside of cells to delete, then hit enter');
scatter(xLoc(:),yLoc(:),10,'.','r');
meshPts_trans = [];
for a = 1:length(meshPts)
    meshPts_trans{a,1} = tformfwd(meshPts{a,1},tform);
    %delete cells with less than 10 points inside. This reduces the need to manually delete cells later on. 
    ptsIdx = [];
    ptsIdx = inpolygon(xLoc(:,1),yLoc(:,1),meshPts_trans{a,1}(:,1),meshPts_trans{a,1}(:,2));
    if sum(ptsIdx) <= 10
        meshPts_trans{a,1} = [];
    end
    if ~isempty(meshPts_trans{a,1})
        scatter(meshPts_trans{a,1}(:,1),meshPts_trans{a,1}(:,2),'.','g');
    end
end
meshPts_trans = meshPts_trans(~cellfun('isempty',meshPts_trans));
xlabel('x (nm)');ylabel('y (nm)');
axis ij;
axis square;

%select cells to delete by picking point inside cell
delPts = [];
n = 0;
but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    if isempty(xi)
        break
    end
    n = n+1;
    text(xi,yi,num2str(n),'color','black','fontsize',13,'fontweight','bold');
    delPts(n,:) = [xi yi];
end

deleteIdxTotal = [];
for a = 1:length(meshPts_trans)
   deletePts = [];
   deletePts = inpolygon(delPts(:,1),delPts(:,2),meshPts_trans{a,1}(:,1),meshPts_trans{a,1}(:,2));
   if sum(deletePts) ~= 0
       meshPts_trans{a,1} = [];
   end
end

%Find the center of the outlines and the center of the localizations inside the cell outlines
%for a finer transformation, then calculate new tform. This produces a more
%finely sampled transformation with many more control points.
meshPts_trans = meshPts_trans(~cellfun('isempty',meshPts_trans));
LocCenter = [];
OutlineCenter = [];
for a = 1:length(meshPts_trans)
    ptsIdx = inpolygon(xLoc(:,1),yLoc(:,1),meshPts_trans{a,1}(:,1),meshPts_trans{a,1}(:,2));
    LocCenter(a,1) = mean(xLoc(ptsIdx,1));
    LocCenter(a,2) = mean(yLoc(ptsIdx,1));
    OutlineCenter(a,1) = mean(meshPts_trans{a,1}(:,1));
    OutlineCenter(a,2) = mean(meshPts_trans{a,1}(:,2));
end
OutlineCenter = OutlineCenter(~any(isnan(LocCenter),2),:);
LocCenter = LocCenter(~any(isnan(LocCenter),2),:);
if length(OutlineCenter(:,1) <= 3)
[tform2, FRE2, TRE2, FRE_full2, TRE_full2] = matlab_transformation(...
    OutlineCenter,LocCenter, 'nonreflective similarity');
else
[tform2, FRE2, TRE2, FRE_full2, TRE_full2] = matlab_transformation(...
    OutlineCenter,LocCenter, 'affine');
end

%Plot the points used for the transformation
figure;
hold on;
meshPts_trans2 = [];
for a = 1:length(meshPts_trans)
meshPts_trans2{a,1} = tformfwd(meshPts_trans{a,1},tform2);
scatter(meshPts_trans2{a,1}(:,1),meshPts_trans2{a,1}(:,2),'.','g');
end
hold on;
scatter(xLoc(:),yLoc(:),10,'.','r');
xlim([min(xLoc(:))-500 max(xLoc(:))+500]);
ylim([min(yLoc(:))-500 max(yLoc(:))+500]);
xlabel('x (nm)');ylabel('y (nm)');
axis ij;
axis square;


%Plot final figure and delete all localizations not located within a cell
%outline
figure; 
hold on; 
for a = 1:length(meshPts_trans)
meshPts_trans2{a,1} = tformfwd(meshPts_trans{a,1},tform2);
scatter(meshPts_trans2{a,1}(:,1),meshPts_trans2{a,1}(:,2),'.','g');
end
for a = 1:length(meshPts_trans2)
    ptsIdx = inpolygon(xLoc(:,1),yLoc(:,1),meshPts_trans2{a,1}(:,1),meshPts_trans2{a,1}(:,2));
    scatter(xLoc(ptsIdx,1),yLoc(ptsIdx,1),'.','r');
 end
axis ij;
axis square;

fiberData1 = [];
for a = 1:length(meshPts_trans2)
    ptsIdx = inpolygon(xLoc(:,1),yLoc(:,1),meshPts_trans2{a,1}(:,1),meshPts_trans2{a,1}(:,2));
    fiberData1(a).xLoc = xLoc(ptsIdx,1);
    fiberData1(a).yLoc = yLoc(ptsIdx,1);
%     fiberData1(a).zLoc = zLoc(ptsIdx,1);
    fiberData1(a).zLoc = zLoc_IndexCorrected(ptsIdx,1);
    fiberData1(a).sigmaX = sigmaX(ptsIdx,1);
    fiberData1(a).sigmaY = sigmaY(ptsIdx,1);
    fiberData1(a).sigmaZ = sigmaZ(ptsIdx,1);
    fiberData1(a).numPhotons = numPhotons(ptsIdx,1);
    fiberData1(a).meanBkgnd = meanBkgnd(ptsIdx,1);
    fiberData1(a).Frame = frameNum(ptsIdx,1);
    fiberData1(a).meshPtsX = meshPts_trans2{a,1}(:,1);
    fiberData1(a).meshPtsY = meshPts_trans2{a,1}(:,2);
    
    fiberData = fiberData1(1,1); 
    c = 1;
    for b = 1:length(fiberData1)
        if ~isempty(fiberData1(b).xLoc)
            fiberData(c) = fiberData1(b);
            c = c+1;
        end
    end
end
fiberData(1).FrameRange = frameRange;
filename = regexprep(filename, '\.[^\.]*$', '');

%save fiberData variable with localizations for each cell
outputFilePrefix = [filename,' Oufti Fibers '];
save([dataPath1 outputFilePrefix ' Fiber Data'],'fiberData','filename'); 

clear all;
close all;
