function Select_Cell()
% select data from a filtered fits (output from scatter3_wl) .mat file and
% a WL image. Saved to a structure. Used to be used for CreS fibers, and
% terminology of variables / dialog boxes will reflect this.
transmitted = 0; % i.e., is it untransformed PAmCherry data?

useTimeColors = 0;
scatterSize = 30;
wlShiftX = 0;
wlShiftY = 0;

% nmPerPixel = 125.78; 
nmPerPixel = 108; 


scrsz = get(0,'ScreenSize');

%% open datafiles

[locFile locPath] = uigetfile({'*.mat';'*.*'},'Open data file with filtered localizations');
if isequal(locFile,0)
    error('User cancelled the program');
end
load([locPath locFile]);
% [whiteLightFile whiteLightPath] = uigetfile({'*.tif';'*.*'},'Open image stack with white light image');
[whiteLightFile whiteLightPath] = uigetfile({'*.tiff';'*.*'},'Open image stack with white light image');


% The uigetfile below allows you to add fibers to an existing structure file.
% The new figures and structure file are saved in a new folder.
% You must copy over the figures and fits if you want to collect everything 
% in the new folder.
[fiberFile fiberPath] = uigetfile({'*.mat';'*.*'},'To add fibers to an existing structure containing CreS fiber data, open it now');



% dealing with a structure array, comment this when not doing so
% xLoc = croppedDataSets(2).xLoc;
% yLoc = croppedDataSets(2).yLoc;
% zLoc = croppedDataSets(2).zLoc;
% sigmaX = croppedDataSets(2).sigmaX;
% sigmaY = croppedDataSets(2).sigmaY;
% sigmaZ = croppedDataSets(2).sigmaZ;
% frameNum = croppedDataSets(2).frameNum;
% numPhotons = croppedDataSets(2).numPhotons;
% meanBkgnd = croppedDataSets(2).meanBkgnd;

% xLoc, etc. are later specified to be the points inside the ROIs
xLocFull = xLoc;
yLocFull = yLoc;
zLocFull = zLoc_IndexCorrected;
sigmaXFull = sigmaX;
sigmaYFull = sigmaY;
sigmaZFull = sigmaZ;
frameNumFull = frameNum;
numPhotonsFull = numPhotons;
meanBkgndFull = meanBkgnd;

locName = locFile(1:length(locFile)-4);
fiberData(1).Name = locName;

% puts all data one directory below the image
outputFilePrefix = [locPath locName ' Fibers ' ...
    datestr(now,'yyyymmdd HHMM') filesep];

mkdir(outputFilePrefix);


%% Choose a desired parameter set for reconstruction
    
    dlg_title = 'Please Input Parameters';
    prompt = {  'Size of points in reconstruction',...
        'Temporal Color Coding',...
        'White light shift X (in nm)',...
        'White light shift Y (in nm)',...
        };
    def = { ...
        num2str(scatterSize), ...
        num2str(useTimeColors), ...
        num2str(wlShiftX), ...
        num2str(wlShiftY), ...
        };
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    
    scatterSize = str2double(inputdialog{1});
    useTimeColors = str2double(inputdialog{2});
    wlShiftX = str2double(inputdialog{3});
    wlShiftY = str2double(inputdialog{4});
    
    anotherfiber = true;
    fiberNum = 1;

% if a previous selection was specified, this sets the counter to begin
% at the number following the last fiber picked.
if ~isequal(fiberFile,0)
    load ([fiberPath fiberFile]);
    ROICenterX = [fiberData.ROICenterX];
    ROICenterY = [fiberData.ROICenterY];
    fiberNum = length(ROICenterX) + 1;
end
    
if transmitted == 1
    ROI_WL = [256, 256, 256, 256];
else
    ROI_WL = [1, 1, 270, 270];
end
     % [x0, y0, xlength, ylength]
     
while anotherfiber == true    
    %% Plot the white light image if specified
    
    close all
   
    if whiteLightFile ~= 0
%         whiteLightInfo = imfinfo([whiteLightPath whiteLightFile]);
        whiteLightInfo = imfinfo([whiteLightFile]);
        
        whiteLight = zeros(whiteLightInfo(1).Height, whiteLightInfo(1).Width);
        % average white light images together to get better SNR
        for a = 1:length(whiteLightInfo)
            whiteLight = whiteLight + double(imread([whiteLightPath whiteLightFile], ...
                'Info', whiteLightInfo));
        end
        % resize white light to the size of the ROI of the single molecule fits
        whiteLight = whiteLight(ROI_WL(2):ROI_WL(2)+ROI_WL(4)-1,ROI_WL(1):ROI_WL(1)+ROI_WL(3)-1);
        % rescale white light image to vary from 0 to 1
        whiteLight = (whiteLight-min(whiteLight(:)))/(max(whiteLight(:))-min(whiteLight(:)));
        [xWL yWL] = meshgrid((ROI_WL(1):ROI_WL(1)+ROI_WL(3)-1) * nmPerPixel + wlShiftX, ...
            (ROI_WL(2):ROI_WL(2)+ROI_WL(4)-1) * nmPerPixel + wlShiftY);
        %         [xWL yWL] = meshgrid((1:(whiteLightInfo(1).Width)) * nmPerPixel + wlShiftX, ...
        %         (1:(whiteLightInfo(1).Height)) * nmPerPixel + wlShiftY);
    end
    
    %% ask user what region to plot in superresolution image
    
    if whiteLightFile ~= 0
        xRange = xWL(1,:);
        yRange = yWL(:,1);
        % pick region that contains background
        figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(xRange,yRange,whiteLight);axis image;colormap gray;
        hold on;
    end
    
    figure('Position', get(0,'Screensize'));
    %scatter3(xLoc,yLoc,zLoc,1,'filled');
    scatter(xLocFull,yLocFull,1,'filled');
    xlim([min(xLocFull(:)) max(xLocFull(:))]);
    ylim([min(yLocFull(:)) max(yLocFull(:))]);
    xlabel('x (nm)');ylabel('y (nm)');
    axis ij;
    axis equal;
    
    if fiberNum > 1
        hold on;
        for a = 1:length(ROICenterX)
            text(ROICenterX(a),ROICenterY(a),num2str(a),'color','red','fontsize',13,'fontweight','bold');
        end
    end    
    hold off;

%     ROI = imrect(gca,[min(xLoc(:)) min(yLoc(:)) max(xLoc(:))-min(xLoc(:)) max(yLoc(:))-min(yLoc(:))]);
      ROI = imrect(gca,[min(xLoc(:)) min(yLoc(:)) 20000 20000]);

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
    
    inROI = xLocFull>=ROI(1) & xLocFull<=ROI(1)+ROI(3) & yLocFull>ROI(2) & yLocFull<=ROI(2)+ROI(4) & numPhotonsFull>0;
    %invalidPoints = xLoc_bad>=ROI(1) & xLoc_bad<=ROI(1)+ROI(3) & yLoc_bad>ROI(2) & yLoc_bad<=ROI(2)+ROI(4) ;
    
    %add IMPOLY
    figure('Position', get(0,'Screensize'));
    scatter(xLocFull,yLocFull,1,'filled');
    xlim([ROI(1) ROI(1)+ROI(3)]);
    ylim([ROI(2) ROI(2)+ROI(4)]);
    xlabel('x (nm)');ylabel('y (nm)');
    axis ij;
    axis equal;
    
    h = impoly(gca);
    fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),...
        get(gca,'YLim'));
    setPositionConstraintFcn(h,fcn);
    h1 = wait(h);
    
    in = inpolygon(xLocFull,yLocFull,h1(:,1),h1(:,2));
    inROI = in;
    %end add IMPOLY
    
    %xLocPix = xLocPix(validPoints);
    %yLocPix = yLocPix(validPoints); 
    xLoc = xLocFull(inROI);
    yLoc = yLocFull(inROI);
    zLoc = zLocFull(inROI);
    numPhotons = numPhotonsFull(inROI);
    meanBkgnd = meanBkgndFull(inROI);
    frameNum = frameNumFull(inROI);
    sigmaX = sigmaXFull(inROI);
    sigmaY = sigmaYFull(inROI);
    sigmaZ = sigmaZFull(inROI);
    
    meanNumPhotons = mean(numPhotons);
    
    
    %% plot 3D scatterplot of localizations with white light
    
    f = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','k','renderer','opengl', 'Toolbar', 'figure');
%     f = figure('Position', get(0,'Screensize'));     
    % draws white light image below localizations
    if whiteLightFile~=0
        %imagesc(xRange,yRange,whiteLight);axis image;colormap gray;hold on;
        [x,y,z] = meshgrid(xRange,yRange,[min(zLoc) max(zLoc)]);
        xslice = []; yslice = []; zslice = min(zLoc);
        h=slice(x,y,z,repmat(whiteLight,[1 1 2]),xslice,yslice,zslice,'nearest');
        set(h,'EdgeColor','none','FaceAlpha',0.75);
        colormap gray; hold on;
    end
    
    if useTimeColors == 0
        scatter3(xLoc,yLoc,zLoc,scatterSize,[1 1 0],'filled');
    else
        %         scatter3(xLoc(a),yLoc(a),zLoc(a),scatterSize,frameNum(1):frameNum(length(frameNum)),'filled')
        markerColors = jet(frameNum(length(frameNum))-frameNum(1)+1);
        %     for a = 1:length(frameNum)
        %         scatter3(xLoc(a)-min(xLoc(:)),yLoc(a)-min(yLoc(:)),zLoc(a),scatterSize,'filled',...
        %             'MarkerFaceColor', markerColors(frameNum(a)-frameNum(1)+1,:),...
        %             'MarkerEdgeColor', markerColors(frameNum(a)-frameNum(1)+1,:));
        %     end
        
        for a = 1:length(frameNum)
            scatter3(xLoc(a),yLoc(a),zLoc(a),scatterSize,'filled',...
                'MarkerFaceColor', markerColors(frameNum(a)-frameNum(1)+1,:),...
                'MarkerEdgeColor', markerColors(frameNum(a)-frameNum(1)+1,:));
        end
    end
    axis vis3d equal;
    %     xlim([min(xLoc(:)) max(xLoc(:))]);
    %     ylim([min(yLoc(:)) max(yLoc(:))]);
    %     xlim([min(xLoc(:)) max(xLoc(:))]-min(xLoc(:)));
    %     ylim([min(yLoc(:)) max(yLoc(:))]-min(yLoc(:)));
    xlim([min(xLoc(:)) max(xLoc(:))]);
    ylim([min(yLoc(:)) max(yLoc(:))]);
    xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');

    title({[num2str(length(xLoc)) ' localizations'];...
        ['Mean Number of Signal Photons = ' num2str(meanNumPhotons) ' per integration time'];...
        ['Mean Number of Background Photons = ' num2str(mean(meanBkgnd)) ' per pixel per integration time'];...
        ['Localization Precision \sigma_x = ' num2str(mean(sigmaX)) ' nm'];...
        ['Localization Precision \sigma_y = ' num2str(mean(sigmaY)) ' nm'];...
        ['Localization Precision \sigma_z = ' num2str(mean(sigmaZ)) ' nm']},...
        'color','w');
    set(gca,'color','k');
    set(gca,'xcolor','w');set(gca,'ycolor','w');set(gca,'zcolor','w');
    
    % waits until the user presses 'Continue'
    h = uicontrol('Position',[20 20 200 40],'String','Continue',...
            'Callback','uiresume(gcbf)');
        uiwait(gcf);
    %% Write stats into structure and possibly restart loop with new fiber

        % Queries whether the fiber should be saved
     dlg_title = 'Acceptable fiber?';
    prompt = {'Would you like to save that fiber?'};
    def =       { 'Yes'  };
    questiondialog = questdlg(prompt,dlg_title, def);
    
    % Handle response
    switch questiondialog
        case 'Yes'
    % saves all info into fiberData and updates ROI markers   
    fiberData(fiberNum).xLoc = xLoc;
    fiberData(fiberNum).yLoc = yLoc;
    fiberData(fiberNum).zLoc = zLoc;
    fiberData(fiberNum).sigmaX = sigmaX;
    fiberData(fiberNum).sigmaY = sigmaY;
    fiberData(fiberNum).sigmaZ = sigmaZ;
    fiberData(fiberNum).numPhotons = numPhotons;
%     ROICenterX(fiberNum) = ROI(1)+ROI(3)/2;
%     ROICenterY(fiberNum) = ROI(2)+ROI(4)/2;
    ROICenterX(fiberNum) = (min(xLoc(:)) + max(xLoc(:)))/2;
    ROICenterY(fiberNum) = (min(yLoc(:)) + max(yLoc(:)))/2;
    
     
    fiberData(fiberNum).ROICenterX = ROICenterX(fiberNum);
    fiberData(fiberNum).ROICenterY = ROICenterY(fiberNum);
    fiberData(fiberNum).Frame = frameNum;
   % saves images from each fiber

    saveas(gcf,[outputFilePrefix locName ' Fiber ' num2str(fiberNum) ' 3D.fig']);
    
    fiberNum = fiberNum + 1;
    
        case 'No'
        case 'Cancel'
            error('User cancelled the program');
    end    
    
        % Queries whether to select a new fiber
    dlg_title = 'New fiber?';
    prompt = {'Would you like to choose another fiber?'};
    def =       { 'Yes'  };
    questiondialog = questdlg(prompt,dlg_title, def);
    
    % Handle response
    switch questiondialog
        case 'Yes'
        case 'No'
            anotherfiber = false;
        case 'Cancel'
            error('User cancelled the program');
    end
end
%% Saves all data
   close all
   
   % draws figure with all selections
    if whiteLightFile ~= 0
        xRange = xWL(1,:);
        yRange = yWL(:,1);
        % pick region that contains background
        figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(xRange,yRange,whiteLight);axis image;colormap gray;
        hold on;
    end

    scatter(xLocFull,yLocFull,1,'filled');
    xlim([min(xLocFull(:)) max(xLocFull(:))]);
    ylim([min(yLocFull(:)) max(yLocFull(:))]);
    xlabel('x (nm)');ylabel('y (nm)');
    axis ij;
    hold on;
        for a = 1:length(ROICenterX)
            text(ROICenterX(a),ROICenterY(a),num2str(a),'color','red','fontsize',13,'fontweight','bold');
        end  
    hold off;
    
    fiberData(1).FrameRange = frameRange;
    
    saveas(gcf,[outputFilePrefix locName ' Fiber selections.fig']);
    save([outputFilePrefix locName ' Fiber Data'],'fiberData'); 

    disp('The fiber selection program exited successfully.');
    disp(['The saved fibers are in the folder ' outputFilePrefix]);
if ~isequal(fiberFile,0)
    disp(['You added fibers to ' fiberPath fiberFile]);

end
end
%%
% output excel spreadsheet
% textHeader = {'frame number' ...  
%     'fiduciary corrected x location (nm)' ...
%     'fiduciary corrected y location (nm)' ...
%     'fiduciary corrected z location (nm)' ...
%     'sigma x (nm)' ...
%     'sigma y (nm)' ...
%     'sigma z (nm)' ...
%     'number of photons' ...
%     'mean background photons' };
% output = [frameNum, xLoc, yLoc, zLoc, sigmaX, sigmaY, sigmaZ, numPhotons, meanBkgnd];
% xlswrite([savePath saveFile(1:length(saveFile)-4) '.xlsx'], [textHeader; ...
%     num2cell(output)], ...
%     'valid PSF fits');

