    
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

%%% TODO:
%%% add an 'import' function of the form below:

% loadNum =1;
% while loadNum ~= 0
% [~,tempPath] = uigetfile;
% if tempPath == 0
% break
% else
% r.fitFilePrefix{loadNum} = tempPath;
% loadNum = loadNum+1;
% end
% end

%%% also add an 'autosave'

%%% allow the user to rearrange threshold and fit files:
% tempRaw = g.smacmRawFile;
% for i = 1:9
% g.smacmRawFile{i} = tempRaw{i+1};
% end
% g.smacmRawFile{10} = tempRaw{1};

% tempFits = g.fitFilePrefix;
% for i = 1:7
% g.fitFilePrefix{i} = tempFits{i+1};
% end
% g.fitFilePrefix{8} = tempFits{1};



function easy_dhpsf_V2()
% easy_dhpsf allows scientific users to extract single-molecule
% localizations in 3D when using the double-helix point spread function
% widefield microscope. Tiff stacks of images are analyzed using template
% matching followed by double-Gaussian fitting to extract estimates of the
% molecule positions.
s.nhaData=0;


%% global file locations
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
s.calFilePrefix = '';
% DHPSF fiducial tracking data location
s.fidFilePrefix = {};
% DHPSF SM fit data location
s.fitFilePrefix = {};
% DHPSF threshold data save location
s.threshFilePrefix = {};
%% global fitting parameters
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
s.nmPerPixel = 108.;  % old value is 125.78; 8B back = 160
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
% s.boxRadius = 100;        % round(7*160/s.nmPerPixel);
% smoothing filter width for identifying DHPSF SMs, units of pixels
s.gaussianFilterSigma = 1.5*160/s.nmPerPixel;
% minimum lateral distance between identified SMs, units of pixels
s.minDistBetweenSMs = 7.5*160/s.nmPerPixel;
%
channelChoices = ['0';'r';'g'];
global r g
%% GUI parameters
figSize = [470 660];
figMargin = 20;
panelMargin = 10;
buttonWidth = 90;

%%  Initialize and hide the GUI as it is being constructed.
hfig = figure('Visible','off','Position',[1,1,figSize],...
    'MenuBar','None','ToolBar','none',...
    'DefaultUIPanelUnits','pixels','DefaultAxesUnits','pixels',...
    'DefaultUIPanelFontName','SegoeUI','DefaultUIPanelFontSize',10,...
    'DefaultUIControlFontName','SegoeUI','DefaultUIControlFontSize',9);

%% Construct the components.
htb = uitoolbar(hfig);
% upper master controls
hbuttonNew = uipushtool(htb,'CData',f_iconRead(fullfile(matlabroot,...
    'toolbox','matlab','icons','file_new.png')),...
    'TooltipString','New project',...
    'ClickedCallback',{@buttonNewProj_Callback});
hbuttonLoad = uipushtool(htb,'CData',f_iconRead(fullfile(matlabroot,...
    'toolbox','matlab','icons','file_open.png')),...
    'TooltipString','Load project',...
    'ClickedCallback',{@buttonLoadProj_Callback});
hbuttonSave = uipushtool(htb,'CData',f_iconRead(fullfile(matlabroot,...
    'toolbox','matlab','icons','file_save.png')),...
    'TooltipString','Save project',...
    'ClickedCallback',{@buttonSaveProj_Callback});
htextProjStatus = uicontrol('Style','text',...
    'Position',[figMargin,625,figSize(1)-2*figMargin,15]);
% individual panels for each module
hpanelSetup = uipanel('Title','Setup',...
    'Position',[figMargin 490 figSize(1)-2*figMargin 115]);
hpanelCal = uipanel('Title','Calibrate DHPSF',...
    'Position',[figMargin 405 figSize(1)-2*figMargin 65]);
hpanelThresh = uipanel('Title','Calibrate SM identification',...
    'Position',[figMargin 275 figSize(1)-2*figMargin 110]);
hpanelFid = uipanel('Title','Track fiduciaries',...
    'Position',[figMargin 190 figSize(1)-2*figMargin 65]);
hpanelFit = uipanel('Title','Localize DHPSF SMs',...
    'Position',[figMargin 105 figSize(1)-2*figMargin 65]);
hpanelOut = uipanel('Title','Output DHPSF SM localizations',...
    'Position',[figMargin figMargin figSize(1)-2*figMargin 65]);
% text boxes for statuses
htextCalStatus = uicontrol('Parent',hpanelCal,'Style','text',...
    'Position',[figSize(1)-120,panelMargin+5,70,15]);
htextFidStatus = uicontrol('Parent',hpanelFid,'Style','text',...
    'Position',[figSize(1)-120,panelMargin+5,70,15]);
htextThreshStatus = uicontrol('Parent',hpanelThresh,'Style','text',...
    'Position',[figSize(1)-120,panelMargin+5,70,15]);
htextFitStatus = uicontrol('Parent',hpanelFit,'Style','text',...
    'Position',[figSize(1)-120,panelMargin+5,70,15]);
% setup controls
% channel selection - see also the popupSetupChannel_Callback function
htextSetupChannel = uicontrol('Parent',hpanelSetup,'Style','text',...
    'String','Channel:',...
    'Position',[panelMargin,panelMargin+40,60,30],...
    'HorizontalAlignment','left');
hpopupSetupChannel = uicontrol('Parent',hpanelSetup,'Style','popupmenu',...
    'Position',[panelMargin+60,panelMargin+45,40,30],...
    'Callback',{@popupSetupChannel_Callback},...
    'String','0|R|G','Enable','on'); % TODO: move this into update code so that the string is set as R|G once using multicolor
htextSetupConv = uicontrol('Parent',hpanelSetup,'Style','text',...
    'String','Conversion gain:',...
    'Position',[panelMargin*2+100,panelMargin+50,70,30],...
    'HorizontalAlignment','left');
heditSetupConv = uicontrol('Parent',hpanelSetup,'Style','edit',...
    'String','26.93',...
    'Position',[panelMargin*2+170,panelMargin+50,50,20],...
    'HorizontalAlignment','left',...
    'Callback',{@editSetupConv_Callback});
htextSetupPixSize = uicontrol('Parent',hpanelSetup,'Style','text',...
    'String','Pixel size (nm):',...
    'Position',[panelMargin*3+220,panelMargin+50,60,30],...
    'HorizontalAlignment','left');
heditSetupPixSize = uicontrol('Parent',hpanelSetup,'Style','edit',...
    'String','120.1',...
    'Position',[panelMargin*3+280,panelMargin+50,50,20],...
    'HorizontalAlignment','left',...
    'Callback',{@editSetupPixSize_Callback});
htextSetupFittingMethod = uicontrol('Parent',hpanelSetup,'Style','text',...
    'String','Fitting Method:',...
    'Position',[panelMargin,panelMargin-5,100,30],...
    'HorizontalAlignment','left');
hpopupSetupFittingMethod = uicontrol('Parent',hpanelSetup,'Style','popupmenu',...
    'Position',[panelMargin+100,panelMargin,150,30],...
    'Callback',{@popupSetupFittingMethod_Callback},...
    'String',{'LSQ with DG model','MLE with DG model'},'Enable','on', ...
    'Value', 1); % TODO: move this into update code so that the string is set as R|G once using multicolor


% calibration controls
hbuttonCalRun = uicontrol('Parent',hpanelCal,'Style','pushbutton',...
    'String','Run','Position',[panelMargin,panelMargin,buttonWidth,25],...
    'Callback',{@buttonCalRun_Callback});
htextCalSel = uicontrol('Parent',hpanelCal,'Style','text',...
    'String','Calibration bead:',...
    'Position',[panelMargin*2+buttonWidth,panelMargin,70,30],...
    'HorizontalAlignment','left');
hpopupCalSel = uicontrol('Parent',hpanelCal,'Style','popupmenu',...
    'Position',[panelMargin*2+buttonWidth+70,panelMargin,50,30],...
    'Callback',{@popupCalSel_Callback});
haxesCalImg = axes('parent',hpanelCal,'position',[panelMargin*2+buttonWidth+130,panelMargin,45,45]);
% fiducial fitting controls
hbuttonFidRun = uicontrol('Parent',hpanelFid,'Style','pushbutton',...
    'String','Run','Position',[panelMargin,panelMargin,buttonWidth,25],...
    'Callback',{@buttonFidRun_Callback});
hcheckFidUse = uicontrol('Parent',hpanelFid,'Style','checkbox',...
    'String','Use fiduciaries','Min',0,'Max',1,...
    'Position',[panelMargin*2+buttonWidth,panelMargin,110,25],...
    'Callback',{@checkFidUse_Callback});
% thresholding controls
hbuttonThreshRun = uicontrol('Parent',hpanelThresh,'Style','pushbutton',...
    'String','Run','Position',[panelMargin,panelMargin+50,buttonWidth,25],...
    'Callback',{@buttonThreshRun_Callback});
htextThreshFileSel = uicontrol('Parent',hpanelThresh,'Style','text',...
    'String','Filename:',...
    'Position',[panelMargin*2+buttonWidth,panelMargin+50,70,30],...
    'HorizontalAlignment','left');
hpopupThreshFileSel = uicontrol('Parent',hpanelThresh,'Style','popupmenu',...
    'Position',[panelMargin*2+buttonWidth+70,panelMargin+50,210,30],...
    'Callback',{@popupThreshFileSel_Callback});
htextThreshSel = uicontrol('Parent',hpanelThresh,'Style','text',...
    'String','Template:',...
    'Position',[panelMargin*2+buttonWidth,panelMargin+20,70,30],...
    'HorizontalAlignment','left');
hpopupThreshSel = uicontrol('Parent',hpanelThresh,'Style','popupmenu',...
    'Position',[panelMargin*2+buttonWidth+70,panelMargin+20,50,30],...
    'Callback',{@popupThreshSel_Callback});
htextThreshVal = uicontrol('Parent',hpanelThresh,'Style','text',...
    'String','Threshold:',...
    'Position',[panelMargin*2+buttonWidth,panelMargin,70,20],...
    'HorizontalAlignment','left');
heditThreshVal = uicontrol('Parent',hpanelThresh,'Style','edit',...
    'String','0.000',...
    'Position',[panelMargin*2+buttonWidth+70,panelMargin,50,20],...
    'HorizontalAlignment','left',...
    'Callback',{@editThreshVal_Callback});
haxesThreshImg = axes('parent',hpanelThresh,'position',[panelMargin*2+buttonWidth+130,panelMargin,45,45]);
hbuttonThreshClone = uicontrol('Parent',hpanelThresh,'Style','pushbutton',...
    'String','Use for all','Position',[panelMargin,panelMargin,buttonWidth,25],...
    'Callback',{@buttonThreshClone_Callback});
% SM fitting controls
hbuttonFitRun = uicontrol('Parent',hpanelFit,'Style','pushbutton',...
    'String','Run','Position',[panelMargin,panelMargin,buttonWidth,25],...
    'Callback',{@buttonFitRun_Callback});
hbuttonFitDebug = uicontrol('Parent',hpanelFit,'Style','pushbutton',...
    'String','Debug','Position',[panelMargin*2+buttonWidth,panelMargin,buttonWidth,25],...
    'Callback',{@buttonFitDebug_Callback});
% data output controls
hbuttonOutExport = uicontrol('Parent',hpanelOut,'Style','pushbutton',...
    'String','Export to csv','Position',[panelMargin,panelMargin,buttonWidth,25],...
    'Callback',{@buttonOutExport_Callback});
hbuttonOutScatter = uicontrol('Parent',hpanelOut,'Style','pushbutton',...
    'String','3D scatterplot',...
    'Position',[panelMargin*2+buttonWidth,panelMargin,buttonWidth,25],...
    'Callback',{@buttonOutScatter_Callback});
hbuttonOutHist = uicontrol('Parent',hpanelOut,'Style','pushbutton',...
    'String','2D histogram','Position',...
    [panelMargin*3+buttonWidth*2,panelMargin,buttonWidth,25],...
    'Callback',{@buttonOutHist_Callback});
hbuttonOutReg = uicontrol('Parent',hpanelOut,'Style','pushbutton',...
    'String','Filter Output','Position',...
    [panelMargin*4+buttonWidth*3,panelMargin,buttonWidth-10,25],...
    'Callback',{@buttonOutReg_Callback});
% align([hpanelCal,hpanelFid,hpanelThresh,hpanelFit,hpanelOut],'Center','None');
newProj;

%% Initialize the GUI.
% Change units to normalized so components resize automatically.
set([hfig,htextProjStatus,...
    hpanelSetup,hpanelCal,hpanelFid,hpanelThresh,hpanelFit,hpanelOut,...
    htextCalStatus,htextFidStatus,htextThreshStatus,htextFitStatus,...
    hbuttonCalRun,hbuttonFidRun,hbuttonThreshRun,hbuttonThreshClone,...
    hbuttonFitRun,hbuttonFitDebug, hbuttonOutExport,hbuttonOutScatter,...
    hbuttonOutHist,hbuttonOutReg,htextSetupConv,heditSetupConv, htextSetupChannel,...
    hpopupSetupChannel,htextSetupFittingMethod,hpopupSetupFittingMethod,htextSetupPixSize, ...
    heditSetupPixSize,htextCalSel,...
    hpopupCalSel,haxesCalImg,hcheckFidUse,...
    htextThreshFileSel,hpopupThreshFileSel,...
    htextThreshSel,hpopupThreshSel,htextThreshVal,heditThreshVal,haxesThreshImg...
    ],'Units','normalized');

% Assign the GUI a name to appear in the window title.
set(hfig,'Name','Easy-DHPSF')
% initialize axes
axis(haxesCalImg,'off');
axis(haxesThreshImg,'off');
% Move the GUI to the center of the screen.
movegui(hfig,'center')
% Make the GUI visible.
set(hfig,'Visible','on');

    %% helper functions for setting the state of the program
    function newProj
        % confirm new project if current project is unsaved
        if s.projStatus(5) == 0 && any(s.projStatus(1:5))
            prompt = {'The current project is not saved. Continue?'};
            questiondialog = questdlg(prompt,'Confirm','Yes','No','Yes');
            % Handle response
            switch questiondialog
                case 'Yes'
                case 'No'
                    return
            end
        end
        projFile = '';
        s.smacmRawFile = {};
        s.smacmDarkFile = {};
        s.smacmSifFile = {};
        s.smacmLocFile = {};
        s.smacmFullLocFile = '';
        s.calFilePrefix = '';
        s.fidFilePrefix = {};
        s.fitFilePrefix = {};
        s.projStatus = false(1,5);
        s.templateImgs = [];
        s.numCalBeads = 1;
        s.calBeadIdx = 1;
        s.useFids = false;
        s.templateIdxs = [];
        s.templateThreshs = [];
        s.templateLocs = zeros(6,5);
        s.threshFileSelIdx = 1;
        s.templateSelIdx = 1;
        s.smacmRawROI = [0 0 0 0];
        s.fittingMethod = 'LSQ with DG model';
        s.channel = '0';
        s.tempThresh = [];
        s.smacmEMGain = 1;
        r = s;
        r.channel = 'r';
        g = s;
        g.channel = 'g';
        updateGUI;
    end
    function loadProj
        % confirm new project if current project is unsaved
        if s.projStatus(5) == 0 && any(s.projStatus(1:5))
            prompt = {'The current project is not saved. Continue?'};
            questiondialog = questdlg(prompt,'Confirm','Yes','No','Yes');
            % Handle response
            switch questiondialog
                case 'Yes'
                case 'No'
                    return
            end
        end
        [projFile, projPath] = uigetfile({'*.mat';'*.*'},'Open project MAT file');
        if isequal(projFile,0)
            projFile = '';
            return;
        end
        projFile = [projPath projFile];
        temp = load(projFile); % load into structure due to restriction on static workspaces
        s = temp.s;
        g = temp.g;
        r = temp.r;
        clear temp
        if ~isfield(s,'nhaData') % to be backwards-compatible
            s.nhaData = 0;
            r.nhaData = 0;
            g.nhaData = 0;
        end
        updateGUI;
    end
    function saveProj
        [tempFile, tempPath] = uiputfile({'*.mat';'*.*'},'Save project MAT file');
        if isequal(tempFile,0)
            clear tempFile tempPath
            return;
        end
        projFile = [tempPath tempFile];
        s.projStatus(5) = true;
        g.projStatus(5) = true;
        r.projStatus(5) = true;
        save(projFile,'s','r','g');
        clear tempFile tempPath
        updateGUI;
    end

    function updateGUI
        if s.channel == 'r'
            r = s;
        elseif s.channel == 'g'
            g = s;
        end
                
        %set(hpopupSetupFittingMethod,'Value',find(fittingMethodChoices==s.fittingMethod));
        set(hpopupSetupChannel,'Value',find(channelChoices==s.channel));
        set(heditSetupConv,'String',num2str(s.conversionGain));
        set(heditSetupPixSize,'String',num2str(s.nmPerPixel));
        if s.projStatus(1)
            set(htextCalStatus,'String','Complete','BackgroundColor','g');
            set(hbuttonThreshRun,'Enable','on');
            popupChoices = num2str(1:s.numCalBeads,'%g|');
            set(hpopupCalSel,'Enable','on',...
                'String',popupChoices(1:length(popupChoices)-1),...
                'Value',s.calBeadIdx);
            imagesc(squeeze(s.templateImgs(round(size(s.templateImgs,1)/2),:,:)),'parent',haxesCalImg);
            axis(haxesCalImg,'image');
            axis(haxesCalImg,'off');
            colormap(haxesCalImg,'hot');
        else
            set(htextCalStatus,'String','Incomplete','BackgroundColor','y');
            set(hbuttonFidRun,'Enable','off');
            set(hbuttonThreshRun,'Enable','off');
            set(hpopupCalSel,'Enable','off','String','1','Value',s.calBeadIdx);
            cla(haxesCalImg);
        end
        
        if s.projStatus(2)
            set(htextThreshStatus,'String','Complete','BackgroundColor','g');
            set(hbuttonFidRun,'Enable','on');
            set(hbuttonFitRun,'Enable','on');
            popupChoices = sprintf('%s|',s.smacmRawFile{:});
            set(hpopupThreshFileSel,'Enable','on',...
                'String',popupChoices(1:length(popupChoices)-1),...
                'Value',s.threshFileSelIdx);
            popupChoices = num2str(1:length(s.templateIdxs),'%g|');
            set(hpopupThreshSel,'Enable','on',...
                'String',popupChoices(1:length(popupChoices)-1),...
                'Value', s.templateSelIdx);
            set(heditThreshVal,'Enable','on','String',num2str(s.templateThreshs(...
                s.threshFileSelIdx,s.templateSelIdx)));
            imagesc(squeeze(s.templateImgs(s.templateIdxs(s.templateSelIdx),:,:)),...
                'parent',haxesThreshImg);
            axis(haxesThreshImg,'image');
            axis(haxesThreshImg,'off');
            colormap(haxesThreshImg,'hot');
        else
            set(htextThreshStatus,'String','Incomplete','BackgroundColor','y');
            set(hbuttonFitRun,'Enable','off');
            set(hpopupThreshFileSel,'Enable','off','String',' ','Value',1);
            set(hpopupThreshSel,'Enable','off','String','1','Value',s.templateSelIdx);
            set(heditThreshVal,'Enable','off','String','0.000','Value',0.000);
            cla(haxesThreshImg);
        end
        
        if  length(s.smacmRawFile)>1 && ~any(s.templateThreshs(s.threshFileSelIdx,:)==0)
            set(hbuttonThreshClone,'Enable','on');
        else
            set(hbuttonThreshClone,'Enable','off');
        end
        
        if s.projStatus(3)
            set(htextFidStatus,'String','Complete','BackgroundColor','g');
            set(hcheckFidUse,'Enable','on','Value',s.useFids);
        else
            set(htextFidStatus,'String','Incomplete','BackgroundColor','y');
            set(hcheckFidUse,'Enable','off');
        end
        if s.projStatus(4)
            set(htextFitStatus,'String','Complete','BackgroundColor','g');
            set(hbuttonOutExport,'Enable','on');
            set(hbuttonOutScatter,'Enable','on');
            set(hbuttonOutHist,'Enable','on');
            set(hbuttonFitDebug,'Enable','on');
            set(hbuttonOutReg,'Enable','on');
        else
            set(htextFitStatus,'String','Incomplete','BackgroundColor','y');
            set(hbuttonOutExport,'Enable','off');
            set(hbuttonOutScatter,'Enable','off');
            set(hbuttonOutHist,'Enable','off');
            set(hbuttonFitDebug,'Enable','off');
            set(hbuttonOutReg,'Enable','off');
        end
        if s.projStatus(5)
            set(htextProjStatus,'String',projFile,...
                'TooltipString',projFile,...
                'BackgroundColor',get(hfig,'Color'));
        else
            if isempty(projFile)
                set(htextProjStatus,'String','Unsaved project',...
                    'TooltipString','','BackgroundColor','y');
            else
                set(htextProjStatus,'String',['*' projFile],...
                    'TooltipString',projFile,'BackgroundColor','y');
            end
        end
    end

%% GUI function callbacks

    % master controls
    function buttonNewProj_Callback(~,~)
        newProj;
    end
    function buttonLoadProj_Callback(~,~)
        loadProj;
    end
    function buttonSaveProj_Callback(~,~)
        saveProj;
    end

    % setup controls
    function popupSetupChannel_Callback(source,~) 
        % if selecting a color after doing some analysis in the '0' channel
        % the user is given the option to move their analysis into their
        % channel of choice\
        % channelChoices = ['0';'r';'g']; (initialized at beginning)
        goMulticolor = 'No: Overwrite';
        if s.channel == '0' && any(s.projStatus)
            goMulticolor = questdlg('You already have non-labeled data. Do you want to import this into your chosen color?', ...
            'Switch to multicolor mode', ...
            'Yes','No: Overwrite','Cancel','Yes');
            if strcmp(goMulticolor,'Cancel')
                return
            end
        end
        s.channel = channelChoices(get(source,'Value'));
        if s.channel == 'r'
            if strcmp(goMulticolor,'Yes')
                r = s;
            else
                % the typical choice, for when switching between channels
                s = r;
            end
        elseif s.channel == 'g'
            if strcmp(goMulticolor,'Yes')
                g = s;
            else
                % the typical choice, for when switching between channels
                s = g;
            end
        end
        updateGUI;
    end
    

    function editSetupConv_Callback(source,~)
        s.conversionGain = str2double(get(source,'String'));
        s.projStatus(5) = false;
        updateGUI;
    end
    

    function editSetupPixSize_Callback(source,~)
        s.nmPerPixel = str2double(get(source,'String'));
        % update all dependent parameters
        s.sigmaBounds = [1.0 1.5]*160/s.nmPerPixel;
        s.lobeDistBounds = [3.5 10]*160/s.nmPerPixel;
        s.boxRadius = 13;        % round(7*160/s.nmPerPixel);
%         s.boxRadius = 100;        % round(7*160/s.nmPerPixel);
        s.gaussianFilterSigma = 1.5*160/s.nmPerPixel;
        s.minDistBetweenSMs = 7.5*160/s.nmPerPixel;
        s.projStatus(5) = false;
        updateGUI;
    end

    function popupSetupFittingMethod_Callback(source,~)
         sels = get(hpopupSetupFittingMethod,'String');
         idx  = get(hpopupSetupFittingMethod,'Value');
         s.fittingMethod =(sels{idx});
     
       % s.fittingMethod = ChannelChoices(get(source,'Value'));
        s.projStatus(5) = false;
        updateGUI;
    end

    % calibration controls
    function buttonCalRun_Callback(~,~)
        [s.calFilePrefix, s.numCalBeads] = ...
            f_calDHPSF_V6(s.conversionGain,s.nmPerPixel,s.boxRadius,s.channel,s.sigmaBounds,s.lobeDistBounds, s.fittingMethod);
%         [s.calFilePrefix, s.numCalBeads] = ...
%             f_calDHPSF_V5debug();
%         %             f_calDHPSF_V3(s.conversionGain,s.nmPerPixel,s.boxRadius,s.channel,s.sigmaBounds,s.lobeDistBounds);
        
        s.projStatus(1) = true;
        s.projStatus(5) = false;
        temp=load([s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],'template');
        s.templateImgs = temp.template;
        clear temp;
        set(hpopupSetupFittingMethod,'Enable','off');
        %set(hpopupSetupFittingMethod,'visible','off');  
        s.projStatus(5) = false;
        updateGUI;
    end
    function popupCalSel_Callback(source,~) 
        s.calBeadIdx = get(source,'Value');
        temp=load([s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],'template');
        s.templateImgs = temp.template;
        clear temp;
        
        s.projStatus(5) = false;
        updateGUI;
    end

    % threshold controls
    function buttonThreshRun_Callback(~,~) 
        [s.templateIdxs,s.smacmRawROI,s.smacmRawFile, s.smacmRawPath, ...
             s.smacmDarkFile,s.smacmDarkPath, s.smacmSifFile, s.smacmSifPath, s.smacmEMGain,...
             s.templateLocs,s.threshFilePrefix,s.nhaData,s.tempThresh] = f_calSMidentification_V6(...
             [s.calFilePrefix 'calibration.mat'],...
             s.calBeadIdx, [s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],...
            s.boxRadius,s.channel,s.sigmaBounds,s.gaussianFilterSigma,s.minDistBetweenSMs, s.fittingMethod);
        if isempty(s.tempThresh) 
            s.templateThreshs = zeros(length(s.smacmRawFile),length(s.templateIdxs));
        else
            s.templateThreshs = s.tempThresh;
        end
        
        s.projStatus(2) = true;
        s.projStatus(5) = false;
        updateGUI;     
    end
    function popupThreshFileSel_Callback(source,~) 
        s.threshFileSelIdx = get(source,'Value');
        updateGUI;
    end
    function popupThreshSel_Callback(source,~) 
        s.templateSelIdx = get(source,'Value');
        updateGUI;
    end
    function editThreshVal_Callback(source,~)
        s.templateThreshs(s.threshFileSelIdx,s.templateSelIdx)...
            = str2double(get(source,'String'));
        
        s.projStatus(5) = false;
        updateGUI;
    end
    function buttonThreshClone_Callback(~,~)
        s.templateThreshs=repmat(s.templateThreshs(s.threshFileSelIdx,:),...
                                        length(s.smacmRawFile),1);
        disp(['The threshold settings for ' s.smacmRawFile{s.threshFileSelIdx},...
               ' are now being used for all files.']);
        updateGUI;
    end

    % fiducial controls: gives the location of the raw fid fit files as a cell array
    function buttonFidRun_Callback(~,~) 
        [s.fidFilePrefix] = f_trackFiducials_V2(...
            s.smacmRawFile, s.smacmRawPath, [s.calFilePrefix 'calibration.mat'],s.calBeadIdx,...
            [s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],...
            s.templateIdxs, s.templateThreshs/100000, s.smacmDarkFile,s.smacmDarkPath, s.smacmSifFile, s.smacmSifPath, s.boxRadius, ...
            s.channel, s.gaussianFilterSigma,s.minDistBetweenSMs,...
            s.lobeDistBounds,s.conversionGain,s.nmPerPixel,s.smacmEMGain,s.templateLocs,s.sigmaBounds, s.fittingMethod);
%   [s.fidFilePrefix] = f_trackFiducials(...
%             s.smacmRawFile, s.smacmRawPath, [s.calFilePrefix 'calibration.mat'],s.calBeadIdx,...
%             [s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],...
%             s.templateIdxs, s.templateThreshs/100000, s.smacmDarkFile, s.smacmSifFile, s.smacmSifPath, s.boxRadius, ...
%             s.channel, s.gaussianFilterSigma,s.minDistBetweenSMs,...
%             s.lobeDistBounds,s.conversionGain,s.nmPerPixel,s.smacmEMGain,s.templateLocs,s.sigmaBounds);
        
        s.projStatus(3) = true;
        s.projStatus(5) = false;
        updateGUI;
    end
    function checkFidUse_Callback(source,~) 
        s.useFids = logical(get(source,'Value'));
        updateGUI;
    end
    
    % SM fitting controls
    function buttonFitRun_Callback(~,~)
        if  any(s.templateThreshs==0)
            msgbox(['One or more of the chosen template thresholds equals 0. '...
                    'Please define all thresholds before template matching.'],...
                    'Define thresholds', 'warn');
            return
        end
%         [s.fitFilePrefix,s.smacmSifFile,s.smacmSifPath] = f_fitSMs(s.smacmRawFile, s.smacmRawPath, ...
%             [s.calFilePrefix 'calibration.mat'],s.calBeadIdx,...
%             s.templateIdxs,s.templateThreshs/10000, s.smacmDarkFile, s.smacmSifFile, s.smacmSifPath, s.boxRadius, ...
%             s.channel,s.sigmaBounds, s.gaussianFilterSigma,s.minDistBetweenSMs,...
%             s.lobeDistBounds,s.conversionGain,s.nmPerPixel,s.smacmEMGain,s.templateLocs,...
%             [s.threshFilePrefix{1} 'threshold output.mat'],s.smacmRawROI,s.nhaData); 

           [s.fitFilePrefix,s.smacmSifFile,s.smacmSifPath] = f_fitSMs_V7(s.smacmRawFile, s.smacmRawPath, ...
            [s.calFilePrefix 'calibration.mat'],s.calBeadIdx,...
            s.templateIdxs,s.templateThreshs/100000, s.smacmDarkFile, s.smacmDarkPath, s.smacmSifFile, s.smacmSifPath, s.boxRadius, ...
            s.channel,s.sigmaBounds, s.gaussianFilterSigma,s.minDistBetweenSMs,...
            s.lobeDistBounds,s.conversionGain,s.nmPerPixel,s.smacmEMGain,s.templateLocs,...
            [s.threshFilePrefix{1} 'threshold output.mat'],s.smacmRawROI,s.nhaData, s.fittingMethod); 
        
         
        s.projStatus(4) = true;
        s.projStatus(5) = false;
        updateGUI;
    end
    function buttonFitDebug_Callback(~,~) 
        [totalPSFfits, numFrames] = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix,s.smacmSifFile, s.smacmSifPath, s.channel);
        if isfield(s,'threshFilePrefix') % backwards compatible before rev 56
        f_debugMoleculeFits(totalPSFfits,numFrames,s.smacmRawFile,...
            s.smacmRawPath,s.smacmDarkFile,s.smacmDarkPath,s.smacmRawROI,s.templateIdxs,...
            s.templateThreshs/100000,[s.threshFilePrefix{1} 'threshold output.mat'],...
            s.gaussianFilterSigma,s.minDistBetweenSMs)
        else
            f_debugMoleculeFits(totalPSFfits,numFrames,s.smacmRawFile,...
            s.smacmRawPath,s.smacmDarkFile,s.smacmDarkPath,s.smacmRawROI);
        end
        updateGUI;
    end
    % output controls
    function buttonOutExport_Callback(~,~) 
        [csvFile, csvPath] = uiputfile({'*.csv';'*.txt';'*.*'},'Save localizations as comma-separated file');
        if isequal(csvFile,0)
            return;
        end
        % makes sure that output has an extension
        if csvFile(end-3) ~= '.'
            csvFile = [csvFile '.csv'];
        end
        % open a file for writing
        [fid,message] = fopen([csvPath csvFile], 'w');
        if ~isempty(message)
            error([csvPath csvFile ': ' message]);
            %return;
        end
        % print a title, followed by a blank line
        fprintf(fid, ['frame num,molecule num,x (nm),y (nm),z (nm),' ...
            'x fid-corrected (nm),y fid-corrected (nm), z fid-corrected (nm),'...
            'photons detected,mean background photons\n']);
        fclose(fid);
        totalPSFfits = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix,s.smacmSifFile, s.smacmSifPath, s.channel);
        goodFit = totalPSFfits(:,17)>0;
        dlmwrite([csvPath csvFile],totalPSFfits(goodFit,[1 2 25 26 27 28 29 30 21 15]),...
            '-append');
        disp('Output file written successfully.');
        updateGUI;
    end
    function buttonOutScatter_Callback(~,~) 
        totalPSFfits = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix,s.smacmSifFile, s.smacmSifPath, s.channel,[s.calFilePrefix 'calibration.mat']);
        f_scatter3(totalPSFfits,s.useFids);
        
        updateGUI;
    end
    function buttonOutHist_Callback(~,~) 
        totalPSFfits = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix,s.smacmSifFile, s.smacmSifPath, s.channel,[s.calFilePrefix 'calibration.mat']);
        f_hist2(totalPSFfits,s.useFids);
        
        updateGUI;
    end
    function buttonOutReg_Callback(~,~)
        % fidTrackX,fidTrackY,fidTrackZ are used to keep fid data separate
        % for 2-color processing. That is, when recombining the channels,
        % the uncorrected localizations are transformed, the tracks are
        % transformed, and finally the transformed localizations and tracks
        % are subtracted to perform fiducial correction. You must still check
        % the 'use fiducials' box (as of rev 48).
        [totalPSFfits, numFrames, fidTrackX, fidTrackY, fidTrackZ,spatialCorr,useCurrent]...
            = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix,s.smacmSifFile, s.smacmSifPath, s.channel,[s.calFilePrefix 'calibration.mat'],s.calBeadIdx);
%         [matFile, matPath] = uiputfile({'*.mat';'*.*'},'Save localizations as old-style .mat file');
%         if isequal(matFile,0)
%             return;
%         end
        %conversionFactor = s.conversionGain /s.smacmEMGain;
        %save([matPath matFile], 'totalPSFfits','numFrames','ROI','conversionFactor');
        
        % these match the numbers in f_fitSMs
%         ampRatioLimit = 0.5;
%         sigmaRatioLimit = 0.4;
f_processFits_vCR(totalPSFfits,numFrames,s.fitFilePrefix,fidTrackX, fidTrackY, fidTrackZ, s.nmPerPixel,spatialCorr,useCurrent,s.calBeadIdx);
%         f_processFits(totalPSFfits,numFrames,ROI,conversionFactor,...
%             s.sigmaBounds,s.lobeDistBounds,ampRatioLimit,sigmaRatioLimit,...
%             s.nmPerPixel);
    end
end
