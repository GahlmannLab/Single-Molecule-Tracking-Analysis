
function shut_pos_gui(data)
%% Initialization
% MCLNanoDriveStage;
MCLNanoDriveStage.ReleaseAllHandles; % remove handles from nanodrive
MCLMicroDriveStage.ReleaseAllHandles;
stage2 = MCLMicroDriveStage;
stage = MCLNanoDriveStage; % create instance of MCLNanoDriveStage
[xM, yM] = stage2.position;
[x,y,z] = stage.GetPosition;


% initialize daq channels
g = daq.createSession('ni');
addDigitalChannel(g, 'Dev1', 'port1/line1','InputOnly'); % [DI,tt]=inputSingleScan(s); create DI(1), DI(2), DI(3)
addDigitalChannel(g, 'Dev1', 'port1/line2','InputOnly');
addDigitalChannel(g, 'Dev1', 'port1/line3','InputOnly');
s1 = daq.createSession('ni');
addDigitalChannel(s1, 'Dev1', 'port0/line0','OutputOnly');
s2 = daq.createSession('ni');
addDigitalChannel(s2, 'Dev1', 'port0/line1','OutputOnly');
s3 = daq.createSession('ni');
addDigitalChannel(s3, 'Dev1', 'port0/line2','OutputOnly');
f1 = daq.createSession('ni');
addDigitalChannel(f1, 'Dev1', 'port0/line3','OutputOnly');
l1 = daq.createSession('ni');
addDigitalChannel(l1, 'Dev1', 'port0/line5','OutputOnly');

o1 = daq.createSession('ni'); % create session for OBIS laser
addDigitalChannel(o1, 'Dev1', 'port0/line6','OutputOnly'); % add a digital channel to send signal to OBIS laser

s = daq.createSession('ni');
addAnalogOutputChannel(s,'Dev1','ao0','Voltage'); % both cameras are queued to the same session (use [a, a] rather than a when triggering)
addAnalogOutputChannel(s,'Dev1','ao1','Voltage');

% fields
ignore = 0;
abort = 0;
FolderName='Z:\Eugene\';
e=0;
llimit = 0;
ulimit = 200;
mllimit = -4000;
mulimit = 4000;
FileName='';
step = 30;
step_size = 50;
frame = 10;
exposure = 0.03;
delayT = 15;
[DI,tt]=inputSingleScan(g);
t= 0.001;
xyStep = 7;
zStep = 5;
xyStepSize = 0.5;
zStepSize = 0.3;
gridFrame = 10;
gridExposure = 50;

% OBIS_exp = 50;
% OBIS_waiting_off = 10;

% % Test run only GUI (w/o functions implemented)
% x=10;
% y=10;
% z=10;
% DI = [1, 1, 1];

%% create the GUI figure
guiFig = figure('Units','pixels','Position',[350 350 350 700],...
    'MenuBar','none','ToolBar','none','Visible','on',...
    'NumberTitle','off', 'Name', 'Shutter/Stage GUI');
defaultBackground = get(0,'defaultUicontrolBackgroundColor');
set(guiFig,'Color',defaultBackground)
% create a handle
handles.output = guiFig;
% create panels for nano-positioner and micro-positioner
panelPosition = [.01 .82 .27 .13];
handles.XYPan = uipanel('Parent',guiFig,'Units','normalized','fontweight','bold',...
    'Position',panelPosition,'Title','XY-Pos','Visible','on');
panelPosition = [.32 .82 .16 .13];
handles.ZPan = uipanel('Parent',guiFig,'Units','normalized','fontweight','bold',...
    'Position',panelPosition,'Title','Z-Pos','Visible','on');
panelPosition = [.01 .63 .47 .15];
handles.Pos = uipanel('Parent',guiFig,'Units','normalized','fontweight','bold',...
    'Position',panelPosition,'Title','Nano-Position','Visible','on');
panelPosition = [.62 .82 .27 .13];
handles.XYPan = uipanel('Parent',guiFig,'Units','normalized','fontweight','bold',...
    'Position',panelPosition,'Title','XY-Pos','Visible','on');
panelPosition = [.52 .63 .47 .15];
handles.Pos = uipanel('Parent',guiFig,'Units','normalized','fontweight','bold',...
    'Position',panelPosition,'Title','Micro-Position','Visible','on');

% draw lines
ax = axes;
set(ax, 'Visible', 'off');
set(ax, 'Position', [0, 0, 1, 1]); %# Stretch the axes over the whole figure.
set(ax, 'Xlim', [0, 1], 'YLim', [0, 1]); %# Switch off autoscaling.

line([0.5, 0.5], [0.585, 1], 'Parent', ax,'Color',[.8 .8 .8]); %# Draw!
line([1, 0], [0.585, 0.585], 'Parent', ax,'Color',[.8 .8 .8]);
line([1, 0], [0.42, 0.42], 'Parent', ax,'Color',[.8 .8 .8]);
line([1, 0], [0.29, 0.29], 'Parent', ax,'Color',[.8 .8 .8]);
line([1, 0], [0.145, 0.145], 'Parent', ax,'Color',[.8 .8 .8]);


% NP label
handles.NP_txt = uicontrol('Units', 'normalized','FontWeight', 'bold', 'Style','text',...
    'String','Nano-Positioning Stage', 'Position',[0.05,0.955,0.4,0.03]);
% MP label
handles.MP_txt = uicontrol('Units', 'normalized','FontWeight', 'bold', 'Style','text',...
    'String','Micro-Positioning Stage', 'Position',[0.55,0.955,0.4,0.03]);
% Z-Cal Label
handles.ZCal_txt = uicontrol('Units', 'normalized','FontWeight', 'bold', 'Style','text',...
    'String','Z-Calibration', 'Position',[0,0.385,1,0.027]);
% Grid Label
handles.Grid_txt = uicontrol('Units', 'normalized','FontWeight', 'bold', 'Style','text',...
    'String','Control-Point Calibration', 'Position',[0,0.255,1,0.027]);

% XY-DIRECTIONS NP
handles.x_LEFT    = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','<','FontWeight', 'bold', 'FontSize', 11, 'Position',[0.04,0.86,0.05,0.04],'CallBack',@xleft_call);
handles.y_UP    = uicontrol('Units','normalized', 'Style','pushbutton',...
    'String','^','FontSize', 16, 'Position',[0.11,0.885,.05,.04],'CallBack',@yup_call);
handles.x_RIGHT    = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','>','FontWeight', 'bold','FontSize', 11, 'Position',[0.18,0.86,0.054,0.04],'CallBack',@xright_call);
handles.y_DOWN    = uicontrol('Units', 'normalized', 'Style','pushbutton',...
    'String','V','FontWeight', 'bold','Position',[0.11,0.835,0.05,0.04],'CallBack',@ydown_call);

% XY-DIRECTIONS MP
handles.x_LEFT1    = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','<','FontWeight', 'bold', 'FontSize', 11, 'Position',[0.65,0.86,0.05,0.04],'CallBack',@xleft1_call);
handles.y_UP1   = uicontrol('Units','normalized', 'Style','pushbutton',...
    'String','^','FontSize', 16, 'Position',[0.72,0.885,.05,.04],'CallBack',@yup1_call);
handles.x_RIGHT1    = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','>','FontWeight', 'bold','FontSize', 11, 'Position',[0.79,0.86,0.054,0.04],'CallBack',@xright1_call);
handles.y_DOWN1   = uicontrol('Units', 'normalized', 'Style','pushbutton',...
    'String','V','FontWeight', 'bold','Position',[0.72,0.835,0.05,0.04],'CallBack',@ydown1_call);

% Z-DIRECTIONS NP
handles.z_UP    = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','^', 'FontSize', 16, 'Position',[0.375,0.885,0.05,0.04],'CallBack',@zup_call);
handles.z_DOWN    = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','V','FontWeight', 'bold', 'Position',[0.375,0.835,0.05,0.04],'CallBack',@zdown_call);

% STEP-SIZE INPUT NP
handles.xy_val = uicontrol('Units', 'normalized','Style','popup',...
    'Position',[0.06,0.77,0.16,0.05],'String',{'20 nm' '0.1 um' '0.5 um','1.0 um','2.0 um','5.0 um','10 um'});
handles.z_val = uicontrol('Units', 'normalized','Style','popup',...
    'Position',[0.32,0.77,0.16,0.05],'String',{'.1 um' '0.5 um','1.0 um','2.0 um','5.0 um','10 um'});

% STEP-SIZE INPUT MP
handles.xy1_val = uicontrol('Units', 'normalized','Style','popup',...
    'Position',[0.67,0.77,0.16,0.05],'String',{'10 um', '20 um', '50 um', '100 um'});

% NANO-POSITIONING STAGE CURRENT POSITION
handles.xpos_txt = uicontrol('Units', 'normalized','Style','text',...
    'String',['X:     ' num2str(x,'%.3f') ' um'], 'Position',[0.05,0.72,0.35,0.03]);
handles.ypos_txt = uicontrol('Units', 'normalized','Style','text',...
    'String',['Y:     ' num2str(y,'%.3f') ' um'], 'Position',[0.05,0.68,0.35,0.03]);
handles.zpos_txt = uicontrol('Units', 'normalized','Style','text',...
    'String',['Z:     ' num2str(z,'%.3f') ' um'], 'Position',[0.05,0.64,0.35,0.03]);

% MICRO-POSITIONING STAGE CURRENT POSITION
handles.xpos_txt1 = uicontrol('Units', 'normalized','Style','text',...
    'String',['X:     ' num2str(xM,'%.3f') ' um'], 'Position',[0.55,0.70,0.35,0.03]);
handles.ypos_txt1 = uicontrol('Units', 'normalized','Style','text',...
    'String',['Y:     ' num2str(yM,'%.3f') ' um'], 'Position',[0.55,0.66,0.35,0.03]);


% REFRESH POSITION NP
handles.refresh_but = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','Refresh Position', 'Position',[0.01,0.59,0.47,0.04],'CallBack',@refresh_call);

% REFRESH POSITION MP
handles.refresh1_but = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','Refresh Position', 'Position',[0.52,0.59,0.47,0.04],'CallBack',@refresh1_call);

% PHASE CONTRAST/MANUAL MODE BUTTONS
handles.manual_but = uicontrol('Units', 'normalized','Style','checkbox',...
    'String','Manual Mode', 'fontweight', 'bold','enable', 'off','value',1, 'Position',[0.02,0.53, 0.3,0.04],'CallBack',@manual_call);
handles.phasec_but = uicontrol('Units', 'normalized','Style','checkbox',...
    'String','Phase-Contrast', 'fontweight', 'bold','enable','off','Position',[0.02,0.48,0.3,0.04],'CallBack',@phasec_call);

% TOGGLE/CHECK BUTTONS *Sync Active must be on for shutters

% if(DI(1) == 1)
%     handles.shut1_chck = uicontrol('Units', 'normalized','Style','togglebutton',...
%         'String','Open 405nm Laser', 'Position',[0.345,0.53,0.3,0.04],'CallBack',@shut1_call);
% else
%     handles.shut1_chck = uicontrol('Units', 'normalized','Style','togglebutton',...
%         'String','Close 405nm Laser', 'Position',[0.345,0.53,0.3,0.04],'CallBack',@shut1_call);
% end
    % Manually control 405 nm Laser
    handles.Ctrl_OBIS_laser = uicontrol('Units', 'normalized','Style','togglebutton',...
        'String','Open 405nm Laser', 'Position',[0.345,0.53,0.3,0.04],'CallBack',@Ctrl_OBIS_laser);
%     % Initialize and set 405
%     handles.Set_OBIS_laser = uicontrol('Units', 'normalized','Style','pushbutton',...
%         'String','Set 405nm Laser', 'Position',[0.68,0.53,0.3,0.04],'CallBack',@Set_OBIS_laser_0);    
        
if(DI(2) == 1)
    handles.shut2_chck = uicontrol('Units', 'normalized','Style','togglebutton',...
        'String','Open 515nm Laser', 'Position',[0.345,0.48,0.3,0.04],'CallBack',@shut2_call);
else
    handles.shut2_chck = uicontrol('Units', 'normalized','Style','togglebutton',...
        'String','Close 515nm Laser', 'Position',[0.345,0.48,0.3,0.04],'CallBack',@shut2_call);
end
if(DI(3) == 1)
    handles.shut3_chck = uicontrol('Units', 'normalized','Style','togglebutton',...
        'String','Open 561nm Laser', 'Position',[0.345,0.43,0.3,0.04],'CallBack',@shut3_call);
else
    handles.shut3_chck = uicontrol('Units', 'normalized','Style','togglebutton',...
        'String','Close 561nm Laser', 'Position',[0.345,0.43,0.3,0.04],'CallBack',@shut3_call);
end
handles.fmirror_chck = uicontrol('Units', 'normalized','Style','togglebutton',...
    'String','Toggle Mirror Up', 'Position',[0.68,0.5,0.3,0.04],'CallBack',@fmir_call);
handles.LED_chck = uicontrol('Units', 'normalized','Style','checkbox',...
    'String','Thorlabs LED', 'Position',[0.68,0.45,0.3,0.05],'CallBack',@led_call);
handles.imgsrc_but = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','Imaging Source', 'Position',[0.02,0.43,0.3,0.04],'CallBack',@imgsrc_call);

% Z-cal
handles.zcal_step    = uicontrol('Units', 'normalized','Style','edit',...
    'String','30', 'Position',[0.58,0.35,0.20,0.03],'CallBack',@zstep_call);
handles.zcal_size   = uicontrol('Units','normalized', 'Style','edit',...
    'String','50', 'Position',[0.26,0.35,.20,.03],'CallBack',@zsize_call);
handles.zcal_frame    = uicontrol('Units', 'normalized','Style','edit',...
    'String','10', 'Position',[0.58,0.305,0.2,0.03],'CallBack',@zframe_call);
handles.zcal_exptime   = uicontrol('Units', 'normalized', 'Style','edit',...
    'String','0.03','Position',[0.26,0.305,0.2,0.03],'CallBack',@ztime_call);

% Z-cal text
handles.Zsize_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','        Step-Size:', 'Position',[0.01,0.345,0.23,0.03]);
handles.Ztime_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','Exposure Time:', 'Position',[0.01,0.3,0.23,0.03]);
handles.Zsize_unit_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','nm   ', 'Position',[0.46,0.345,0.08,0.03]);
handles.Ztime_unit_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','secs', 'Position',[0.46,0.3,0.08,0.03]);
handles.Zstep_unit_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','steps  ', 'Position',[0.78,0.345,0.11,0.03]);
handles.Zframe_unit_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','frames', 'Position',[0.78,0.3,0.11,0.03]);

% Grid Text
handles.xyStepSize_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','XY Step-Size:   ', 'Position',[0.02,0.22,0.23,0.03]);
handles.zStepSize_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','Z Step-Size:      ', 'Position',[0.02,0.185,0.23,0.03]);
handles.gridExposure_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','Exposure Time:  ', 'Position',[0.02,0.15,0.23,0.03]);
handles.xyStep_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','XY Steps:  ', 'Position',[0.52,0.22,0.23,0.03]);
handles.zStep_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','Z Steps:     ', 'Position',[0.52,0.185,0.23,0.03]);
handles.gridFrame_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','Frames:      ', 'Position',[0.52,0.15,0.23,0.03]);
% Grid
handles.xyStepSize_val = uicontrol('Units', 'normalized', 'Style','edit',...
    'String','0.5', 'Position',[0.26,0.225,0.15,0.03],'CallBack',@xyStepSize_call);
handles.zStepSize_val = uicontrol('Units', 'normalized', 'Style','edit',...
    'String','0.3', 'Position',[0.26,0.19,0.15,0.03],'CallBack',@zStepSize_call);
handles.gridExposure_val = uicontrol('Units', 'normalized', 'Style','edit',...
    'String','50', 'Position',[0.26,0.155,0.15,0.03],'CallBack',@gridExposure_call);
handles.xyStep_val    = uicontrol('Units', 'normalized','Style','edit',...
    'String','7', 'Position',[0.72,0.225,0.15,0.03],'CallBack',@xyStep_call);
handles.zStep_val   = uicontrol('Units','normalized', 'Style','edit',...
    'String','5', 'Position',[0.72,0.19,.15,.03],'CallBack',@zStep_call);
handles.gridFrame_val = uicontrol('Units', 'normalized', 'Style','edit',...
    'String','10', 'Position',[0.72,0.155,0.15,0.03],'CallBack',@gridFrame_call);

% Grid Units
handles.xyStep_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','um ', 'Position',[0.415,0.22,0.06,0.03]);
handles.zStep_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','um ', 'Position',[0.415,0.185,0.06,0.03]);
handles.gridFrame_txt = uicontrol('Units', 'normalized', 'Style','text',...
    'String','sec', 'Position',[0.415,0.15,0.06,0.03]);

% BUTTONS
handles.grid_but = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','GO', 'Position',[0.91,0.18,0.07,0.05],'CallBack',@grid_call);
handles.zcal_but = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','GO', 'Position',[0.91,0.315,0.07,0.05],'CallBack',@zcal_call);
handles.loadfile_but = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','LOAD', 'Position',[0.04,0.10,0.2,0.04],'CallBack',@load_call);
handles.start_but = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','START', 'Position',[0.28,0.10,0.2,0.04],'CallBack',@start_call);
handles.abort_but = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','ABORT', 'Position',[0.52,0.10,0.2,0.04],'CallBack',@abort_call);
handles.exit_but = uicontrol('Units', 'normalized','Style','pushbutton',...
    'String','EXIT', 'Position',[0.76,0.10,0.2,0.04],'CallBack',@exit_call);

% STATUS OF SEQUENCE
handles.file_txt = uicontrol('Units', 'normalized','Style','text',...
    'String',['FILE LOADED: ' 'No File Loaded'], 'Position',[0.00,0.015,1,0.07]);
handles.sequence_txt = uicontrol('Units', 'normalized','Style','text',...
    'String',['SEQUENCE STATUS: ' 'Sequence Not Running'], 'fontweight', 'normal', 'Position',[0.00,0.01,1,0.03]);

%% timer for refreshing position
tmr = timer('Name','Reminder',...
    'Period',1,... % Update the time every 60 seconds.
    'StartDelay',0,...
    'TasksToExecute',inf,... % number of times to update
    'ExecutionMode','fixedSpacing',...
    'TimerFcn',{@updater});
start(tmr);
set(guiFig, 'deletefcn', {@deleter})
    function updater(fileTab, eventdata)
        [x,y,z] = stage.GetPosition;
        set(handles.xpos_txt, 'string', ['X:     ' num2str(x,'%.3f') ' um']);
        set(handles.ypos_txt, 'string', ['Y:     ' num2str(y,'%.3f') ' um']);
        set(handles.zpos_txt, 'string', ['Z:     ' num2str(z,'%.3f') ' um']);
        [xM,yM] = stage2.position;
        set(handles.xpos_txt1, 'string', ['X:     ' num2str(xM,'%.3f') ' um']);
        set(handles.ypos_txt1, 'string', ['Y:     ' num2str(yM,'%.3f') ' um']);
    end
    function deleter(fileTab, eventdata)
        stop(tmr);
        delete(tmr);
    end
%% Nanocyte directional control (nano-positioner)
    function xleft1_call(fileTab,eventdata)
        str = get(handles.xy1_val,'String');
        val = get(handles.xy1_val,'Value'); % retrieve selected step-size
        switch str{val};
            case '10 um'
                [xM, yM] = stage2.position;
                if(xM-10<mllimit) % if step goes beyond limits, return
                    return
                end
                stage2.Left(10);
                pause(0.1);
            case '20 um'
                [xM, yM] = stage2.position;
                if(xM-20<mllimit)
                    return
                end
                stage2.Left(20);
                pause(0.1);
            case '50 um'
                [xM, yM] = stage2.position;
                if(xM-50<mllimit)
                    return
                end
                stage2.Left(50);
                pause(0.1);
            case '100 um'
                [xM, yM] = stage2.position;
                if(xM-100<mllimit)
                    return
                end
                stage2.Left(100);
                pause(0.1);
        end
    end
    function xright1_call(fileTab,eventdata)
        str = get(handles.xy1_val,'String');
        val = get(handles.xy1_val,'Value');
        switch str{val};
            case '10 um'
                [xM, yM] = stage2.position;
                if(xM+10>mulimit)
                    return
                end
                stage2.Right(10);
                pause(0.1);
            case '20 um'
                [xM, yM] = stage2.position;
                if(xM+20>mulimit)
                    return
                end
                stage2.Right(20);
                pause(0.1);
            case '50 um'
                [xM, yM] = stage2.position;
                if(xM+50>mulimit)
                    return
                end
                stage2.Right(50);
                pause(0.1);
            case '100 um'
                [xM, yM] = stage2.position;
                if(xM+100>mulimit)
                    return
                end
                stage2.Right(100);
                pause(0.1);
        end
    end
    function yup1_call(fileTab,eventdata)
        str = get(handles.xy1_val,'String');
        val = get(handles.xy1_val,'Value');
        switch str{val};
            case '10 um'
                [xM, yM] = stage2.position;
                if(yM+10>mulimit)
                    return
                end
                stage2.Down(10);
                pause(0.1);
            case '20 um'
                [xM, yM] = stage2.position;
                if(yM+20>mulimit)
                    return
                end
                stage2.Down(20);
                pause(0.1);
            case '50 um'
                [xM, yM] = stage2.position;
                if(yM+50>mulimit)
                    return
                end
                stage2.Down(50);
                pause(0.1);
            case '100 um'
                [xM, yM] = stage2.position;
                if(yM+100>mulimit)
                    return
                end
                stage2.Down(100);
                pause(0.1);
        end
    end
    function ydown1_call(fileTab,eventdata)
        str = get(handles.xy1_val,'String');
        val = get(handles.xy1_val,'Value');
        switch str{val};
            case '10 um'
                [xM, yM] = stage2.position;
                if(yM-10<mllimit)
                    return
                end
                stage2.Up(10);
                pause(0.1);
            case '20 um'
                [xM, yM] = stage2.position;
                if(yM-20<mllimit)
                    return
                end
                stage2.Up(20);
                pause(0.1);
            case '50 um'
                [xM, yM] = stage2.position;
                if(yM-50<mllimit)
                    return
                end
                stage2.Up(50);
                pause(0.1);
            case '100 um'
                [xM, yM] = stage2.position;
                if(yM-100<mllimit)
                    return
                end
                stage2.Up(100);
                pause(0.1);
        end
    end
    function xleft_call(fileTab,eventdata)
        str = get(handles.xy_val,'String');
        val = get(handles.xy_val,'Value');
        switch str{val};
            case '20 nm'
                if(x-0.02<llimit)
                    return
                end
                x = x - 0.02;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.1 um'
                if(x-0.1<llimit)
                    return
                end
                x = x - 0.1;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.5 um'
                if(x-0.5<llimit)
                    return
                end
                x = x - 0.5;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '1.0 um'
                if(x-1<llimit)
                    return
                end
                x = x - 1.0;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '2.0 um'
                if(x-2<llimit)
                    return
                end
                x = x - 2;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '5.0 um'
                if(x-5<llimit)
                    return
                end
                x = x - 5;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '10 um'
                if(x-10<llimit)
                    return
                end
                x = x - 10;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
        end
    end
    function yup_call(fileTab,eventdata)
        str = get(handles.xy_val, 'String');
        val = get(handles.xy_val,'Value');
        switch str{val};
            case '20 nm'
                if(y+0.02>ulimit)
                    return
                end
                y = y + 0.02;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.1 um'
                if(y+0.1>ulimit)
                    return
                end
                y = y + 0.1;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.5 um'
                if(y+0.5>ulimit)
                    return
                end
                y = y + 0.5;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '1.0 um'
                if(y+1>ulimit)
                    return
                end
                y = y + 1.0;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '2.0 um'
                if(y+2>ulimit)
                    return
                end
                y = y + 2;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '5.0 um'
                if(y+5>ulimit)
                    return
                end
                y = y + 5;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '10 um'
                if(y+10>ulimit)
                    return
                end
                y = y + 10;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
        end
    end
    function xright_call(fileTab,eventdata)
        str = get(handles.xy_val, 'String');
        val = get(handles.xy_val,'Value');
        switch str{val};
            case '20 nm'
                if(x+0.02>ulimit)
                    return
                end
                x = x + 0.02;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.1 um'
                if(x+0.1>ulimit)
                    return
                end
                x = x + 0.1;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.5 um'
                if(x+0.5>ulimit)
                    return
                end
                x = x + 0.5;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '1.0 um'
                if(x+1>ulimit)
                    return
                end
                x = x + 1.0;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '2.0 um'
                if(x+2>ulimit)
                    return
                end
                x = x + 2;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '5.0 um'
                if(x+5>ulimit)
                    return
                end
                x = x + 5;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '10 um'
                if(x+10>ulimit)
                    return
                end
                x = x + 10;
                stage.SetX(x);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
        end
    end
    function ydown_call(fileTab,eventdata)
        str = get(handles.xy_val, 'String');
        val = get(handles.xy_val,'Value');
        switch str{val};
            case '20 nm'
                if(y-0.02<llimit)
                    return
                end
                y = y - 0.02;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.1 um'
                if(y-0.1<llimit)
                    return
                end
                y = y - 0.1;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.5 um'
                if(y-0.5<llimit)
                    return
                end
                y = y  - 0.5;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '1.0 um'
                if(y-1<llimit)
                    return
                end
                y = y  - 1.0;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '2.0 um'
                if(y-2<llimit)
                    return
                end
                y = y - 2;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '5.0 um'
                if(y-5<llimit)
                    return
                end
                y = y - 5;
                stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '10 um'
                if(y-10<llimit)
                    return
                end
                y = y  - 10;
               stage.SetY(y);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
        end
    end
    function zup_call(fileTab,eventdata)
        str = get(handles.z_val, 'String');
        val = get(handles.z_val,'Value');
        switch str{val};
            case '.1 um'
                if(z+0.1>ulimit)
                    return
                end
                z = z + 0.1;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.5 um'
                if(z+0.5>ulimit)
                    return
                end
                z = z + 0.5;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '1.0 um'
                if(z+1>ulimit)
                    return
                end
                z = z + 1.0;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '2.0 um'
                if(z+2>ulimit)
                    return
                end
                z = z + 2;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '5.0 um'
                if(z+5>ulimit)
                    return
                end
                z = z + 5;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '10 um'
                if(z+10>ulimit)
                    return
                end
                z = z + 10;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
        end
    end
    function zdown_call(fileTab,eventdata)
        str = get(handles.z_val, 'String');
        val = get(handles.z_val,'Value');
        switch str{val};
            case '.1 um'
                if(z-0.1<llimit)
                    return
                end
                z = z - 0.1;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '0.5 um'
                if(z-0.5<llimit)
                    return
                end
                z = z - 0.5;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '1.0 um'
                if(z-1<llimit)
                    return
                end
                z = z - 1.0;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '2.0 um'
                if(z-2<llimit)
                    return
                end
                z = z - 2;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '5.0 um'
                if(z-5<llimit)
                    return
                end
                z = z - 5;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
            case '10 um'
                if(z-10<llimit)
                    return
                end
                z = z - 10;
                stage.SetZ(z);
                pause(0.1);
                [a,b,c]=stage.GetPosition;
        end
    end
%% Refresh relative location of nano-stage
    function refresh_call(fileTab,eventdata)
        [x,y,z] = stage.GetPosition;
        set(handles.xpos_txt, 'String', ['X:     ' num2str(x,'%.3f') ' um']);
        set(handles.ypos_txt, 'String', ['Y:     ' num2str(y,'%.3f') ' um']);
        set(handles.zpos_txt, 'String', ['Z:     ' num2str(z,'%.3f') ' um']);
    end
%% Manual Mode/Phase Contrast
    function phasec_call(fileTab,eventdata) % issue with flip mirror being up or down
        [DI,t1]=inputSingleScan(g);
%         if(DI(1) == 0)
%             outputSingleScan(s1,0);
%             tic;
%             while toc < 0.01
%             end
%             [DI,t1]=inputSingleScan(g);
%             if(DI(1) == 0)
%                 outputSingleScan(s1,1);
%                 tic;
%                 while toc < 0.01
%                 end
%             end
%             set(handles.shut1_chck, 'String', 'Open 405nm Laser');
%         elseif(DI(1)==1)
%         end
        if(DI(2) == 0)
            outputSingleScan(s2,0);
            tic;
            while toc < 0.01
            end
            [DI,t2]=inputSingleScan(g);
            if(DI(2) == 0)
                outputSingleScan(s2,1);
                tic;
                while toc < 0.01
                end
            end
            set(handles.shut2_chck, 'String', 'Open 515nm Laser');
        elseif(DI(2)==1)
        end
        if(DI(3) == 0)
            outputSingleScan(s3,0);
            tic;
            while toc < 0.01
            end
            [DI,t3]=inputSingleScan(g);
            if(DI(3) == 0)
                outputSingleScan(s3,1);
                tic;
                while toc < 0.01
                end
            end
            set(handles.shut3_chck, 'String', 'Open 561nm Laser');
        elseif(DI(3)==1)
        end
        set(handles.shut1_chck, 'enable', 'off');
        set(handles.shut2_chck, 'enable', 'off');
        set(handles.shut3_chck, 'enable', 'off');
        set(handles.manual_but,'Enable','on');
        set(handles.manual_but,'value',0);
        set(handles.phasec_but,'Enable','off');
    end
    function manual_call(fileTab,eventdata)
        set(handles.phasec_but,'enable','on');
        set(handles.phasec_but,'value',0);
        set(handles.manual_but,'enable','off');
        set(handles.shut1_chck, 'enable', 'on');
        set(handles.shut2_chck, 'enable', 'on');
        set(handles.shut3_chck, 'enable', 'on');
        set(handles.fmirror_chck, 'enable', 'on');
        set(handles.LED_chck, 'enable', 'on');
    end
    function imgsrc_call(fileTab,eventdata)
        imagingSourceCamera;
    end
%% Triggering devices
%     function shut1_call(fileTab,eventdata)
%         [DI,t1]=inputSingleScan(g);
%         if(DI(1) == 0)                        % input of 0 means shutter is open
%             outputSingleScan(s1,0);
%             pause(0.1);
%             [DI,t1]=inputSingleScan(g);
%             if(DI(1) == 0)
%                 outputSingleScan(s1,1);
%                 pause(0.1);
%             end
%             set(handles.shut1_chck, 'String', 'Open 405nm Laser');
%         elseif(DI(1)==1)                      % input of 1 means shutter is closed
%             outputSingleScan(s1,1);
%             pause(0.1);
%             [DI,t1]=inputSingleScan(g);
%             if(DI(1) == 1)
%                 outputSingleScan(s1,0);
%                 pause(0.1);
%             end
%             set(handles.shut1_chck, 'String', 'Close 405nm Laser');
%         end
%     end
    
    function Ctrl_OBIS_laser(fileTab,eventdata)
        if (get(handles.Ctrl_OBIS_laser,'Value') == get(handles.Ctrl_OBIS_laser,'Max'))
            outputSingleScan(o1,1);
        else
            outputSingleScan(o1,0);
        end
        if(strcmp(get(handles.Ctrl_OBIS_laser, 'String'), 'Open 405nm Laser') == 0)
            set(handles.Ctrl_OBIS_laser, 'String', 'Open 405nm Laser');
        else
            set(handles.Ctrl_OBIS_laser, 'String', 'Close 405nm Laser');
        end
    end
    
%     function Set_OBIS_laser_0(Exposure_Time,Waiting_Time)
%         Set_OBIS_laser;
%         function Set_OBIS_laser
%             prompt = {'Exposure Time (ms):','Waiting after turn off (ms):'};
%             dlg_title = 'Initialize OBIS Laser';
%             defaultimes = {'50','10'};
%             answer = inputdlg(prompt,dlg_title,[1 50],defaultimes);
%             OBIS_exp = str2double(answer{1});
%             OBIS_waiting_off = str2double(answer{2});
%         end
%         outputSingleScan(o1,1);
%         tic;
%         while toc < (0.01)
%         end
%         outputSingleScan(o1,0);
%         % Waiting time after turn off the OBIS laser
%         tic;
%         while toc > (0.01)
%         end
%     end

    function shut2_call(fileTab,eventdata)
        [DI,t2]=inputSingleScan(g);
        if(DI(2) == 0)
            outputSingleScan(s2,0);
            pause(0.1);
            [DI,t2]=inputSingleScan(g);
            if(DI(2) == 0)
                outputSingleScan(s2,1);
                pause(0.1);
            end
            set(handles.shut2_chck, 'String', 'Open 515nm Laser');
        elseif(DI(2)==1)
            outputSingleScan(s2,1);
            pause(0.1);
            [DI,t2]=inputSingleScan(g);
            if(DI(2) == 1)
                outputSingleScan(s2,0);
                pause(0.1);
            end
            set(handles.shut2_chck, 'String', 'Close 515nm Laser');
        end
    end
    function shut3_call(fileTab,eventdata)
        [DI,t3]=inputSingleScan(g);
        if(DI(3) == 0)
            outputSingleScan(s3,0);
            pause(0.1);
            [DI,t3]=inputSingleScan(g);
            if(DI(3) == 0)
                outputSingleScan(s3,1);
                pause(0.1);
            end
            set(handles.shut3_chck, 'String', 'Open 561nm Laser');
        elseif(DI(3)==1)
            outputSingleScan(s3,1);
            pause(0.1);
            [DI,t3]=inputSingleScan(g);
            if(DI(3) == 1)
                outputSingleScan(s3,0);
                pause(0.1);
            end
            set(handles.shut3_chck, 'String', 'Close 561nm Laser');
        end
    end
    function fmir_call(fileTab,eventdata)
        if (get(handles.fmirror_chck,'Value') == get(handles.fmirror_chck,'Max'))
            outputSingleScan(f1,1);
        else
            outputSingleScan(f1,0);
        end
        if(strcmp(get(handles.fmirror_chck, 'String'), 'Toggle Mirror Up') == 0)
            set(handles.fmirror_chck, 'String', 'Toggle Mirror Up');
        else
            set(handles.fmirror_chck, 'String', 'Toggle Mirror Down');
        end
    end
    function led_call(fileTab,eventdata)
        if (get(handles.LED_chck,'Value') == get(handles.LED_chck,'Max'))
            outputSingleScan(l1,1);
        else
            outputSingleScan(l1,0);
        end
    end
%% Load excel file with instructions *see Z:\Eugene\input.xlsx
    function load_call(fileTab,eventdata)
        % excel file
        stop(tmr); % pause timer
        [FileName,PathName,FilterIndex] = uigetfile('*.xls;*.xlsx','Select Sequence Instructions');
        cd(PathName);
        e = xlsread(FileName);
        set(handles.file_txt, 'String', ['FILE LOADED: ' FileName]);
        cd('Z:\Eugene\development R161_FromLidkeLab');
        start(tmr);
    end
%% Run the sequence
    function start_call(fileTab,eventdata)
        FolderName = uigetdir('Z:\','Select Output Save Directory'); % select where to save log files
        
        choice = questdlg('Do you want to use stage instructions from input?', ... % Construct a questdlg with three options
            'Stage Instructions', ...
            'Yes','No, use current position', 'Yes');
        
        switch choice  % Handle response
            case 'Yes'
                ignore = 0;
            case 'No, use current position'
                ignore = 1;
        end
        
        [ff gg hh] = stage.GetPosition;
        s.Rate = 1000;
        if(strcmp(FileName,'') == 1) % if no file loaded, return
            return
        end
        
        stop(tmr); % stop timer until after sequence is completed
        
        [nn m] = size(e); % length of instructions
        loop = e(1,1); % number of loops of instructions
        cam1_frms = 0; % total frames from camera 1 per loop
        for q = 5:nn
            cam1_frms = e(q, 8) + cam1_frms;
        end
        cam2_frms = 0; % total frames from camera 2 per loop
        for p = 5:nn
            cam2_frms = e(p, 10) + cam2_frms;
        end
        
        totlines = cam1_frms + nn - 4; % total lines for output file [instructions*(frames+1)]
        zout = zeros(totlines*loop,1,'double'); % output array for z-pos
        shtr1 = ones(totlines*loop,1,'single'); % output array for shutter 1
        shtr2 = ones(totlines*loop,1,'single'); % output array for shutter 2
        shtr3 = ones(totlines*loop,1,'single'); % output array for shutter 3
        totlines2 = cam2_frms + nn - 4;
        zout2 = zeros(totlines2*loop,1,'double'); % same as above for camera 2
        shtr1_2 = ones(totlines2*loop,1,'single');
        shtr2_2 = ones(totlines2*loop,1,'single');
        shtr3_2 = ones(totlines2*loop,1,'single');
        [x,y,z] = stage.GetPosition;
        
        if(ignore == 0)
            stage.SetPosition(ff,gg,e(5,6)); % set starting position
            ON = tic;
            while toc(ON) < (0.1) % pause
            end
        else
        end
        
        set(handles.sequence_txt, 'String', ['SEQUENCE STATUS: ' 'Sequence Running'],'fontweight','bold');
        pause(0.0001);
        
        output_ind = 1; % create index variable to fill output array
        output_ind_2 = 1; % create index variable to fill camera 2 output array
        for j = 1:loop
            currentLoop = j
            
            [DI,t3]=inputSingleScan(g); % initialize shutters (turn off all shutters)
            if((DI(2) == 1) && (DI(3) == 1))
%             elseif DI(1) == 0
%                 outputSingleScan(s1,0);
            elseif DI(2) == 0
                outputSingleScan(s2,0);
            elseif DI(3) == 0
                outputSingleScan(s3,0);
            end
            
            for ii = 5:nn % loop through the instructions
                if(abort~=0)
                    abort = 0;
                    %% output log file if aborted
                    [zoutsize bla] = size(zout);
                    zoutind=1;
                    zoutnew=zeros(1,1,'single');
                    shtr1new=zeros(1,1,'single');
                    shtr2new=zeros(1,1,'single');
                    shtr3new=zeros(1,1,'single');
                    while(zout(zoutind,1) ~= 0 && zoutind ~= zoutsize)
                        zoutind = zoutind + 1;
                    end
                    if(zoutind == zoutsize) % create output file for cam1 after aborted sequence
                    else
                        for u=1:zoutind-1
                            zoutnew(u,1) = zout(u,1);
                            shtr1new(u,1) = shtr1(u,1);
                            shtr2new(u,1) = shtr2(u,1);
                            shtr3new(u,1) = shtr3(u,1);
                        end
                    end
                    [zoutsize2 bla] = size(zout2);
                    zoutind2=1;
                    zout2new=zeros(1,1,'single');
                    shtr1_2new=zeros(1,1,'single');
                    shtr2_2new=zeros(1,1,'single');
                    shtr3_2new=zeros(1,1,'single');
                    while(zout2(zoutind2,1) ~= 0 && zoutind2 ~= zoutsize2)
                        zoutind2 = zoutind2 + 1;
                    end
                    if(zoutind2 == zoutsize2) % create output file for cam2 after aborted sequence
                    else
                        for u=1:zoutind2-1
                            zout2new(u,1) = zout2(u,1);
                            shtr1_2new(u,1) = shtr1_2(u,1);
                            shtr2_2new(u,1) = shtr2_2(u,1);
                            shtr3_2new(u,1) = shtr3_2(u,1);
                        end
                    end
                    output = horzcat(shtr1new, shtr2new, shtr3new, zoutnew); % concatonate shutter and z position data
                    output2 = horzcat(shtr1_2new,shtr2_2new, shtr3_2new,zout2new);
                    output_filename = 'output_cam1'; % NAME THE OUTPUT FILE
                    output_filename2 = 'output_cam2';
                    
                    cd(FolderName); %% output actual locations of nanostage to .txt file if sequence if aborted
                    c= clock;
                    fix(c);
                    if(exist([output_filename '.dat'], 'file')==2) % if file exists, create new file
                        output_filename = [output_filename '_' int2str(c(1)) int2str(c(2)) int2str(c(3)) '_' int2str(c(4)) int2str(c(5)) int2str(c(6))];
                    end
                    fid = fopen([output_filename '.dat'], 'a');
                    fprintf(fid, '%d %d %d %f\r\n', output');
                    fclose(fid);
                    if(exist([output_filename2 '.dat'], 'file')==2) % if file exists, create new file
                        output_filename2 = [output_filename2 '_' int2str(c(1)) int2str(c(2)) int2str(c(3)) '_' int2str(c(4)) int2str(c(5)) int2str(c(6))];
                    end
                    fid = fopen([output_filename2 '.dat'], 'a');
                    fprintf(fid, '%d %d %d %f\r\n', output2');
                    fclose(fid);
                    cd('Z:\Eugene\development R161_FromLidkeLab');
                    ignore = 0;
                    
                    set(handles.sequence_txt, 'String', ['SEQUENCE STATUS: ' 'Sequence Not Running'],'fontweight','normal');
                    start(tmr);
                    return; % sequence aborted
                end
                
                exp_time = e(ii,7);
                exp_time_2 = e(ii,9);
                
                if(ignore == 0)
                    stage.SetPosition(ff,gg,e(ii,6)); % increment to new step
                    ON = tic;
                    while toc(ON) < (0.1)
                    end
                else
                end
                
                shtr1(output_ind,1) = -1; % indicate transition in output file (cam1)
                shtr2(output_ind,1) = -1;
                shtr3(output_ind,1) = -1;
                zout(output_ind, 1) = z;
                shtr1_2(output_ind_2,1) = -1; % indicate transition in output file (cam2)
                shtr2_2(output_ind_2,1) = -1;
                shtr3_2(output_ind_2,1) = -1;
                zout2(output_ind_2, 1) = z;
                
                ins = e(ii,8) * (delayT+exp_time); % instruction size for cam 1
                ins2 = e(ii,10) * (delayT+exp_time_2); % instruction size for cam 2 [# of frames * (delay between triggers (~11ms) + exposure (10-50ms))]
                arr = zeros(ins,1,'double');
                arr2 = zeros(ins2,1,'double');
                
                if(ins == 0) % pause
                else
                    ind = 1;
                    while(ind < ins)
                        for i=1:delayT
                            if(ind < ins)
                                arr(ind,1)=0;
                                ind=ind +1;
                            end
                        end
                        output_ind = output_ind + 1;
                        [x,y,z] = stage.GetPosition;
                        shtr1(output_ind,1) = 1;
                        shtr2(output_ind,1) = 1;
                        shtr3(output_ind,1) = 1;
                        zout(output_ind,1) = z; % take the z position and store for output
                        for j=i:(i+(exp_time-1))
                            if(ind < ins)
                                arr(ind,1)=5;
                                ind=ind+1;
                            end
                        end
                    end
                end
                if(ins2 == 0)
                else
                    ind2 = 1;
                    while(ind2 < ins2)
                        for in=1:delayT
                            if(ind2 < ins2)
                                arr2(ind2,1)=0;
                                ind2=ind2 +1;
                            end
                        end
                        output_ind_2 = output_ind_2 + 1;
                        [x,y,z] = stage.GetPosition;
                        shtr1_2(output_ind_2,1) = 1;
                        shtr2_2(output_ind_2,1) = 1;
                        shtr3_2(output_ind_2,1) = 1;
                        zout2(output_ind_2,1) = z; % take the z position and store for output
                        for j=in:(in+(exp_time_2-1))
                            if(ind2 < ins2)
                                arr2(ind2,1)=5;
                                ind2=ind2+1;
                            end
                        end
                    end
                end
                output_ind = output_ind + 1;
                output_ind_2 = output_ind_2 + 1;
                [n m] = size(arr2);
                [t y] = size(arr);
                while(n>t) % depending on which camera collects more frames, fill the instructions of the other camera with 0 to match the size
                    [t y] = size(arr);
                    arr(t+1,1) = 0;
                    [n m] = size(arr2);
                    [t y] = size(arr);
                end
                while(n<t)
                    [n m] = size(arr2);
                    arr2(n+1,1) = 0;
                    [t y] = size(arr);
                    [n m] = size(arr2);
                end
                
                %% adjust shutter according to input
                % emit 405
                if e(ii,1) == 1
                    outputSingleScan(o1,1);
                    % On time of the OBIS laser
                    tic;
                    while toc < (e(ii,11)/1000)
                    end
                    outputSingleScan(o1,0);
                    % Waiting time after turn off the OBIS laser
                    if e(ii,12) == 0
                    else
                        tic;
                        while toc < (e(ii,12)/1000)
                        end
                    end
                end
                % emit 514
                if e(ii,2)==1
                    outputSingleScan(s2,1);
                    tic;
                    while toc < 0.01
                    end
                    %                     [D2,t3]=inputSingleScan(g);
                    %                     tic;
                    %                     while toc < 0.01
                    %                     end
                    
                    %                     if(D2(2) == 1)
                    %                         outputSingleScan(s2,0);
                    %                         tic;
                    %                         while toc < 0.01
                    %                         end
                    %                     end
                end
                
                % emit 561
                if e(ii,3)== 1
                    outputSingleScan(s3,1);
                    tic;
                    while toc < 0.01
                    end
                    %                     [D3,t3]=inputSingleScan(g);
                    %                     tic;
                    %                     while toc < 0.01
                    %                     end
                    %                     if(D3(3) == 0)
                    %                         outputSingleScan(s3,1);
                    %                         tic;
                    %                         while toc < 0.01
                    %                         end
                    %                     end
                end
                
                % exposure camera
                if (e(ii,8)==0 && e(ii,10)==0)
                else
                    queueOutputData(s,[arr arr2]);
                    ON = tic;
                    while toc(ON) < (0.11)
                    end
                    s.startForeground;
                    ON = tic;
                    while toc(ON) < (0.11)
                    end
                end
                
                % after exposure cameras, turn off corresponding lasers
                [DI,t3]=inputSingleScan(g);
                if((DI(2) == 1) && (DI(3) == 1))
%                 elseif DI(1) == 0
%                     outputSingleScan(s1,0);
                elseif DI(2) == 0
                    outputSingleScan(s2,0); % turn off 514
                    tic;
                    while toc < 0.01
                    end
                elseif DI(3) == 0
                    outputSingleScan(s3,0); % turn off 561
                    tic;
                    while toc < 0.01
                    end
                end                
            end
        end
        
        [DI,t1]=inputSingleScan(g);
%         if(DI(1) == 0)
%             outputSingleScan(s1,0);
%             tic;
%             while toc < 0.01
%             end
%             [DI,t1]=inputSingleScan(g);
%             if(DI(1) == 0)
%                 outputSingleScan(s1,1);
%                 tic;
%                 while toc < 0.01
%                 end
%             end
%             set(handles.shut1_chck, 'String', 'Open 405nm Laser');
%         elseif(DI(1)==1)
%         end
        if(DI(2) == 0)
            outputSingleScan(s2,0);
            tic;
            while toc < 0.01
            end
            [DI,t2]=inputSingleScan(g);
            if(DI(2) == 0)
                outputSingleScan(s2,1);
                tic;
                while toc < 0.01
                end
            end
            set(handles.shut2_chck, 'String', 'Open 515nm Laser');
        elseif(DI(2)==1)
        end
        if(DI(3) == 0)
            outputSingleScan(s3,0);
            tic;
            while toc < 0.01
            end
            [DI,t3]=inputSingleScan(g);
            if(DI(3) == 0)
                outputSingleScan(s3,1);
                tic;
                while toc < 0.01
                end
            end
            set(handles.shut3_chck, 'String', 'Open 561nm Laser');
        elseif(DI(3)==1)
        end
        
        output = horzcat(shtr1, shtr2, shtr3, zout); % concatonate shutter and z position data
        output2 = horzcat(shtr1_2,shtr2_2, shtr3_2,zout2);
        output_filename = 'output_cam1'; % NAME THE OUTPUT FILE
        output_filename2 = 'output_cam2';
        
        %% output actual locations of nanostage to .txt file
        cd(FolderName);
        c= clock;
        fix(c);
        if(exist([output_filename '.dat'], 'file')==2) % if file exists, create new file
            output_filename = [output_filename '_' int2str(c(1)) int2str(c(2)) int2str(c(3)) '_' int2str(c(4)) int2str(c(5)) int2str(c(6))];
        end
        fid = fopen([output_filename '.dat'], 'a');
        fprintf(fid, '%d %d %d %f\r\n', output');
        fclose(fid);
        if(exist([output_filename2 '.dat'], 'file')==2) % if file exists, create new file
            output_filename2 = [output_filename2 '_' int2str(c(1)) int2str(c(2)) int2str(c(3)) '_' int2str(c(4)) int2str(c(5)) int2str(c(6))];
        end
        fid = fopen([output_filename2 '.dat'], 'a');
        fprintf(fid, '%d %d %d %f\r\n', output2');
        fclose(fid);
        cd('Z:\Eugene\development R161_FromLidkeLab');
        ignore = 0;
        set(handles.sequence_txt, 'String', ['SEQUENCE STATUS: ' 'Sequence Not Running'],'fontweight','normal');
        start(tmr);
    end

%% Abort GUI
    function exit_call(fileTab,eventdata)
        [DI,t1]=inputSingleScan(g);
%         if(DI(1) == 0)
%             outputSingleScan(s1,0);
%             tic;
%             while toc < 0.01
%             end
%             [DI,t1]=inputSingleScan(g);
%             if(DI(1) == 0)
%                 outputSingleScan(s1,1);
%                 tic;
%                 while toc < 0.01
%                 end
%             end
%             set(handles.shut1_chck, 'String', 'Open 405nm Laser'); % close all shutters
%         elseif(DI(1)==1)
%         end
        if(DI(2) == 0)
            outputSingleScan(s2,0);
            tic;
            while toc < 0.01
            end
            [DI,t2]=inputSingleScan(g);
            if(DI(2) == 0)
                outputSingleScan(s2,1);
                tic;
                while toc < 0.01
                end
            end
            set(handles.shut2_chck, 'String', 'Open 515nm Laser');
        elseif(DI(2)==1)
        end
        if(DI(3) == 0)
            outputSingleScan(s3,0);
            tic;
            while toc < 0.01
            end
            [DI,t3]=inputSingleScan(g);
            if(DI(3) == 0)
                outputSingleScan(s3,1);
                tic;
                while toc < 0.01
                end
            end
            set(handles.shut3_chck, 'String', 'Open 561nm Laser');
        elseif(DI(3)==1)
        end
        
        outputSingleScan(l1,0); % turn off LED
        
        MCLNanoDriveStage.ReleaseAllHandles; % remove handles from nanodrive
        MCLMicroDriveStage.ReleaseAllHandles;
        close(handles.output);
        clear;
    end
    function abort_call(fileTab,eventdata)
        abort = 1;
        s.stop;
        
        [DI,t1]=inputSingleScan(g);
%         if(DI(1) == 0)
%             outputSingleScan(s1,0);
%             tic;
%             while toc < 0.01
%             end
%             [DI,t1]=inputSingleScan(g);
%             if(DI(1) == 0)
%                 outputSingleScan(s1,1);
%                 tic;
%                 while toc < 0.01
%                 end
%             end
%             set(handles.shut1_chck, 'String', 'Open 405nm Laser'); % close shutters
%         elseif(DI(1)==1)
%         end
        if(DI(2) == 0)
            outputSingleScan(s2,0);
            tic;
            while toc < 0.01
            end
            [DI,t2]=inputSingleScan(g);
            if(DI(2) == 0)
                outputSingleScan(s2,1);
                tic;
                while toc < 0.01
                end
            end
            set(handles.shut2_chck, 'String', 'Open 515nm Laser');
        elseif(DI(2)==1)
        end
        if(DI(3) == 0)
            outputSingleScan(s3,0);
            tic;
            while toc < 0.01
            end
            [DI,t3]=inputSingleScan(g);
            if(DI(3) == 0)
                outputSingleScan(s3,1);
                tic;
                while toc < 0.01
                end
            end
            set(handles.shut3_chck, 'String', 'Open 561nm Laser');
        elseif(DI(3)==1)
        end
    end

%% Z-calibration (retrieve information)
    function zstep_call(fileTab,eventdata)
        step = str2double(get(handles.zcal_step,'String'));
    end
    function zsize_call(fileTab,eventdata)
        step_size = str2double(get(handles.zcal_size,'String'));
    end
    function zframe_call(fileTab,eventdata)
        frame = str2double(get(handles.zcal_frame,'String'));
    end
    function ztime_call(fileTab,eventdata)
        exposure = str2double(get(handles.zcal_exptime,'String'));
    end

%% Z-calibration run
    function zcal_call(fileTab,eventdata)
        FolderName = uigetdir('Z:\','Select Output Save Directory'); % select save directory for log file
        s.Rate = 1000;
        stop(tmr);
        [ff,gg,hh] = stage.GetPosition;
        step_size = step_size/1000;
        if((hh-step_size*step)<0 || (hh+step_size*step)>200) % ensure z-calibration stays within limits
            display('Z-calibration goes out of bounds');
            return;
        end
        if(step==0 || step_size==0 || frame ==0|| exposure == 0)
            return
        end
        exposure = exposure * 1000; % convert into miliseconds
        cc1 = zeros(step,1,'double');
        cc2 = zeros(step*2,1,'double');
        cc3 = zeros(step,1,'double');
        for m = 1:step % create instructions for z-calibration (start from origin and increment down)
            cc1(m,1) = hh - step_size*m;
        end
        for n = 1:step*2 % from the lowest z-position, increment up past the origin to the highes z-position
            cc2(n,1) = (hh-step*step_size) + step_size*n;
        end
        for p = 1:step % return to the origin
            cc3(p,1) = (hh+step*step_size) - step_size*p;
        end
        ccfinal = vertcat(cc1, cc2, cc3);
        %create instructions for shutters
        shtr1 = ones(step*4*(frame+1),1,'single');
        shtr2 = ones(step*4*(frame+1),1,'single');
        shtr3 = ones(step*4*(frame+1),1,'single');
        % execute instructions
        zout=zeros(4*step*(frame+1),1,'double'); % output for z positions
        yout=zeros(4*step*(frame+1),1,'double'); % output for z positions
        xout=zeros(4*step*(frame+1),1,'double'); % output for z positions
        stage.SetPosition(ff,gg,hh); % set starting position
        ON = tic;
        while toc(ON) < (0.1)
        end
        set(handles.sequence_txt, 'String', ['SEQUENCE STATUS: ' 'Sequence Running'],'fontweight','bold');
        pause(0.0001);
        output_ind = 1;
        for i = 1:(size(ccfinal))
            if(abort~=0) % abort z-calibration
                abort = 0;
                set(handles.sequence_txt, 'String', ['SEQUENCE STATUS: ' 'Sequence Not Running'],'fontweight','normal');
                step_size = str2double(get(handles.zcal_size,'String'));
                exposure = str2double(get(handles.zcal_exptime,'String'));
                start(tmr);
                return;
            end
            [x,y,z] = stage.GetPosition;
            shtr1(output_ind,1) = -1;
            shtr2(output_ind,1) = -1;
            shtr3(output_ind,1) = -1;
            zout(output_ind, 1) = z;
            xout(output_ind,1) = x;
            yout(output_ind,1) = y;
            
            stage.SetZ(ccfinal(i,1)); % increment to new step
            ON = tic;
            while toc(ON) < (0.1)
            end
            
            ins = (frame * (delayT+exposure)); % size of instructions array
            a = zeros(ins,1,'double'); % create instructions
            ind = 1;
            while(ind < ins)
                for d=1:delayT
                    if(ind < ins)
                        a(ind,1)=0;
                        ind=ind +1;
                    end
                end
                output_ind = output_ind + 1;
                [x,y,z] = stage.GetPosition;
                shtr1(output_ind,1) = 1;
                shtr2(output_ind,1) = 1;
                shtr3(output_ind,1) = 1;
                zout(output_ind,1) = z; % take the z position and store for output
                xout(output_ind,1) = x; 
                yout(output_ind,1) = y; 
                for j=d:(d+(exposure-1))
                    if(ind < ins)
                        a(ind,1)=5;
                        ind=ind+1;
                    end
                end
            end
            
            [pp qq] = size(a); % fill in the last cells with instructions
            a(pp,1) = 5;
            for o=1:delayT
                a(pp+o,1) = 0;
            end
            
            output_ind = output_ind + 1;
            queueOutputData(s,[a a]);
            ON = tic;
            while toc(ON) < (0.11)
            end
            s.startForeground; % begin acquisition for that step
            ON = tic;
            while toc(ON) < (0.11)
            end
        end
        output = horzcat(shtr1, shtr2, shtr3, zout, xout, yout); % concatonate shutter and z position data
        output_filename = 'output_zcal'; % NAME THE OUTPUT FILE
        
        %% output actual locations of nanostage to .txt file
        cd(FolderName);
        c= clock;
        fix(c);
        if(exist([output_filename '.dat'], 'file')==2) % if file exists, create new file
            output_filename = [output_filename '_' int2str(c(1)) int2str(c(2)) int2str(c(3)) '_' int2str(c(4)) int2str(c(5)) int2str(c(6))];
        end
        fid = fopen([output_filename '.dat'], 'a');
        fprintf(fid, '%d %d %d %d %d %f\r\n', output');
        fclose(fid);
        cd('Z:\Eugene\development R161_FromLidkeLab');
        set(handles.sequence_txt, 'String', ['SEQUENCE STATUS: ' 'Sequence Not Running'],'fontweight','normal');
        step_size = str2double(get(handles.zcal_size,'String'));
        exposure = str2double(get(handles.zcal_exptime,'String'));
        start(tmr); % resume timer (update positions)
    end

%% Refresh position (manually retrieve position of micro-stage)
    function refresh1_call(fileTab,eventdata)
        [xM,yM] = stage2.position;
        set(handles.xpos_txt1, 'String', ['X:     ' num2str(xM,'%.3f') ' um']);
        set(handles.ypos_txt1, 'String', ['Y:     ' num2str(yM,'%.3f') ' um']);
    end
%% Grid sequence analysis
% initialize parameters
    function xyStep_call(fileTab,eventdata)
        xyStep = str2double(get(handles.xyStep_val,'String'));
    end
    function zStep_call(fileTab,eventdata)
        zStep = str2double(get(handles.zStep_val,'String'));
    end
    function xyStepSize_call(fileTab,eventdata)
        xyStepSize = str2double(get(handles.xyStepSize_val,'String'));
    end
    function zStepSize_call(fileTab,eventdata)
        zStepSize = str2double(get(handles.zStepSize_val,'String'));
    end
    function gridFrame_call(fileTab,eventdata)
        gridFrame = str2double(get(handles.gridFrame_val,'String'));
    end
    function gridExposure_call(fileTab,eventdata)
        gridExposure = str2double(get(handles.gridExposure_val,'String'));
    end
    function grid_call(fileTab,eventdata)
        FolderName = uigetdir('Z:\','Select Output Save Directory'); % select directory to save log file
        set(handles.sequence_txt, 'String', ['SEQUENCE STATUS: ' 'Sequence Running'],'fontweight','bold');
        pause(0.0001);
        
        s.Rate = 1000;
        stop(tmr);
        
        [xor,yor,zor] = stage.GetPosition;
        [x, y, z] = stage.GetPosition;
        %do not change these variables
        delayGrid = gridExposure + 1;
        distance_x = xor + xyStep*xyStepSize;
        distance_y = yor + xyStep*xyStepSize;
        xi = xor;
        yi = yor;
        
        %z- instructions
        zamp = zStepSize * zStep;
        zstart = zor - zamp;
        zend = zor + zamp;
        zinstr = zeros(((zStep*2)+1),1,'double');
        zind = 1;
        for k = zstart:0.3:zend
            zinstr(zind,1) = k;
            zind = zind + 1;
        end
        
        stage.SetPosition(xor, yor, zstart);
        tic;
        while toc < (0.1)
        end
        
        forwardarr = zeros(2,9,'double'); % instructions for snaking through grid
        backwardarr = zeros(2,9,'double'); % instructions for returning to the starting position, going backwards
        ind=1;
        forwardarr(1,ind) =  xi;                      % origin
        forwardarr(2,ind) =  yi;
        
        while (xi ~= distance_x || yi ~=distance_y)        %instructions for snake loop
            if(xi == xor)
                while(xi~=distance_x)
                    xi = xi+xyStepSize;
                    ind = ind+1;
                    forwardarr(1,ind) =  xi;
                    forwardarr(2,ind) =  yi;
                end
                if(xi==distance_x && yi==distance_y)
                    break
                end
                if(xi == distance_x && yi ~=distance_y)
                    yi = yi+xyStepSize;
                    ind = ind+1;
                    forwardarr(1,ind) =  xi;
                    forwardarr(2,ind) =  yi;
                end
                if(xi==distance_x &&yi==distance_y)
                    break
                end
            end
            
            if(xi == distance_x)
                while( xi~=xor)
                    xi = xi-xyStepSize;
                    ind = ind+1;
                    forwardarr(1,ind) =  xi;
                    forwardarr(2,ind) =  yi;
                end
                if(xi==distance_x &&yi==distance_y)
                    break
                end
                
                if(xi == xor && yi ~=distance_y)
                    yi = yi+xyStepSize;
                    ind = ind+1;
                    forwardarr(1,ind) =  xi;
                    forwardarr(2,ind) =   yi;
                end
                if(xi==distance_x && yi==distance_y)
                    break
                end
            end
        end
        
        [zz, ind2g] = size(forwardarr);
        while(mod(sqrt(ind2g),1)~=0)                           % finish the grid in case of 'even' NxM
            forwardarr(1,ind2g+1)=forwardarr(1,ind2g) - xyStepSize;
            forwardarr(2,ind2g+1)=forwardarr(2,ind2g);
            [zz, ind2g] = size(forwardarr);
        end
        
        for i=1:ind2g                                          % copy reverse instructions for backwards
            backwardarr(1,i) = forwardarr(1,ind2g-i+1);
            backwardarr(2,i) = forwardarr(2,ind2g-i+1);
        end
        %         forwardarr
        
        %% Run through the grid
        zcal = 1;
        [zcal_size, o] = size(zinstr);
        
        result_ind=1;
        
        b = 0; % meaning index is at origin
        % b = 1 meaning index is at (xend, yend)
        results_grid = zeros(3,1,'double');
        ins = (gridFrame * (delayGrid+gridExposure));
        a = zeros(ins,1,'double');
        ind = 1;
        while(ind < ins)
            for d=1:delayGrid
                if(ind < ins)
                    a(ind,1)=0;
                    ind=ind +1;
                end
            end
            for j=d:(d+(gridExposure-1))
                if(ind < ins)
                    a(ind,1)=5;
                    ind=ind+1;
                end
            end
        end
        [pp, qq] = size(a);
        a(pp,1) = 5;
        for o=1:delayGrid
            a(pp+o,1) = 0;
        end
        %
        %         xyStep = 7;
        %         zStep = 5;
        %         xyStepSize = 0.5;
        %         zStepSize = 0.3;
        %         gridFrame = 10;
        %         gridExposure = 50;
        
        
        %create instructions for shutters
        [DI,t3]=inputSingleScan(g);
        if(DI(1) == 0)
            shtr1 = ones(((xyStep+1)*(xyStep+1)*((zStep*2)+1)*(gridFrame+1)),1,'single');
        else
            shtr1 = zeros(((xyStep+1)*(xyStep+1)*((zStep*2)+1)*(gridFrame+1)),1,'single');
        end
        if(DI(2) == 0)
            shtr2 = ones(((xyStep+1)*(xyStep+1)*((zStep*2)+1)*(gridFrame+1)),1,'single');
        else
            shtr2 = zeros(((xyStep+1)*(xyStep+1)*((zStep*2)+1)*(gridFrame+1)),1,'single');
        end
        if(DI(3) == 0)
            shtr3 = ones(((xyStep+1)*(xyStep+1)*((zStep*2)+1)*(gridFrame+1)),1,'single');
        else
            shtr3 = zeros(((xyStep+1)*(xyStep+1)*((zStep*2)+1)*(gridFrame+1)),1,'single');
        end
        
        % execute instructions
        zout=zeros(((xyStep+1)*(xyStep+1)*((zStep*2)+1)*(gridFrame+1)),1,'double'); % output for z positions
        xout=zeros(((xyStep+1)*(xyStep+1)*((zStep*2)+1)*(gridFrame+1)),1,'double'); % output for z positions
        yout=zeros(((xyStep+1)*(xyStep+1)*((zStep*2)+1)*(gridFrame+1)),1,'double'); % output for z positions
        
        output_ind = 1;
        while(zcal ~= zcal_size + 1)
            wait(s);
            if(abort~=0)
                abort = 0;
                set(handles.sequence_txt, 'String', ['SEQUENCE STATUS: ' 'Sequence Not Running'],'fontweight','normal');
                start(tmr);
                return;
            end
            [fgrid,ggrid,h] = stage.GetPosition;
            stage.SetZ(zinstr(zcal,1));
            tic;
            while toc < (0.1)
            end
            
            if(b == 0)
                [o, forwardsz] = size(forwardarr);
                for f = 1:forwardsz
                    setX = forwardarr(1,f);
                    setY = forwardarr(2,f);
                    %            stage.xset(setX);
                    %            stage.yset(setY);
                    [i o p] = stage.GetPosition;
                    shtr1(output_ind,1) = -1; % indicate transition in output file
                    shtr2(output_ind,1) = -1;
                    shtr3(output_ind,1) = -1;
                    zout(output_ind, 1) = p;
                    xout(output_ind, 1) = i;
                    yout(output_ind, 1) = o;
                    output_ind = output_ind + 1;
                    
                    stage.SetPosition(setX,setY, zinstr(zcal,1));
                    tic;
                    while toc < (0.1)
                    end
                    
                    [i o p] = stage.GetPosition;
                    results_grid(1,result_ind) = i;
                    results_grid(2,result_ind) = o;
                    results_grid(3,result_ind) = p;
                    result_ind = result_ind + 1
                    
                    for zpo = 1:gridFrame
                        [i o p] = stage.GetPosition;
                        zout(output_ind, 1) = p;
                        xout(output_ind, 1) = i;
                        yout(output_ind, 1) = o;
                        output_ind = output_ind + 1;
                    end
                    %take 20 frames at 20 ms exposure time.
                    queueOutputData(s,[a a]);
                    tic;
                    while toc < (delayGrid/100)
                    end
                    s.startBackground;
                    while(~s.IsDone)
                        wait(s);
                    end
                    tic;
                    while toc < (delayGrid/100)
                    end
                    %end frame collection
                    
                end
                zcal = zcal + 1;
                b = 1;
                continue
            end
            
            if(b == 1)
                [o, backwardsz] = size(backwardarr);
                for f = 1:backwardsz
                    setX = backwardarr(1,f);
                    setY = backwardarr(2,f);
                    [i, o, p] = stage.GetPosition;
                    shtr1(output_ind,1) = -1; % indicate transition in output file
                    shtr2(output_ind,1) = -1;
                    shtr3(output_ind,1) = -1;
                    xout(output_ind, 1) = i;
                    yout(output_ind, 1) = o;
                    zout(output_ind, 1) = p;
                    output_ind = output_ind + 1;
                    %            stage.xset(setX);
                    %            stage.yset(setY);
                    stage.SetPosition(setX,setY, zinstr(zcal,1));
                    tic;
                    while toc < (0.1)
                    end
                    [i, o, p] = stage.GetPosition;
                    results_grid(1,result_ind) = i;
                    results_grid(2,result_ind) = o;
                    results_grid(3,result_ind) = p;
                    
                    result_ind = result_ind + 1
                    
                    for zpo = 1:gridFrame
                        [i, o, p] = stage.GetPosition;
                        zout(output_ind, 1) = p;
                        xout(output_ind, 1) = i;
                        yout(output_ind, 1) = o;
                        output_ind = output_ind + 1;
                    end
                    
                    %take 20 frames at 20 ms exposure time.
                    queueOutputData(s,[a a]);
                    tic;
                    while toc < (delayGrid/100)
                    end
                    s.startBackground;
                    while(~s.IsDone)
                        wait(s);
                    end
                    tic;
                    while toc < (delayGrid/100)
                    end
                    %end frame collection
                end
                zcal = zcal + 1;
                b = 0;
            end
            
        end
        
        output = horzcat(shtr1, shtr2, shtr3, xout, yout, zout); % concatonate shutter and z position data
        output_filename = 'output_controlpoint'; % NAME THE OUTPUT FILE
        
        %% output actual locations of nanostage to .txt file
        cd(FolderName);
        c= clock;
        fix(c);
        if(exist([output_filename '.dat'], 'file')==2) % if file exists, create new file
            output_filename = [output_filename '_' int2str(c(1)) int2str(c(2)) int2str(c(3)) '_' int2str(c(4)) int2str(c(5)) int2str(c(6))];
        end
        fid = fopen([output_filename '.dat'], 'a');
        fprintf(fid, '%d %d %d %f %f %f\r\n', output');
        fclose(fid);
        cd('Z:\Eugene\development R161_FromLidkeLab');
        set(handles.sequence_txt, 'String', ['SEQUENCE STATUS: ' 'Sequence Not Running'],'fontweight','normal');
        start(tmr); % resume timer to refresh positions
    end

end
