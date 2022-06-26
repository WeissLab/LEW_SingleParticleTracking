%% LEW SingleParticleTracking 4
% Lucien Weiss - WeissLab.ca 2021
% Created in the MoernerLab 2015
% Keyboard Shortcuts:
% Right Arrow - forward one frame
% Left Arrow - back one frame
% Up Arrow - Change Molecule
% Down Arrow - Change Molecule
% a - Autorun backward
% s - Stop Autorun
% d - Autorun forward
% p - cycle through plotting (all data, recent data, and hiding data)
% t - show/hide data - not yet incorporated]
% n - new molecule
% x - toggle tracking mode;


function [] = LEW_MULTICOLOR_TRACKING_V4 %#ok<FNDEF>
%% LEW_MULTICOLOR_TRACKING_v4
% Purpose: This code enables one channel to be displayed and tracked while
% showing a single image from a separate file in the other two channels.
% Version 4 changes the outputs from a single matrix to a structure with an
% element for each molecule. It also more compact in the sense that there
% is only one figure window.
figure(321)
close;
warning off all

%% VERSION INFO
    Version = '4.6';
    % Version 4.6
    % Uplodated several functions
    
    % Version 4.5
    % n - asymmetric rotating fit
    
    % Version 4.4
    % n - new molecule
    % x - toggle tracking mode;
    
    % Version 4.3
    % Pause on bad localization
    % Stop Tracking on bad localization
    
    % Version 4.2
        % fixed DHPSF
    % Version 4.1
        % fixed values showing up (position units need to be listed before position)
        % fixed scaling issue (would only change upon loading the movie before)
        % fixed clearing frames
        % changed key shortcuts

%% STORED VARIABLES
    DATA = struct;
    list_of_molecules = {};
    current_molecule=[];
    current_molecule_name='';
    FILEDATA = struct;
    max_frames = 1;
    scaling = [0 1; 0 1;0 1];     
    roi = [0,0];
    frame_limits = [1 inf];
    amplitude_limits = [0 2^16];
    amplitude_range = diff(amplitude_limits);
    background_limits = [0 2^16];
    background_range = diff(background_limits);
    sigma_limits = [.1 5];
    sigma_range = diff(sigma_limits);
    
%% PLACEHOLDER VARIABLES
    frame = 1;
    lower_frame_bin = 1;
    upper_frame_bin = 2;
    color_channel_loaded = [0 0 0]; 
    red_channel_image_scaled = [0 1; 0 0];
    green_channel_image_scaled = [];
    blue_channel_image_scaled = [];
    toggle_autorun = [0 0];
    toggle_plotting = 1;
    toggle_show_display = 1;
    Image_Region = [];
    
%% PROGRAM SETTINGS
%     fitting_options = optimset('FunValCheck','on','MaxIter', 1000 , 'TolFun',1e-4,'TolX',1e-4,'Display','none','Diagnostics','off','Jacobian','on');%optimset('FunValCheck','on');
    fitting_options = optimset('FunValCheck','on','Diagnostics','off','MaxIter', 1000 ,'Jacobian','off','Hessian','on', 'Display', 'off','TolFun',1e-4,'TolX',1e-4);
    fitting_region = 5; % 11x11 'UseParallel','off'
    border_region = ones(fitting_region*2+1); border_region(2:end-1,2:end-1)=nan;
    [regional_indices_row,regional_indices_col] = ndgrid(-fitting_region:fitting_region,-fitting_region:fitting_region);
    

%% GUI SETTINGS
    BG.GUI = [.9 .9 .9];
    BG.CONTROL = [.9 .9 .9];
    BG.FRAME = [.9 .9 .9];
    BG.DATA = [.9 .9 .9];
    BG.LOCALIZATION = [.9 .9 .9];
    FP.DATA_BACKGROUND = [.9 .9 .9];
    FP.DATA = [0 0 0];
    list_of_tracking_modes = {'Symmetric Gaussian';'Restricted Symmetric Gaussian';'Asymmetric Gaussian';'Asymmetric Rotating Gaussian';'Two Gaussians'};
    list_of_variables_for_hist = {'row','col','amplitude','sigma','background','signal_above_BG','fitted_signal','residual','localization_info','frame_cleared','frame_locked','frame_enabled'};

%% CREATE GUI
    TRACKING_FIGURE = figure(321);
    set(TRACKING_FIGURE,'position',[100 100 800 500],'Color',BG.GUI,'name','LEW TRACKING CODE','NumberTitle','off','WindowKeyPressFcn',{@KEYPRESS});
    set(0,'defaultuicontrolunits','normalized');
    % Control Panel
    CONTROL_PANEL = uipanel('parent',TRACKING_FIGURE,'BackgroundColor',BG.CONTROL,'units','normalized','Position',[0.01 0.01 .23 .98]);
        % Load Data
        uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','String','Load Data','units','normalized','Position', [0,.95,.5,.05],'callback',{@LOAD_DATA});
        function LOAD_DATA(~,~)
            [Loaded_File_Name, Loaded_File_Path] = uigetfile('*_LEW_tracking_v*.mat');
            Loaded_Data = load([Loaded_File_Path Loaded_File_Name]);
            DATA = Loaded_Data.DATA;
            set(tracking_mode_popup,'value',Loaded_Data.current_tracking_mode);
            color_channel_loaded = Loaded_Data.color_channel_loaded;
            FILEDATA.red_channel_filename = Loaded_Data.FILEDATA_filenames.red_channel_filename;
            FILEDATA.red_channel_filelocation = Loaded_File_Path;
            set(loadredchannel,'string',FILEDATA.red_channel_filename(max(1,end-19):end-4));
            FILEDATA.red_channel_fileinfo = imfinfo([Loaded_File_Path FILEDATA.red_channel_filename]);
            max_frames = length(FILEDATA.red_channel_fileinfo);
            if color_channel_loaded(2) == 1
                FILEDATA.green_channel_filename = Loaded_Data.FILEDATA_filenames.green_channel_filename;
                FILEDATA.green_channel_filelocation = Loaded_File_Path;
                set(loadgreenchannel,'string',FILEDATA.green_channel_filename(max(1,end-19):end-4));
                FILEDATA.green_channel_fileinfo = imfinfo([Loaded_File_Path FILEDATA.green_channel_filename]);
                green_channel_image=single(imread([FILEDATA.green_channel_filelocation FILEDATA.green_channel_filename],'Index',1,'Info',FILEDATA.green_channel_fileinfo));
                green_channel_image_scaled = mat2gray(green_channel_image,double([(1+scaling(2,1))*min(green_channel_image(:)), max(green_channel_image(:))*scaling(2,2)]));

            end
            if color_channel_loaded(3) == 1
                FILEDATA.blue_channel_filename = Loaded_Data.FILEDATA_filenames.blue_channel_filename;
                FILEDATA.blue_channel_filelocation = Loaded_File_Path;
                set(loadbluechannel,'string',FILEDATA.blue_channel_filename(max(1,end-19):end-4));
                FILEDATA.blue_channel_fileinfo = imfinfo([Loaded_File_Path FILEDATA.blue_channel_filename]);
                blue_channel_image=single(imread([FILEDATA.blue_channel_filelocation FILEDATA.blue_channel_filename],'Index',1,'Info',FILEDATA.blue_channel_fileinfo));
                blue_channel_image_scaled = mat2gray(blue_channel_image,double([(1+scaling(3,1))*min(blue_channel_image(:)), max(blue_channel_image(:))*scaling(3,2)]));
            end
            INITIALIZE_IMAGE;
            list_of_molecules = Loaded_Data.list_of_molecules;
            current_molecule = 1;
            current_molecule_name = 'track 1';
            Unselected = list_of_molecules;
            Unselected(any(cellfun(@isempty,Unselected),2),:) = [];
            set(Unselected_molecules,'String',Unselected,'Value',1)
            set(Current_molecule_popup,'string',list_of_molecules,'value',1);
        end

        % Save Data
        uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','String','Save Data','Position', [.5,.95,.5,.05],'units','normalized','callback',{@SAVE_DATA});
        function SAVE_DATA(~,~)
            current_time = clock;
            imagename = FILEDATA.red_channel_filename(1:end-4);
            savetime = [num2str(current_time(1)) num2str(current_time(2), '%02i') num2str(current_time(3), '%02i') '_' num2str(mod(current_time(4)-1,12)+1, '%02i') '.' num2str(current_time(5),'%02i')];
            filename = [imagename '_saved_at_' savetime '_LEW_tracking_v' Version '.mat'];
            FILEDATA_filenames.red_channel_filename = FILEDATA.red_channel_filename;
            if color_channel_loaded(2) == 1
                FILEDATA_filenames.green_channel_filename = FILEDATA.green_channel_filename;
            end
            if color_channel_loaded(3) == 1
                FILEDATA_filenames.blue_channel_filename = FILEDATA.blue_channel_filename; 
            end
            current_tracking_mode = get(tracking_mode_popup,'value'); 
            save(filename,'DATA','FILEDATA_filenames','current_tracking_mode','color_channel_loaded','list_of_molecules');
        end
        % Image Read
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Load Images','Position', [0,.9,.63,.04],'units','normalized','backgroundcolor',BG.CONTROL);
        loadredchannel = uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','String','Load Movie (red)','Position', [0,.85,.6,.05],'units','normalized','callback',{@LOAD_RED});
        function  LOAD_RED(~,~)
            [RED_Movie_Name, RED_Movie_Path] = uigetfile({'*.tif';'*.*'},'Select movie for tracking (Red channel');
            if RED_Movie_Path == 0
                return;
            end
            color_channel_loaded(1)=1;
            FILEDATA.red_channel_filename = RED_Movie_Name;
            FILEDATA.red_channel_filelocation = RED_Movie_Path;
            set(loadredchannel,'string',RED_Movie_Name(max(1,end-19):end-4));
            FILEDATA.red_channel_fileinfo = imfinfo([RED_Movie_Path RED_Movie_Name]);
            max_frames = length(FILEDATA.red_channel_fileinfo);
            INITIALIZE_IMAGE;
            ADD_MOLECULE;
        end
        loadgreenchannel = uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','String','Load Image (green)','Position', [0,.80,.6,.05],'units','normalized','callback',{@LOAD_GREEN});
        function  LOAD_GREEN(~, ~)
            [GREEN_Movie_Name, GREEN_Movie_Path] = uigetfile({'*.tif';'*.*'},'Select Green Channel');
            if GREEN_Movie_Path == 0
                return;
            end
            color_channel_loaded(2)=1;
            FILEDATA.green_channel_filename = GREEN_Movie_Name;
            FILEDATA.green_channel_filelocation = GREEN_Movie_Path;
            set(loadgreenchannel,'string',GREEN_Movie_Name(max(1,end-19):end-4));
            FILEDATA.green_channel_fileinfo = imfinfo([GREEN_Movie_Path GREEN_Movie_Name]);
            green_channel_image=single(imread([FILEDATA.green_channel_filelocation FILEDATA.green_channel_filename],'Index',1,'Info',FILEDATA.green_channel_fileinfo));
            green_channel_image_scaled = mat2gray(green_channel_image,double([(1+scaling(2,1))*min(green_channel_image(:)), max(green_channel_image(:))*scaling(2,2)]));
            INITIALIZE_IMAGE;
        end
        loadbluechannel = uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','String','Load Image (blue)','Position', [0,.75,.6,.05],'units','normalized','callback',{@LOAD_BLUE});
        function  LOAD_BLUE(~, ~)
            % Sets the green channel image.
            [BLUE_Movie_Name, BLUE_Movie_Path] = uigetfile({'*.tif';'*.*'},'Select Blue Channel');
            if BLUE_Movie_Path == 0
                return;
            end
            color_channel_loaded(3)=1;
            FILEDATA.blue_channel_filename = BLUE_Movie_Name;
            FILEDATA.blue_channel_filelocation = BLUE_Movie_Path;
            set(loadbluechannel,'string',BLUE_Movie_Name(max(1,end-19):end-4));
            FILEDATA.blue_channel_fileinfo = imfinfo([BLUE_Movie_Path BLUE_Movie_Name]);
            blue_channel_image=single(imread([FILEDATA.blue_channel_filelocation FILEDATA.blue_channel_filename],'Index',1,'Info',FILEDATA.blue_channel_fileinfo));
            blue_channel_image_scaled = mat2gray(blue_channel_image,double([(1+scaling(3,1))*min(blue_channel_image(:)), max(blue_channel_image(:))*scaling(3,2)]));
            INITIALIZE_IMAGE;
        end
        % Color Scaling
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Scaling','Position', [.63,.9,.37,.04],'units','normalized','backgroundcolor',BG.CONTROL);
        cp.red.min = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',scaling(1,1),'Position', [.63,.85,.15,.05],'units','normalized','callback',{@COLOR_SET});
        cp.green.min  = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',scaling(2,1),'Position', [.63,.80,.15,.05],'units','normalized','callback',{@COLOR_SET});
        cp.blue.min  = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',scaling(3,1),'Position', [.63,.75,.15,.05],'units','normalized','callback',{@COLOR_SET});
        cp.red.max  = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',scaling(1,2),'Position', [.81,.85,.15,.05],'units','normalized','callback',{@COLOR_SET});
        cp.green.max  = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',scaling(2,2),'Position', [.81,.80,.15,.05],'units','normalized','callback',{@COLOR_SET});
        cp.blue.max  = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',scaling(3,2),'Position', [.81,.75,.15,.05],'units','normalized','callback',{@COLOR_SET});
        function COLOR_SET(~,~)
            % Images are scaled from the minimum value to (max value/scaling factor)
            scaling(1,1) = str2double(get(cp.red.min,'string'));
            scaling(2,1) = str2double(get(cp.green.min,'string'));
            scaling(3,1) = str2double(get(cp.blue.min,'string'));
            scaling(1,2) = str2double(get(cp.red.max,'string'));
            scaling(2,2) = str2double(get(cp.green.max,'string'));
            scaling(3,2) = str2double(get(cp.blue.max,'string'));
            % Need to reset images
            if color_channel_loaded(1) == 1
                red_channel_image_scaled = mat2gray(red_channel_image_scaled,double([(1+scaling(1,1))*min(red_channel_image_scaled(:)), max(red_channel_image_scaled(:))*scaling(1,2)]));
            end
            if color_channel_loaded(2) == 1
                green_channel_image_scaled = mat2gray(green_channel_image_scaled,double([(1+scaling(2,1))*min(green_channel_image_scaled(:)), max(green_channel_image_scaled(:))*scaling(2,2)]));
            end
            if color_channel_loaded(3) == 1
                blue_channel_image_scaled = mat2gray(blue_channel_image_scaled,double([(1+scaling(3,1))*min(blue_channel_image_scaled(:)), max(blue_channel_image_scaled(:))*scaling(3,2)]));
            end
            UPDATE_IMAGE;
        end
        % Fit Type
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Fit type','Position', [0,.7,.3,.05],'units','normalized','backgroundcolor',BG.CONTROL);
        tracking_mode_popup = uicontrol('parent',CONTROL_PANEL,'Style','popup','String',list_of_tracking_modes,'Position', [.3,.7,.7,.05],'units','normalized','callback',{@CHANGE_TRACKING_MODE});
        function CHANGE_TRACKING_MODE(~,~)
            switch char(list_of_tracking_modes(get(tracking_mode_popup,'value')))
                case 'Symmetric Gaussian'
                    list_of_variables_for_hist = {'row','col','sigma','amplitude','background','signal_above_BG','fitted_signal','residual','localization_info','frame_cleared','frame_locked','frame_enabled'};
                    sigma_limits = [.1 3];
                    fitting_region = 5; % 11x11
                    
                case 'Restricted Symmetric Gaussian'
                    list_of_variables_for_hist = {'row','col','sigma','amplitude','background','signal_above_BG','fitted_signal','residual','localization_info','frame_cleared','frame_locked','frame_enabled'};
                    sigma_limits = [.1 3];
                    fitting_region = 4; % 9x9
                case 'Asymmetric Gaussian'
                    list_of_variables_for_hist = {'row','col','sigma_row','sigma_col','amplitude','background','z','signal_above_BG','fitted_signal','residual','localization_info','frame_cleared','frame_locked','frame_enabled'};
                    sigma_limits = [.1 5];
                    fitting_region = 7; % 15x15
                case 'Asymmetric Rotating Gaussian'
                    list_of_variables_for_hist = {'row','col','sigma_row','sigma_col','theta','amplitude','background','z','signal_above_BG','fitted_signal','residual','localization_info','frame_cleared','frame_locked','frame_enabled'};
                    sigma_limits = [.25 4];
                    fitting_region = 7; % 15x15
                case 'Two Gaussians'
                    sigma_limits = [.1 3];
                    list_of_variables_for_hist = {'row','col','row_1','col_1','row_2','col_2','sigma_1','sigma_2','amplitude_1','amplitude_2','background','theta','z','signal_above_BG','fitted_signal','residual','localization_info','frame_cleared','frame_locked','frame_enabled'};
                    fitting_region = 9; % 19x19
            end
            border_region = ones(fitting_region*2+1); border_region(2:end-1,2:end-1)=nan;
            set(lp.raw_data_axes,'xlim',[.5 fitting_region*2+1.5],'ylim',[.5 fitting_region*2+1.5]);
            set(lp.fit_data_axes,'xlim',[.5 fitting_region*2+1.5],'ylim',[.5 fitting_region*2+1.5]);
            set(lp.residual_data_axes,'xlim',[.5 fitting_region*2+1.5],'ylim',[.5 fitting_region*2+1.5]);
            set(cp.sigma_min,'string',sigma_limits(1));
            set(cp.sigma_max,'string',sigma_limits(2));
            CHANGE_LIMITS(0,0);
            [regional_indices_row,regional_indices_col] = ndgrid(-fitting_region:fitting_region,-fitting_region:fitting_region);
            set(dp.data_selection,'string',list_of_variables_for_hist);
        end
        % Localization
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Parameter','Position', [0,.68,.4,.03],'units','normalized','backgroundcolor',BG.CONTROL);
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Min','Position', [.4,.68,.3,.03],'units','normalized','backgroundcolor',BG.CONTROL);
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Max','Position', [.7,.68,.3,.03],'units','normalized','backgroundcolor',BG.CONTROL);
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Frame','Position', [0,.65,.4,.03],'units','normalized','backgroundcolor',BG.CONTROL);
        cp.frame_min = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',frame_limits(1),'Position', [.4,.65,.3,.03],'units','normalized','callback',{@CHANGE_LIMITS});
        cp.frame_max = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',frame_limits(2),'Position', [.7,.65,.3,.03],'units','normalized','callback',{@CHANGE_LIMITS});
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Amplitude','Position', [0,.62,.4,.03],'units','normalized','backgroundcolor',BG.CONTROL);
        cp.amplitude_min = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',amplitude_limits(1),'Position', [.4,.62,.3,.03],'units','normalized','callback',{@CHANGE_LIMITS});
        cp.amplitude_max = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',2^16,'Position', [.7,.62,.3,.03],'units','normalized','callback',{@CHANGE_LIMITS});
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Sigma','Position', [0,.59,.4,.03],'units','normalized','backgroundcolor',BG.CONTROL);
        cp.sigma_min = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',sigma_limits(1),'Position', [.4,.59,.3,.03],'units','normalized','callback',{@CHANGE_LIMITS});
        cp.sigma_max = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',sigma_limits(2),'Position', [.7,.59,.3,.03],'units','normalized','callback',{@CHANGE_LIMITS});
        uicontrol('parent',CONTROL_PANEL,'Style','text','String','Background','Position', [0,.56,.4,.03],'units','normalized','backgroundcolor',BG.CONTROL);
        cp.background_min = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',background_limits(1),'Position', [.4,.56,.3,.03],'units','normalized','callback',{@CHANGE_LIMITS});
        cp.background_max = uicontrol('parent',CONTROL_PANEL,'Style','edit','String',2^16,'Position', [.7,.56,.3,.03],'units','normalized','callback',{@CHANGE_LIMITS});
        function CHANGE_LIMITS(~,~)
            frame_limits  = [str2double(get(cp.frame_min,'string')), str2double(get(cp.frame_max,'string'))];
            DATA.(genvarname(current_molecule_name)).min_frame_allowed = frame_limits(1); %#ok<*DEPGENAM>
            DATA.(genvarname(current_molecule_name)).max_frame_allowed = frame_limits(2);            
            amplitude_limits = [str2double(get(cp.amplitude_min,'string')), str2double(get(cp.amplitude_max,'string'))];
            amplitude_range = diff(amplitude_limits);
            background_limits = [str2double(get(cp.background_min,'string')), str2double(get(cp.background_max,'string'))];
            background_range = diff(background_limits);
            sigma_limits = [str2double(get(cp.sigma_min,'string')), str2double(get(cp.sigma_max,'string'))];
            sigma_range = diff(sigma_limits);
        end
        uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','Position', [0,.50,1,.05],'units','normalized','value',0,'string','Clear Frames Outside Limits','callback',{@CLEAR_FRAMES_OUTSIDE_LIMITS});
        cp.toggle_tracking = uicontrol('parent',CONTROL_PANEL,'Style','togglebutton','Position', [0,.35,1,.05],'units','normalized','value',1,'string','Toggle track mode');

        function CLEAR_FRAMES_OUTSIDE_LIMITS(~,~)
            outside_frame_limits = find(DATA.(genvarname(current_molecule_name)).frame<frame_limits(1) | DATA.(genvarname(current_molecule_name)).frame>frame_limits(2));
            outside_background_limits = find(DATA.(genvarname(current_molecule_name)).background<background_limits(1) | DATA.(genvarname(current_molecule_name)).background>background_limits(2));
            switch char(list_of_tracking_modes(get(tracking_mode_popup,'value')))
                case 'Symmetric Gaussian'
                    outside_amplitude_limits = find(DATA.(genvarname(current_molecule_name)).amplitude<amplitude_limits(1) | DATA.(genvarname(current_molecule_name)).amplitude>amplitude_limits(2));
                    outside_sigma_limits = find(DATA.(genvarname(current_molecule_name)).sigma<sigma_limits(1) | DATA.(genvarname(current_molecule_name)).sigma>sigma_limits(2));
                case 'Restricted Symmetric Gaussian'
                    outside_amplitude_limits = find(DATA.(genvarname(current_molecule_name)).amplitude<amplitude_limits(1) | DATA.(genvarname(current_molecule_name)).amplitude>amplitude_limits(2));
                    outside_sigma_limits = find(DATA.(genvarname(current_molecule_name)).sigma<sigma_limits(1) | DATA.(genvarname(current_molecule_name)).sigma>sigma_limits(2));
                case 'Asymmetric Gaussian'
                    outside_amplitude_limits = find(DATA.(genvarname(current_molecule_name)).amplitude<amplitude_limits(1) | DATA.(genvarname(current_molecule_name)).amplitude>amplitude_limits(2));
                    outside_sigma_limits = [find(DATA.(genvarname(current_molecule_name)).sigma_row<sigma_limits(1) | DATA.(genvarname(current_molecule_name)).sigma_row>sigma_limits(2)), ...
                                            find(DATA.(genvarname(current_molecule_name)).sigma_col<sigma_limits(1) | DATA.(genvarname(current_molecule_name)).sigma_col>sigma_limits(2))];
                case 'Asymmetric Rotating Gaussian'
                    outside_amplitude_limits = find(DATA.(genvarname(current_molecule_name)).amplitude<amplitude_limits(1) | DATA.(genvarname(current_molecule_name)).amplitude>amplitude_limits(2));
                    outside_sigma_limits = [find(DATA.(genvarname(current_molecule_name)).sigma_row<sigma_limits(1) | DATA.(genvarname(current_molecule_name)).sigma_row>sigma_limits(2)), ...
                                            find(DATA.(genvarname(current_molecule_name)).sigma_col<sigma_limits(1) | DATA.(genvarname(current_molecule_name)).sigma_col>sigma_limits(2))];
                
                case 'Two Gaussians'
                    outside_amplitude_limits = [find(DATA.(genvarname(current_molecule_name)).amplitude_1<amplitude_limits(1) | DATA.(genvarname(current_molecule_name)).amplitude_1>amplitude_limits(2)), ...
                                                find(DATA.(genvarname(current_molecule_name)).amplitude_2<amplitude_limits(1) | DATA.(genvarname(current_molecule_name)).amplitude_2>amplitude_limits(2))];
                    outside_sigma_limits = [find(DATA.(genvarname(current_molecule_name)).sigma_1<sigma_limits(1) | DATA.(genvarname(current_molecule_name)).sigma_1>sigma_limits(2)), ...
                                            find(DATA.(genvarname(current_molecule_name)).sigma_2<sigma_limits(1) | DATA.(genvarname(current_molecule_name)).sigma_2>sigma_limits(2))];
            end
            outside_limits = unique([outside_frame_limits,outside_amplitude_limits,outside_background_limits,outside_sigma_limits]);
            CLEAR_FRAME(0,0,outside_limits)
        end


        %set_roi
        uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','String','Set ROI','Position', [0,.45,.25,.05],'units','normalized','callback',{@SET_ROI});
        %delete_outer_locs 
        uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','String',{'Delete Outsiders'},'Position', [.25,.45,.42,.05],'units','normalized','callback',{@DELETE_OUTSIDERS},'fontsize',7);
        %view_roi
        uicontrol('parent',CONTROL_PANEL,'Style','togglebutton','String','Show ROI','Position', [.67,.45,.33,.05],'units','normalized','callback',{@SHOW_ROI});
        function [] = SET_ROI(~,~)
            set(gcf,'CurrentAxes',fp.current_frame_axes)
            [~, x, y] = roipoly;
            roi=[x y];
            set(fp.plot_roi,'XData',roi(:,1),'YData', roi(:,2),'visible','on');
        end
        function [] = DELETE_OUTSIDERS(~,~)
            inside_roi = inpolygon(DATA.(genvarname(current_molecule_name)).col,DATA.(genvarname(current_molecule_name)).row,roi(:,1),roi(:,2));
            outside_roi=(inside_roi-1)*-1;
            fs=find(outside_roi.*(1:numel(outside_roi))>0);
            CLEAR_FRAME(0,0,fs)
        end
        function [] = SHOW_ROI(~,~)
            if strcmp(get(fp.plot_roi,'visible'),'off')==1
                set(fp.plot_roi,'visible','on');
            else
                set(fp.plot_roi,'visible','off');
            end
        end


        % clear_frame
        uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','Position', [0,.40,.5,.05],'units','normalized','value',0,'string','Clear Frame','callback',{@CLEAR_FRAME,NaN});
        function CLEAR_FRAME(~,~,frame_to_clear)
            if isnan(frame_to_clear)
                frame_to_clear = frame;
            end
            DATA.(genvarname(current_molecule_name)).frame_cleared(frame_to_clear) = 1;
            DATA.(genvarname(current_molecule_name)).row(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).col(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).row_1(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).col_1(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).row_2(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).col_2(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).amplitude(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).amplitude_1(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).amplitude_2(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).sigma(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).sigma_row(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).sigma_col(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).sigma_1(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).sigma_2(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).background(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).theta(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).z(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).signal_above_BG(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).fitted_signal(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).residual(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).localization_info(frame_to_clear) = NaN;
            DATA.(genvarname(current_molecule_name)).frame_enabled(frame_to_clear) = 0;
            UPDATE_IMAGE;
            DISPLAY_FITTING_HISTOGRAM
        end

        % lock_frame
        uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','Position', [.5,.40,.5,.05],'units','normalized','value',0,'string','Lock Frame','callback',{@LOCK_FRAME});
        function LOCK_FRAME(~,~)
            DATA.(genvarname(current_molecule_name)).frame_locked(frame) = 1;
            DATA.(genvarname(current_molecule_name)).frame_enabled(frame) = 0;
        end

        % Add Molecule
        uicontrol('parent',CONTROL_PANEL,'Style','pushbutton','String','Add molecule','Position', [0,.25,.5,.05],'units','normalized','callback',{@ADD_MOLECULE});
        function ADD_MOLECULE(~,~)
            molecule_name = ['track ' num2str(size(list_of_molecules,1)+1)];
            list_of_molecules = [list_of_molecules; molecule_name];             
            DATA.(genvarname(molecule_name)).name = molecule_name;
            current_time = clock;
            DATA.(genvarname(molecule_name)).date_created = [num2str(current_time(1)) '/' num2str(current_time(2),'%02i') '/' num2str(current_time(3),'%02i') ' ' num2str(mod(current_time(4)-1,12)+1,'%02i') ':' num2str(current_time(5),'%02i')];
            DATA.(genvarname(molecule_name)).saved = 0;
            DATA.(genvarname(molecule_name)).deleted = 0;
            DATA.(genvarname(molecule_name)).localization_warning = 0;
            DATA.(genvarname(molecule_name)).frame = (1:max_frames);
            DATA.(genvarname(molecule_name)).row = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).col = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).row_1 = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).col_1 = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).row_2 = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).col_2 = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).amplitude = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).amplitude_1 = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).amplitude_2 = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).sigma = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).sigma_row = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).sigma_col = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).sigma_1 = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).sigma_2 = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).background = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).theta = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).z = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).signal_above_BG = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).fitted_signal = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).residual = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).localization_info = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).frame_cleared = NaN(1,max_frames);
            DATA.(genvarname(molecule_name)).frame_locked = zeros(1,max_frames);
            DATA.(genvarname(molecule_name)).frame_enabled = zeros(1,max_frames);
            DATA.(genvarname(molecule_name)).min_frame_allowed = 1;
            DATA.(genvarname(molecule_name)).max_frame_allowed = max_frames;
            DATA.(genvarname(molecule_name)).min_frame_used = NaN;
            DATA.(genvarname(molecule_name)).max_frame_used = NaN;

            Unselected = [cellstr(get(Unselected_molecules,'String')); molecule_name];
            Unselected(any(cellfun(@isempty,Unselected),2),:) = [];
            set(Unselected_molecules,'String',Unselected,'Value',1)
            set(Current_molecule_popup,'string',list_of_molecules,'value',size(list_of_molecules,1));
            SELECT_MOLECULE;
        end
        Current_molecule_popup = uicontrol('parent',CONTROL_PANEL,'Style','popup','String',list_of_molecules,'Position', [.5,.25,.5,.05],'units','normalized','callback',{@SELECT_MOLECULE});
        function SELECT_MOLECULE(~,~)
            Current_value = get(Current_molecule_popup,'value');
            Current_selection = char(list_of_molecules(Current_value));
            track_number = sscanf(Current_selection,'%s%d');
            current_molecule_name = Current_selection;
            current_molecule = track_number(end);
            set(cp.frame_min,'string',DATA.(genvarname(current_molecule_name)).min_frame_allowed);
            set(cp.frame_max,'string',DATA.(genvarname(current_molecule_name)).max_frame_allowed);
            UPDATE_IMAGE
        end

        % Select Molecules
        uicontrol('parent',CONTROL_PANEL,'style','text','position',[0  .2 .5 .05],'units','normalized','string','Not tracking','fontsize',10);
        uicontrol('parent',CONTROL_PANEL,'style','text','position',[.5 .2 .5 .05],'units','normalized','string','Tracking','fontsize',10);
        Unselected_molecules = uicontrol('parent',CONTROL_PANEL,'style','listbox','position',[0 0 .5 .2],'units','normalized','max',1,'Callback', {@TRACK_MOLECULES,1});
        Selected_molecules   = uicontrol('parent',CONTROL_PANEL,'style','listbox','position',[.5 0 .5 .2],'units','normalized','max',1,'Callback', {@TRACK_MOLECULES,-1});
        function TRACK_MOLECULES(~,~,input)
            if strcmp(get(TRACKING_FIGURE,'SelectionType'),'open')
                Unselected = cellstr(get(Unselected_molecules,'String'));
                Selected = cellstr(get(Selected_molecules, 'string'));
                switch input 
                    case 1
                        clicked = get(Unselected_molecules,'Value');
                        Current_selection = Unselected(clicked);
                        Unselected(clicked,:) = [];
                        Selected = [Selected; Current_selection];
                        track_number = sscanf(char(Current_selection),'%s%d');
                        set(Current_molecule_popup,'string',list_of_molecules,'value',track_number(end));
                        SELECT_MOLECULE;
                    case -1
                        clicked = get(Selected_molecules,'Value');
                        Current_selection = Selected(clicked);
                        Unselected = [Unselected; Current_selection];
                        Selected(clicked,:) = [];
                end
                Unselected(any(cellfun(@isempty,Unselected),2),:) = [];
                Selected(any(cellfun(@isempty,Selected),2),:) = [];
                set(Unselected_molecules,'String', sort(Unselected),'Value',1)
                set(Selected_molecules,'String', sort(Selected),'Value',1)
            end
        end

    % Current Frame
    FRAME_PANEL = uipanel('parent',TRACKING_FIGURE,'BackgroundColor',BG.FRAME,'units','normalized','Position',[0.26 0.26 .48 .63]);
        fp.current_frame_axes = axes('Units','normal', 'Position', [0 0 1 1], 'Parent', FRAME_PANEL,'HitTest','on','ButtonDownFcn',@CLICK);
        fp.current_frame = image(ones(10,10,3),'Parent',fp.current_frame_axes,'HitTest','on','ButtonDownFcn',@CLICK);
        set(fp.current_frame_axes,'XLimMode','auto','YLimMode','auto')
        hold on;
        fp.plot_positions = plot(NaN,NaN,'w','Parent',fp.current_frame_axes,'HitTest','on','ButtonDownFcn',@CLICK);
        fp.plot_recent_position = plot(NaN,NaN,'w','Parent',fp.current_frame_axes,'HitTest','off','visible','off');
        fp.scatter_positions = scatter(NaN,NaN,'y','Parent',fp.current_frame_axes,'HitTest','on','ButtonDownFcn',@CLICK);
        fp.plot_roi = plot(NaN,NaN,'w','Parent',fp.current_frame_axes,'HitTest','off','visible','off');
        
function INITIALIZE_IMAGE(~,~)
    if color_channel_loaded(1) == 1
        red_channel_image=single(imread([FILEDATA.red_channel_filelocation FILEDATA.red_channel_filename],'Index',frame,'Info',FILEDATA.red_channel_fileinfo));
        red_channel_image_scaled = mat2gray(red_channel_image,double([(1+scaling(1,1))*min(red_channel_image(:)), max(red_channel_image(:))*scaling(1,2)]));
    elseif color_channel_loaded(2) == 1
        red_channel_image_scaled = green_channel_image_scaled*0;
    elseif color_channel_loaded(3) == 1
        red_channel_image_scaled = blue_channel_image_scaled*0;
    end
    if  color_channel_loaded == [1 0 0]
        green_channel_image_scaled = red_channel_image_scaled*0;
        blue_channel_image_scaled = green_channel_image_scaled;
    elseif color_channel_loaded(2:3) == [1 0]
        blue_channel_image_scaled = red_channel_image_scaled*0;
    elseif color_channel_loaded(2:3) == [0 1]
        green_channel_image_scaled = red_channel_image_scaled*0;
    end
    set(fp.current_frame,'CData',cat(3,red_channel_image_scaled,green_channel_image_scaled,blue_channel_image_scaled));
    figure(321);
    axis tight;
end

    MOLECULE_FIT_PANEL = uipanel('parent',TRACKING_FIGURE,'BackgroundColor',BG.FRAME,'units','normalized','Position',[0.26 0.88 .48 .11]);
        fp.data_frame_text = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[0  .6 .2 .4],'units','normalized','string','frame','fontsize',10,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_frame_number = uicontrol('parent',MOLECULE_FIT_PANEL,'style','edit','position',[0  .0 .2 .6],'units','normalized','string',frame,'fontsize',10,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA,'callback',{@ENTER_FRAME});
        function ENTER_FRAME(~,~)
            CHANGE_FRAME(str2double(get(fp.data_frame_number,'string')));
        end
        fp.data_col_position_text = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.2  .6 .15 .4],'units','normalized','string','col','fontsize',10,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_col_position = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.2  .0 .15 .6],'units','normalized','string',NaN,'fontsize',10,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_row_position_text = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.35  .6 .15 .4],'units','normalized','string','row','fontsize',10,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_row_position = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.35  .0 .15 .6],'units','normalized','string',NaN,'fontsize',10,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_sigma_text = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.5  .6 .1 .4],'units','normalized','string','sigma','fontsize',10,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_sigma = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.5  .0 .1 .6],'units','normalized','string', num2str(NaN),'fontsize',10,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_amplitude_text = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.60  .6 .2 .4],'units','normalized','string','amplitude','fontsize',9,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_amplitude = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.60  .0 .2 .6],'units','normalized','string',NaN,'fontsize',9,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_background_text = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.80  .6 .20 .4],'units','normalized','string','background','fontsize',9,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);
        fp.data_background = uicontrol('parent',MOLECULE_FIT_PANEL,'style','text','position',[.80  .0 .20 .6],'units','normalized','string',NaN,'fontsize',9,'BackgroundColor',FP.DATA_BACKGROUND,'ForegroundColor',FP.DATA);

    % Localization Data
    DATA_PANEL = uipanel('parent',TRACKING_FIGURE,'BackgroundColor',BG.DATA,'units','normalized','Position',[0.26 0.01 .48 .23]);
        dp.data_plot_axes = axes('Parent', DATA_PANEL,'Units','normal', 'OuterPosition', [.3 0 .5 1]);
        dp.data_plot = plot([0 0 15 15 3 3 3 3 0 0 15 15 3 3 9 9 3 3 15 15 0 0],'Parent', dp.data_plot_axes);
        dp.data_selection = uicontrol('parent',DATA_PANEL,'style','popup','string',list_of_variables_for_hist,'position',[0,.8,.3,.1],'units','normalized','callback',{@DISPLAY_FITTING_HISTOGRAM});
        dp.data_min_max_frames_text = uicontrol('parent',DATA_PANEL,'style','text','string','Min/Max frame','units','normal','position',[0 .2 .3 .1],'BackgroundColor',FP.DATA_BACKGROUND);
        dp.data_min_max_frames = uicontrol('parent',DATA_PANEL,'style','text','string','NaN / NaN','units','normal','position',[0 .0 .3 .2],'BackgroundColor',FP.DATA_BACKGROUND);
        dp.data_hist_axes = axes('Parent', DATA_PANEL,'Units','normal', 'OuterPosition', [.75 0 .25 1]);
        dp.data_hist = barh([0 0],'Parent', dp.data_hist_axes);
        linkaxes([dp.data_hist_axes,dp.data_plot_axes],'y');
        set(gca,'yticklabel',[],'ylim',[-1 16]);

    % Localization Image
    LOCALIZATION_PANEL = uipanel('parent',TRACKING_FIGURE,'BackgroundColor',BG.LOCALIZATION,'units','normalized','Position',[0.76 0.01 .23 .98]);
        lp.raw_data_axes = axes('Units','normal', 'Position', [0 0.67 1 .3],'Parent', LOCALIZATION_PANEL);
        lp.raw_data = image(regional_indices_row,'Parent',lp.raw_data_axes,'HitTest','off');
        set(gca,'xtick',[])

        lp.fit_data_axes  = axes('Units','normal', 'Position', [0 0.34 1 .3],'Parent', LOCALIZATION_PANEL);
        lp.fit_data = image(regional_indices_row,'Parent',lp.fit_data_axes,'HitTest','off');
        set(gca,'xtick',[])

        lp.residual_data_axes     = axes('Units','normal', 'Position', [0 0.01 1 .3],'Parent', LOCALIZATION_PANEL);
        lp.residual_data = image(regional_indices_row,'Parent',lp.residual_data_axes,'HitTest','off');
        set(gca,'xtick',[])


%% CONTROL GUI
function KEYPRESS(~,eventdata)
    % this function controls all the keyboard shortcuts
    if numel(get(Selected_molecules,'string')) == 0
        frame_change = 3;
    else
        frame_change = 1;
    end
    try eventdata.Source.String
        switch eventdata.Source.String
            case 'Toggle track mode'
                eventdata.Character = 'x';
        end
    catch
    end
    
    switch eventdata.Character
        case 29 % Right Arrow
            CHANGE_FRAME(frame+1);
        case 28 % Left Arrow
            CHANGE_FRAME(frame-1);
        case 31 % Down Arrow (increase molecule number)
            entries = length(get(Current_molecule_popup,'String'));
            value   = get(Current_molecule_popup,'Value');
            if value<entries
                set(Current_molecule_popup,'Value',value+1);
                SELECT_MOLECULE;
            end
        case 30 % Up Arrow (decrease molecule number)
            value   = get(Current_molecule_popup,'Value');
            if value>1
                set(Current_molecule_popup,'Value',value-1);
                SELECT_MOLECULE;
            end
        case 'a' % Run Backward
            toggle_autorun = [mod(toggle_autorun(1)+1,2),0];
            AUTORUN('reverse',frame_change)
            
        case 's' % Run Forward
            toggle_autorun = [0, 0];
            drawnow;
        case 'd' % Run Forward
            toggle_autorun = [0, mod(toggle_autorun(2)+1,2)];
            AUTORUN('forward',frame_change)
        case 'p' % change plotting data
            toggle_plotting = mod(toggle_plotting+1,3);
            CHANGE_PLOTTING_DATA(toggle_plotting)
        case 't' % show/hide display data
            toggle_show_display = mod(toggle_show_display+1,2);
            set(Current_molecule_popup,'Value',toggle_show_display);
        case 'x'
            toggle_tracking_mode = mod(get(cp.toggle_tracking,'Value')+1,2);
            set(cp.toggle_tracking,'Value',toggle_tracking_mode);
        case 'n'
            ADD_MOLECULE(0,0);
    end
end

function [] = AUTORUN(direction,speed)
    switch direction
        case 'reverse'
            while toggle_autorun == [1, 0] 
                CHANGE_FRAME(frame-speed);
                drawnow;
            end
        case 'forward'
            while toggle_autorun == [0, 1] 
                CHANGE_FRAME(frame+speed);
                drawnow;
            end
    end
end

function [] = CHANGE_PLOTTING_DATA(plotting_mode)
    if plotting_mode == 0
        set(fp.plot_positions,'visible','off')
        set(fp.plot_recent_position,'visible','off')
    elseif plotting_mode == 1
        set(fp.plot_positions,'visible','on')
        set(fp.plot_recent_position,'visible','off')
    elseif plotting_mode == 2
        set(fp.plot_positions,'visible','off')
        set(fp.plot_recent_position,'visible','on')
    end
end

function [] = CHANGE_FRAME(frame_input)
    % This function changes the frame
    newframe = min(frame_input,max_frames);
    newframe = max(newframe,1);
    if frame ~= newframe
        frame = newframe;
        lower_frame_bin = max(frame-2,1);
        upper_frame_bin = min(frame+2,max_frames);
        PREPARE_TO_LOCALIZE
        UPDATE_IMAGE
    else
        toggle_autorun = [0 0];
    end
end
function CLICK(~,~)
    if isempty(current_molecule)
        ADD_MOLECULE;
    end
    molecule_name = current_molecule_name;
    Unselected = cellstr(get(Unselected_molecules,'String'));
    Selected = cellstr(get(Selected_molecules, 'string'));     
    new_clicked = find(strcmp(molecule_name,Unselected)==1, 1);
    if isempty(new_clicked)~=1
        Current_selection = Unselected(new_clicked);
        Unselected(new_clicked,:) = [];
        Selected = [Selected; Current_selection];
    end
    Unselected(any(cellfun(@isempty,Unselected),2),:) = [];
    Selected(any(cellfun(@isempty,Selected),2),:) = [];
    set(Unselected_molecules,'String', sort(Unselected),'Value',1)
    set(Selected_molecules,'String', sort(Selected),'Value',1)
    input_click = get(fp.current_frame_axes,'Currentpoint');
    DATA.(genvarname(molecule_name)).frame_cleared(frame) = 0;
    DATA.(genvarname(molecule_name)).frame_locked(frame) = 0;
    DATA.(genvarname(molecule_name)).frame_enabled(frame) = 1;
    switch char(list_of_tracking_modes(get(tracking_mode_popup,'value')))
        case 'Symmetric Gaussian'
            LOCALIZE_SYMMETRIC_GAUSSIAN(molecule_name,input_click([1,3]));
        case 'Restricted Symmetric Gaussian'
            LOCALIZE_RESTRICTED_GAUSSIAN(molecule_name,input_click([1,3]));
        case 'Asymmetric Gaussian'
            LOCALIZE_ASYMMETRIC_GAUSSIAN(molecule_name,input_click([1,3]));
        case 'Asymmetric Rotating Gaussian'
            LOCALIZE_ASYMMETRIC_ROTATING_GAUSSIAN(molecule_name,input_click([1,3]));
        case 'Two Gaussians'
            minRowROI = max([1,round(input_click(3))-fitting_region]);
            minColROI = max([1,round(input_click(1))-fitting_region]);
            maxRowROI = min([size(red_channel_image_scaled,1),round(input_click(3))+fitting_region]);
            maxColROI = min([size(red_channel_image_scaled,2),round(input_click(1))+fitting_region]);
%             [rows,cols] = ndgrid(minRowROI:maxRowROI,minColROI:maxColROI);
            ROI = red_channel_image_scaled(minRowROI:maxRowROI,minColROI:maxColROI);
            [rowMax, colMax] = ind2sub(size(ROI),find(imregionalmax(ROI))); % value
            sortedPeaks = sortrows([rowMax,colMax,ROI(sub2ind(size(ROI),rowMax,colMax))],-3);%             keydown = waitforbuttonpress;
%             if (keydown == 0)
%                 input_click2 = get(fp.current_frame_axes,'Currentpoint'); % click on second point
%                 LOCALIZE_2X_GAUSSIANS(molecule_name,input_click([1,3]),input_click2([1,3]));
%             else
            input_click = [sortedPeaks(1,1),sortedPeaks(1,2)]+[minRowROI, minColROI]-1;
            input_click2 = [sortedPeaks(2,1),sortedPeaks(2,2)]+[minRowROI, minColROI]-1; % if the mouse button isn't pressed, just use the original point
            LOCALIZE_2X_GAUSSIANS(molecule_name,input_click([2,1]),input_click2([2,1]));
%             end
    end   
    UPDATE_IMAGE;
    DISPLAY_FITTING_HISTOGRAM;
end

function UPDATE_IMAGE(~,~)
    red_channel_image=single(imread([FILEDATA.red_channel_filelocation FILEDATA.red_channel_filename],'Index',frame,'Info',FILEDATA.red_channel_fileinfo));
    red_channel_image_scaled = mat2gray(red_channel_image,double([(1+scaling(1,1))*min(red_channel_image(:)), max(red_channel_image(:))*scaling(1,2)]));
    set(fp.current_frame,'CData',cat(3,red_channel_image_scaled,green_channel_image_scaled,blue_channel_image_scaled));
    set(fp.data_frame_number,'string',frame);
    set(fp.plot_positions,'xdata',DATA.(genvarname(current_molecule_name)).col,'ydata',DATA.(genvarname(current_molecule_name)).row);
    set(fp.plot_recent_position,'xdata',DATA.(genvarname(current_molecule_name)).col(lower_frame_bin:upper_frame_bin),'ydata',DATA.(genvarname(current_molecule_name)).row(lower_frame_bin:upper_frame_bin));
    position_data = NaN(size(list_of_molecules,1),2); n = 0;
    for molecule = 1:numel(list_of_molecules)
        molecule_name = char(list_of_molecules(molecule));
        n = n+1;
        position_data(n,1) = DATA.(genvarname(molecule_name)).col(frame);
        position_data(n,2) = DATA.(genvarname(molecule_name)).row(frame);
    end
    set(fp.scatter_positions,'xdata',position_data(:,1),'ydata',position_data(:,2));

    switch char(list_of_tracking_modes(get(tracking_mode_popup,'value')))
        case 'Symmetric Gaussian'
            current_row = num2str(DATA.(genvarname(current_molecule_name)).row(frame),'%.2f');
            current_col = num2str(DATA.(genvarname(current_molecule_name)).col(frame),'%.2f');
            current_sigma = num2str(DATA.(genvarname(current_molecule_name)).sigma(frame),'%.2f');
            current_amplitude = num2str(DATA.(genvarname(current_molecule_name)).amplitude(frame),'%.1f');
        case 'Restricted Symmetric Gaussian'
            current_row = num2str(DATA.(genvarname(current_molecule_name)).row(frame),'%.2f');
            current_col = num2str(DATA.(genvarname(current_molecule_name)).col(frame),'%.2f');
            current_sigma = num2str(DATA.(genvarname(current_molecule_name)).sigma(frame),'%.2f');
            current_amplitude = num2str(DATA.(genvarname(current_molecule_name)).amplitude(frame),'%.1f');
        case 'Asymmetric Gaussian'
            current_row = num2str(DATA.(genvarname(current_molecule_name)).row(frame),'%.2f');
            current_col = num2str(DATA.(genvarname(current_molecule_name)).col(frame),'%.2f');
            current_sigma = {num2str(DATA.(genvarname(current_molecule_name)).sigma_row(frame),'%.2f');num2str(DATA.(genvarname(current_molecule_name)).sigma_col(frame),'%.2f')};
            current_amplitude = num2str(DATA.(genvarname(current_molecule_name)).amplitude(frame),'%.1f');
        case 'Asymmetric Rotating Gaussian'
            current_row = num2str(DATA.(genvarname(current_molecule_name)).row(frame),'%.2f');
            current_col = num2str(DATA.(genvarname(current_molecule_name)).col(frame),'%.2f');
            current_sigma = {num2str(DATA.(genvarname(current_molecule_name)).sigma_row(frame),'%.2f');num2str(DATA.(genvarname(current_molecule_name)).sigma_col(frame),'%.2f')};
            current_amplitude = num2str(DATA.(genvarname(current_molecule_name)).amplitude(frame),'%.1f');
        case 'Two Gaussians'
            current_row = {num2str(DATA.(genvarname(current_molecule_name)).row_1(frame),'%.2f');num2str(DATA.(genvarname(current_molecule_name)).row_2(frame),'%.2f')};
            current_col = {num2str(DATA.(genvarname(current_molecule_name)).col_1(frame),'%.2f');num2str(DATA.(genvarname(current_molecule_name)).col_2(frame),'%.2f')};
            current_sigma = {num2str(DATA.(genvarname(current_molecule_name)).sigma_1(frame),'%.2f');num2str(DATA.(genvarname(current_molecule_name)).sigma_2(frame),'%.2f')};
            current_amplitude = {num2str(DATA.(genvarname(current_molecule_name)).amplitude_1(frame),'%.2f');num2str(DATA.(genvarname(current_molecule_name)).amplitude_2(frame),'%.2f')};
    end
    set(fp.data_row_position,'string',current_row);
    set(fp.data_col_position,'string',current_col);
    set(fp.data_sigma,'string',current_sigma);
    set(fp.data_amplitude,'string',current_amplitude)
    set(fp.data_background,'string',num2str(DATA.(genvarname(current_molecule_name)).background(frame),'%.1f'))

    if ~isnan(DATA.(genvarname(current_molecule_name)).row(frame))
        set(dp.data_plot,'xdata',0,'ydata',1)
        DISPLAY_FITTING_HISTOGRAM;
    end
end

function DISPLAY_FITTING_HISTOGRAM(~,~)
    current_variable = char(cellstr(list_of_variables_for_hist(get(dp.data_selection,'value'))));
    ydata = DATA.(genvarname(current_molecule_name)).(genvarname(current_variable));
    xdata = DATA.(genvarname(current_molecule_name)).frame;
    ydata_min = min(ydata)*0.99-.01;
    ydata_max = max(ydata)*1.01+.01;
    edges = linspace(ydata_min,ydata_max,20);
    min_frame = DATA.(genvarname(current_molecule_name)).min_frame_used;
    max_frame = DATA.(genvarname(current_molecule_name)).max_frame_used;
    set(dp.data_min_max_frames,'string',[num2str(min_frame) ' / ' num2str(max_frame)]);
    set(dp.data_plot,'xdata',xdata,'ydata',ydata)
    set(dp.data_plot_axes,'YLim',[ydata_min,ydata_max])
    hdata = histcounts(ydata,edges);
    set(dp.data_hist,'ydata',hdata,'xdata',edges);
end

function PREPARE_TO_LOCALIZE
    if get(cp.toggle_tracking,'Value')
        list_of_selected_molecules = get(Selected_molecules,'string');
    else
        list_of_selected_molecules = {};
    end
    
    switch char(list_of_tracking_modes(get(tracking_mode_popup,'value')))
        case 'Symmetric Gaussian'
            for molecule = 1:numel(list_of_selected_molecules)
                molecule_name = char(list_of_selected_molecules(molecule));
                 if DATA.(genvarname(molecule_name)).frame_locked(frame) == 0 && ...
                   DATA.(genvarname(molecule_name)).frame_enabled(frame) == 1 && ...
                   ge(frame,DATA.(genvarname(molecule_name)).min_frame_allowed) == 1 && ...
                   le(frame,DATA.(genvarname(molecule_name)).max_frame_allowed) == 1
                    input_click(1) = mean(DATA.(genvarname(molecule_name)).col(lower_frame_bin:upper_frame_bin),'omitnan');
                    input_click(2) = mean(DATA.(genvarname(molecule_name)).row(lower_frame_bin:upper_frame_bin),'omitnan');
                    LOCALIZE_SYMMETRIC_GAUSSIAN(molecule_name,input_click);
                end                 
            end       
        case 'Restricted Symmetric Gaussian'
            for molecule = 1:numel(list_of_selected_molecules)
                molecule_name = char(list_of_selected_molecules(molecule));
                 if DATA.(genvarname(molecule_name)).frame_locked(frame) == 0 && ...
                   DATA.(genvarname(molecule_name)).frame_enabled(frame) == 1 && ...
                   ge(frame,DATA.(genvarname(molecule_name)).min_frame_allowed) == 1 && ...
                   le(frame,DATA.(genvarname(molecule_name)).max_frame_allowed) == 1
                    input_click(1) = mean(DATA.(genvarname(molecule_name)).col(frame),'omitnan');
                    input_click(2) = mean(DATA.(genvarname(molecule_name)).row(frame),'omitnan');
                    LOCALIZE_RESTRICTED_GAUSSIAN(molecule_name,input_click);
                end                 
            end  
        case 'Asymmetric Gaussian'
            for molecule = 1:numel(list_of_selected_molecules)
                molecule_name = char(list_of_selected_molecules(molecule));
                 if DATA.(genvarname(molecule_name)).frame_locked(frame) == 0 && ...
                   DATA.(genvarname(molecule_name)).frame_enabled(frame) == 1 && ...
                   ge(frame,DATA.(genvarname(molecule_name)).min_frame_allowed) == 1 && ...
                   le(frame,DATA.(genvarname(molecule_name)).max_frame_allowed) == 1
                    input_click(1) = mean(DATA.(genvarname(molecule_name)).col(lower_frame_bin:upper_frame_bin),'omitnan');
                    input_click(2) = mean(DATA.(genvarname(molecule_name)).row(lower_frame_bin:upper_frame_bin),'omitnan');
                    LOCALIZE_ASYMMETRIC_GAUSSIAN(molecule_name,input_click);
                end                 
            end 
        case 'Asymmetric Rotating Gaussian'
            for molecule = 1:numel(list_of_selected_molecules)
                molecule_name = char(list_of_selected_molecules(molecule));
                 if DATA.(genvarname(molecule_name)).frame_locked(frame) == 0 && ...
                   DATA.(genvarname(molecule_name)).frame_enabled(frame) == 1 && ...
                   ge(frame,DATA.(genvarname(molecule_name)).min_frame_allowed) == 1 && ...
                   le(frame,DATA.(genvarname(molecule_name)).max_frame_allowed) == 1
                    input_click(1) = mean(DATA.(genvarname(molecule_name)).col(lower_frame_bin:upper_frame_bin),'omitnan');
                    input_click(2) = mean(DATA.(genvarname(molecule_name)).row(lower_frame_bin:upper_frame_bin),'omitnan');
                    LOCALIZE_ASYMMETRIC_ROTATING_GAUSSIAN(molecule_name,input_click);
                end                 
            end 
        case 'Two Gaussians'
            for molecule = 1:numel(list_of_selected_molecules)
                molecule_name = char(list_of_selected_molecules(molecule));
                if DATA.(genvarname(molecule_name)).frame_locked(frame) == 0 && ...
                   DATA.(genvarname(molecule_name)).frame_enabled(frame) == 1 && ...
                   ge(frame,DATA.(genvarname(molecule_name)).min_frame_allowed) == 1 && ...
                   le(frame,DATA.(genvarname(molecule_name)).max_frame_allowed) == 1
                    input_click(1) = mean(DATA.(genvarname(molecule_name)).col_1(lower_frame_bin:upper_frame_bin),'omitnan');
                    input_click(2) = mean(DATA.(genvarname(molecule_name)).row_1(lower_frame_bin:upper_frame_bin),'omitnan');
                    input_click2(1) = mean(DATA.(genvarname(molecule_name)).col_2(lower_frame_bin:upper_frame_bin),'omitnan');
                    input_click2(2) = mean(DATA.(genvarname(molecule_name)).row_2(lower_frame_bin:upper_frame_bin),'omitnan');
                    LOCALIZE_2X_GAUSSIANS(molecule_name,input_click,input_click2);
                end  
                
            end
    end
end
%% LOCALIZATION CODE
function LOCALIZE_SYMMETRIC_GAUSSIAN(molecule_name,input_click)
    input_pixels = round(input_click);
    pixel_limits = {[round(input_pixels(2))-fitting_region,input_pixels(2)+fitting_region],[input_pixels(1)-fitting_region,round(input_pixels(1))+fitting_region]};
    % Determine if it's okay
    Image_Region = double(imread([FILEDATA.red_channel_filelocation FILEDATA.red_channel_filename],'Index',frame,'Info',FILEDATA.red_channel_fileinfo,'PixelRegion',pixel_limits));
    % Prepare Guesses
    BG_guess = (mean(border_region(:).*Image_Region(:),'omitnan')-background_limits(1))/background_range;
    AMP_guess = (max(Image_Region(:))-mean(border_region(:).*Image_Region(:),'omitnan')-amplitude_limits(1))/amplitude_range;
    % Prepare lower and upper bounds
    [Localized, resnorm, residual, exitflag] = lsqnonlin(@SYMMETRIC_GAUSSIAN_FIT, [.5, .5, AMP_guess, .5, BG_guess], [0 0 0 0 0], [1, 1, 1, 1, 1], fitting_options);
    fit_localization = (Localized(5)*background_range+background_limits(1)) + (Localized(3)*amplitude_range+amplitude_limits(1))*exp((-(regional_indices_row-((Localized(1)-.5)*(fitting_region+1))).^2-(regional_indices_col-((Localized(2)-.5)*(fitting_region+1))).^2)/(2*(Localized(4)*sigma_range+sigma_limits(1))^2));
    if strcmp(molecule_name,current_molecule_name)==1
        fit_region_image = repmat(mat2gray(Image_Region),[1,1,3]);
        fit_localization_image = repmat(mat2gray(fit_localization),[1,1,3]);
        residual_in_gray = mat2gray(reshape(residual,[fitting_region*2+1,fitting_region*2+1]));
        set(lp.raw_data,'CData',fit_region_image);
        set(lp.fit_data,'CData',fit_localization_image);
        set(lp.residual_data,'CData',cat(3,residual_in_gray,residual_in_gray,residual_in_gray));
    end
    SYMMETRIC_LOCALIZTION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag);
end
function [Delta, Jacobian] = SYMMETRIC_GAUSSIAN_FIT(guess)
    Row = (guess(1)-.5)*fitting_region;
    Col = (guess(2)-.5)*fitting_region;
    Amplitude = guess(3)*amplitude_range+amplitude_limits(1);
    Sigma = guess(4)*sigma_range+sigma_limits(1);
    Background = guess(5)*background_range+background_limits(1);
    Guess_Image = Background + Amplitude*exp((-(regional_indices_row-Row).^2-(regional_indices_col-Col).^2)/(2*Sigma^2));
    Delta = - Image_Region(:) + Guess_Image(:);
    if nargout > 1
        Jacobian= zeros(numel(Image_Region),5);
        Jacobian(:,1) = Amplitude*exp((-(regional_indices_row(:)-Row).^2-(regional_indices_col(:)-Col).^2)/(2*Sigma^2)) .* (regional_indices_row(:)-Row)/(Sigma^2);
        Jacobian(:,2) = Amplitude*exp((-(regional_indices_row(:)-Row).^2-(regional_indices_col(:)-Col).^2)/(2*Sigma^2)) .* (regional_indices_col(:)-Col)/(Sigma^2);
        Jacobian(:,3) =           exp((-(regional_indices_row(:)-Row).^2-(regional_indices_col(:)-Col).^2)/(2*Sigma^2));
        Jacobian(:,4) = Amplitude*exp((-(regional_indices_row(:)-Row).^2-(regional_indices_col(:)-Col).^2)/(2*Sigma^2)) .* ((regional_indices_row(:)-Row).^2+(regional_indices_col(:)-Col).^2) / Sigma^3;
        Jacobian(:,5) = ones(numel(Image_Region),1);
    end
end

function SYMMETRIC_LOCALIZTION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag)
    if Localized(4)*sigma_range+sigma_limits(1) < .0
        if DATA.(genvarname(molecule_name)).localization_warning == 0
           % IF BAD, BUT NO WARNING ADD A WARNING
           DATA.(genvarname(molecule_name)).localization_warning = 1;
        else
            % IF BAD, WITH WARNING, REMOVE TRACK
            Unselected = cellstr(get(Unselected_molecules,'String'));
            Selected = cellstr(get(Selected_molecules, 'string'));   
            new_clicked = find(strcmp(molecule_name,Selected)==1, 1);
            Current_selection = Selected(new_clicked);
            Unselected = [Unselected; Current_selection];
            Selected(new_clicked,:) = [];
            Unselected(any(cellfun(@isempty,Unselected),2),:) = [];
            Selected(any(cellfun(@isempty,Selected),2),:) = [];
            set(Unselected_molecules,'String', sort(Unselected),'Value',1)
            set(Selected_molecules,'String', sort(Selected),'Value',1)
            if isempty(Selected) == 1
                % IF ONLY TRACK, STOP RUNNING
                toggle_autorun = [0 0];
            end
        end
    else
        % IF GOOD, REMOVE A WARNING
        DATA.(genvarname(molecule_name)).localization_warning = 0;
        DATA.(genvarname(molecule_name)).signal_above_BG(frame) = sum(Image_Region(:)-Localized(5)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).fitted_signal(frame) = sum(fit_localization(:)-Localized(5)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).residual(frame) = resnorm;
        DATA.(genvarname(molecule_name)).row(frame) = input_pixels(2)+(Localized(1)-.5)*fitting_region;
        DATA.(genvarname(molecule_name)).col(frame) = input_pixels(1)+(Localized(2)-.5)*fitting_region;
        DATA.(genvarname(molecule_name)).amplitude(frame) = Localized(3)*amplitude_range+amplitude_limits(1);
        DATA.(genvarname(molecule_name)).background(frame) = Localized(5)*background_range+background_limits(1);
        DATA.(genvarname(molecule_name)).sigma(frame) = Localized(4)*sigma_range+sigma_limits(1);
        DATA.(genvarname(molecule_name)).frame_enabled(lower_frame_bin:upper_frame_bin) = 1;
        DATA.(genvarname(molecule_name)).localization_info(frame) = exitflag;
        DATA.(genvarname(molecule_name)).min_frame_used = min(frame,(DATA.(genvarname(molecule_name)).min_frame_used));
        DATA.(genvarname(molecule_name)).max_frame_used = max(frame,(DATA.(genvarname(molecule_name)).max_frame_used));
    end
end
function LOCALIZE_RESTRICTED_GAUSSIAN(molecule_name,input_click) %#ok<INUSD>
%     input_pixels = round(input_click);
    input_pixels(1) = DATA.(genvarname(molecule_name)).col(frame);
    input_pixels(2) = DATA.(genvarname(molecule_name)).row(frame);
    pixel_limits = {[round(input_pixels(2))-fitting_region,round(input_pixels(2))+fitting_region],[round(input_pixels(1))-fitting_region,round(input_pixels(1))+fitting_region]};
    % Determine if it's okay
    Image_Region = double(imread([FILEDATA.red_channel_filelocation FILEDATA.red_channel_filename],'Index',frame,'Info',FILEDATA.red_channel_fileinfo,'PixelRegion',pixel_limits));
    % Prepare Guesses
    %%%%%%%%%%%%%% 
    amplitude_limits = round([.5*DATA.(genvarname(molecule_name)).amplitude(frame),2*DATA.(genvarname(molecule_name)).amplitude(frame)]);
    set(cp.amplitude_min,'string',num2str(amplitude_limits(1)));
    set(cp.amplitude_max,'string',num2str(amplitude_limits(2)));
    amplitude_range = diff(amplitude_limits);
    AMP_guess = (DATA.(genvarname(molecule_name)).amplitude(frame)-amplitude_limits(1))/amplitude_range;
    
    background_limits = round([.995*DATA.(genvarname(molecule_name)).background(frame),1.005*DATA.(genvarname(molecule_name)).background(frame)]);
    set(cp.background_min,'string',num2str(background_limits(1)));
    set(cp.background_max,'string',num2str(background_limits(2)));
    background_range = diff(background_limits);
    BG_guess = (DATA.(genvarname(molecule_name)).background(frame)-background_limits(1))/background_range;

    %%%%%%%%%%%%%
    % Prepare lower and upper bounds
    [Localized, resnorm, residual, exitflag] = lsqnonlin(@RESTRICTED_GAUSSIAN_FIT, [(input_pixels(2)-round(input_pixels(2)))/(2*fitting_region+1)+.5, (input_pixels(1)-round(input_pixels(1)))/(2*fitting_region+1)+.5, AMP_guess, .5, BG_guess], [.375 .375 0 0 0], [.625, .625, 1, 1, 1], fitting_options);
    fit_localization = (Localized(5)*background_range+background_limits(1)) + (Localized(3)*amplitude_range+amplitude_limits(1))*exp((-(regional_indices_row-((Localized(1)-.5)*(fitting_region+1))).^2-(regional_indices_col-((Localized(2)-.5)*(fitting_region+1))).^2)/(2*(Localized(4)*sigma_range+sigma_limits(1))^2));
    if strcmp(molecule_name,current_molecule_name)==1
        fit_region_image = repmat(mat2gray(Image_Region),[1,1,3]);
        fit_localization_image = repmat(mat2gray(fit_localization),[1,1,3]);
        residual_in_gray = mat2gray(reshape(residual,[fitting_region*2+1,fitting_region*2+1]));
        set(lp.raw_data,'CData',fit_region_image);
        set(lp.fit_data,'CData',fit_localization_image);
        set(lp.residual_data,'CData',cat(3,residual_in_gray,residual_in_gray,residual_in_gray));
    end
    RESTRICTED_LOCALIZATION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag);
end

function [Delta, Jacobian] = RESTRICTED_GAUSSIAN_FIT(guess)
    Row = (guess(1)-.5)*fitting_region;
    Col = (guess(2)-.5)*fitting_region;
    Amplitude = guess(3)*amplitude_range+amplitude_limits(1);
    Sigma = guess(4)*sigma_range+sigma_limits(1);
    Background = guess(5)*background_range+background_limits(1);
    Guess_Image = Background + Amplitude*exp((-(regional_indices_row-Row).^2-(regional_indices_col-Col).^2)/(2*Sigma^2));
    Delta = - Image_Region(:) + Guess_Image(:);
    if nargout > 1
        Jacobian= zeros(numel(Image_Region),5);
        Jacobian(:,1) = Amplitude*exp((-(regional_indices_row(:)-Row).^2-(regional_indices_col(:)-Col).^2)/(2*Sigma^2)) .* (regional_indices_row(:)-Row)/(Sigma^2);
        Jacobian(:,2) = Amplitude*exp((-(regional_indices_row(:)-Row).^2-(regional_indices_col(:)-Col).^2)/(2*Sigma^2)) .* (regional_indices_col(:)-Col)/(Sigma^2);
        Jacobian(:,3) =           exp((-(regional_indices_row(:)-Row).^2-(regional_indices_col(:)-Col).^2)/(2*Sigma^2));
        Jacobian(:,4) = Amplitude*exp((-(regional_indices_row(:)-Row).^2-(regional_indices_col(:)-Col).^2)/(2*Sigma^2)) .* ((regional_indices_row(:)-Row).^2+(regional_indices_col(:)-Col).^2) / Sigma^3;
    end
end

function RESTRICTED_LOCALIZATION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag)
    if Localized(4)*sigma_range+sigma_limits(1) < .0
        if DATA.(genvarname(molecule_name)).localization_warning == 0
           % IF BAD, BUT NO WARNING ADD A WARNING
           DATA.(genvarname(molecule_name)).localization_warning = 1;
        else
            % IF BAD, WITH WARNING, REMOVE TRACK
            Unselected = cellstr(get(Unselected_molecules,'String'));
            Selected = cellstr(get(Selected_molecules, 'string'));   
            new_clicked = find(strcmp(molecule_name,Selected)==1, 1);
            Current_selection = Selected(new_clicked);
            Unselected = [Unselected; Current_selection];
            Selected(new_clicked,:) = [];
            Unselected(any(cellfun(@isempty,Unselected),2),:) = [];
            Selected(any(cellfun(@isempty,Selected),2),:) = [];
            set(Unselected_molecules,'String', sort(Unselected),'Value',1)
            set(Selected_molecules,'String', sort(Selected),'Value',1)
            if isempty(Selected) == 1
                % IF ONLY TRACK, STOP RUNNING
                toggle_autorun = [0 0];
            end
        end
    else
        % IF GOOD, REMOVE A WARNING
        DATA.(genvarname(molecule_name)).localization_warning = 0;
        DATA.(genvarname(molecule_name)).signal_above_BG(frame) = sum(Image_Region(:)-Localized(5)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).fitted_signal(frame) = sum(fit_localization(:)-Localized(5)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).residual(frame) = resnorm;
        DATA.(genvarname(molecule_name)).row(frame) = input_pixels(2)+(Localized(1)-.5)*fitting_region;
        DATA.(genvarname(molecule_name)).col(frame) = input_pixels(1)+(Localized(2)-.5)*fitting_region;
        DATA.(genvarname(molecule_name)).amplitude(frame) = Localized(3)*amplitude_range+amplitude_limits(1);
        DATA.(genvarname(molecule_name)).background(frame) = Localized(5)*background_range+background_limits(1);
        DATA.(genvarname(molecule_name)).sigma(frame) = Localized(4)*sigma_range+sigma_limits(1);
        DATA.(genvarname(molecule_name)).frame_enabled(lower_frame_bin:upper_frame_bin) = 1;
        DATA.(genvarname(molecule_name)).localization_info(frame) = exitflag;
        DATA.(genvarname(molecule_name)).min_frame_used = min(frame,(DATA.(genvarname(molecule_name)).min_frame_used));
        DATA.(genvarname(molecule_name)).max_frame_used = max(frame,(DATA.(genvarname(molecule_name)).max_frame_used));
    end
end

function LOCALIZE_ASYMMETRIC_GAUSSIAN(molecule_name,input_click)
    input_pixels = round(input_click);
    pixel_limits = {[input_pixels(2)-fitting_region,input_pixels(2)+fitting_region],[input_pixels(1)-fitting_region,input_pixels(1)+fitting_region]};
    % Determine if it's okay
    Image_Region = double(imread([FILEDATA.red_channel_filelocation FILEDATA.red_channel_filename],'Index',frame,'Info',FILEDATA.red_channel_fileinfo,'PixelRegion',pixel_limits));
    % Prepare Guesses
    BG_guess = (mean(border_region(:).*Image_Region(:),'omitnan')-background_limits(1))/background_range;
    AMP_guess = (max(Image_Region(:))-mean(border_region(:).*Image_Region(:),'omitnan')-amplitude_limits(1))/amplitude_range;
    % Prepare lower and upper bounds
    tic
    [Localized, resnorm, residual, exitflag] = lsqnonlin(@ASYMMETRIC_GAUSSIAN_FIT, [.5, .5, AMP_guess, .5, .5, BG_guess], [0 0 0 0 0 0], [1, 1, 1, 1, 1, 1], fitting_options);
    toc
    fit_localization = (Localized(6)*background_range+background_limits(1)) + (Localized(3)*amplitude_range+amplitude_limits(1))*exp((-(regional_indices_row-((Localized(1)-.5)*(fitting_region+1))).^2)/(2*(Localized(4)*sigma_range+sigma_limits(1))^2)-((regional_indices_col-((Localized(2)-.5)*(fitting_region+1))).^2)/(2*(Localized(5)*sigma_range+sigma_limits(1))^2));
    if strcmp(molecule_name,current_molecule_name)==1
        fit_region_image = repmat(mat2gray(Image_Region),[1,1,3]);
        fit_localization_image = repmat(mat2gray(fit_localization),[1,1,3]);
        residual_in_gray = mat2gray(reshape(residual,[fitting_region*2+1,fitting_region*2+1]));
        set(lp.raw_data,'CData',fit_region_image);
        set(lp.fit_data,'CData',fit_localization_image);
        set(lp.residual_data,'CData',cat(3,residual_in_gray,residual_in_gray,residual_in_gray));
    end
    ASYMMETRIC_LOCALIZATION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag);
end

function [Delta] = ASYMMETRIC_GAUSSIAN_FIT(guess)
    Row = (guess(1)-.5)*fitting_region;
    Col = (guess(2)-.5)*fitting_region;
    Amplitude = guess(3)*amplitude_range+amplitude_limits(1);
    Sigma_row = guess(4)*sigma_range+sigma_limits(1);
    Sigma_col = guess(5)*sigma_range+sigma_limits(1);
    Background = guess(6)*background_range+background_limits(1);
    Guess_Image = Background + Amplitude*exp(-((regional_indices_row-Row).^2)/(2*Sigma_row^2)-((regional_indices_col-Col).^2)/(2*Sigma_col^2));
    Delta = - Image_Region(:) + Guess_Image(:);
end

function ASYMMETRIC_LOCALIZATION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag)
    if Localized(4)*sigma_range+sigma_limits(1) < .5 || Localized(5)*sigma_range+sigma_limits(1) < .5
        if DATA.(genvarname(molecule_name)).localization_warning == 0
           % IF BAD, BUT NO WARNING ADD A WARNING
           DATA.(genvarname(molecule_name)).localization_warning = 1;
        else
            % IF BAD, WITH WARNING, REMOVE TRACK
            Unselected = cellstr(get(Unselected_molecules,'String'));
            Selected = cellstr(get(Selected_molecules, 'string'));   
            new_clicked = find(strcmp(molecule_name,Selected)==1, 1);
            Current_selection = Selected(new_clicked);
            Unselected = [Unselected; Current_selection];
            Selected(new_clicked,:) = [];
            Unselected(any(cellfun(@isempty,Unselected),2),:) = [];
            Selected(any(cellfun(@isempty,Selected),2),:) = [];
            set(Unselected_molecules,'String', sort(Unselected),'Value',1)
            set(Selected_molecules,'String', sort(Selected),'Value',1)
            if isempty(Selected) == 1
                % IF ONLY TRACK, STOP RUNNING
                toggle_autorun = [0 0];
            end
        end
    else
        % IF GOOD, REMOVE A WARNING
        DATA.(genvarname(molecule_name)).localization_warning = 0;
        DATA.(genvarname(molecule_name)).signal_above_BG(frame) = sum(Image_Region(:)-Localized(6)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).fitted_signal(frame) = sum(fit_localization(:)-Localized(6)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).residual(frame) = resnorm;
        DATA.(genvarname(molecule_name)).row(frame) = input_pixels(2)+(Localized(1)-.5)*fitting_region;
        DATA.(genvarname(molecule_name)).col(frame) = input_pixels(1)+(Localized(2)-.5)*fitting_region;
        DATA.(genvarname(molecule_name)).amplitude(frame) = Localized(3)*amplitude_range+amplitude_limits(1);
        DATA.(genvarname(molecule_name)).background(frame) = Localized(6)*background_range+background_limits(1);
        DATA.(genvarname(molecule_name)).sigma_row(frame) = Localized(4)*sigma_range+sigma_limits(1);
        DATA.(genvarname(molecule_name)).sigma_col(frame) = Localized(5)*sigma_range+sigma_limits(1);
        DATA.(genvarname(molecule_name)).frame_enabled(lower_frame_bin:upper_frame_bin) = 1;
        DATA.(genvarname(molecule_name)).localization_info(frame) = exitflag;
        DATA.(genvarname(molecule_name)).min_frame_used = min(frame,(DATA.(genvarname(molecule_name)).min_frame_used));
        DATA.(genvarname(molecule_name)).max_frame_used = max(frame,(DATA.(genvarname(molecule_name)).max_frame_used));
    end
end

function LOCALIZE_ASYMMETRIC_ROTATING_GAUSSIAN(molecule_name,input_click)
    input_pixels = round(input_click);
    pixel_limits = {[input_pixels(2)-fitting_region,input_pixels(2)+fitting_region],[input_pixels(1)-fitting_region,input_pixels(1)+fitting_region]};
    % Determine if it's okay
    Image_Region = double(imread([FILEDATA.red_channel_filelocation FILEDATA.red_channel_filename],'Index',frame,'Info',FILEDATA.red_channel_fileinfo,'PixelRegion',pixel_limits));
    % Prepare Guesses
    BG_guess = (mean(border_region(:).*Image_Region(:),'omitnan')-background_limits(1))/background_range;
    AMP_guess = (max(Image_Region(:))-mean(border_region(:).*Image_Region(:),'omitnan')-amplitude_limits(1))/amplitude_range;
    % Prepare lower and upper bounds
    [Localized, resnorm, residual, exitflag] = lsqnonlin(@ASYMMETRIC_ROTATING_GAUSSIAN_FIT, [.5, .5, AMP_guess, .5, .5, .5, BG_guess], [0 0 0 0 0 0 0], [1, 1, 1, 1, 1, 1, 1], fitting_options);
    
    theta_limits = [-pi/4 pi/4];
    theta_range =  diff(theta_limits);
    Row = (Localized(1)-.5)*(fitting_region+1);
    Col = (Localized(2)-.5)*(fitting_region+1);
    Sigma_row = Localized(4)*sigma_range+sigma_limits(1);
    Sigma_col = Localized(5)*sigma_range+sigma_limits(1);
    Theta = Localized(6)*theta_range+theta_limits(1);

    a = ( cos(Theta)^2 / (2*Sigma_row^2)) + (sin(Theta)^2 / (2*Sigma_col^2));
    b = (-sin(2*Theta) / (4*Sigma_row^2)) + (sin(2*Theta) / (4*Sigma_col^2));
    c = ( sin(Theta)^2 / (2*Sigma_row^2)) + (cos(Theta)^2 / (2*Sigma_col^2));
    fit_localization = (Localized(7)*background_range+background_limits(1)) +...
                       (Localized(3)*amplitude_range+amplitude_limits(1))* ...
                        exp(-(a*(regional_indices_row-Row).^2 + 2*b*(regional_indices_row-Row).*(regional_indices_col-Col)+c*(regional_indices_col-Col).^2));
%     figure(2)
    
%     fit_localization = (Localized(7)*background_range+background_limits(1)) + (Localized(3)*amplitude_range+amplitude_limits(1))*exp((-(regional_indices_row-((Localized(1)-.5)*(fitting_region+1))).^2)/(2*(Localized(4)*sigma_range+sigma_limits(1))^2)-((regional_indices_col-((Localized(2)-.5)*(fitting_region+1))).^2)/(2*(Localized(5)*sigma_range+sigma_limits(1))^2));
    if strcmp(molecule_name,current_molecule_name)==1
        fit_region_image = repmat(mat2gray(Image_Region),[1,1,3]);
        fit_localization_image = repmat(mat2gray(fit_localization),[1,1,3]);
        residual_in_gray = mat2gray(reshape(residual,[fitting_region*2+1,fitting_region*2+1]));
        set(lp.raw_data,'CData',fit_region_image);
        set(lp.fit_data,'CData',fit_localization_image);
        set(lp.residual_data,'CData',cat(3,residual_in_gray,residual_in_gray,residual_in_gray));
    end
    ASYMMETRIC_ROTATING_LOCALIZATION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag);
end

function [Delta] = ASYMMETRIC_ROTATING_GAUSSIAN_FIT(guess)
    theta_limits = [-pi/4 pi/4];
    theta_range =  diff(theta_limits);
    Row = (guess(1)-.5)*fitting_region;
    Col = (guess(2)-.5)*fitting_region;
    Amplitude = guess(3)*amplitude_range+amplitude_limits(1);
    Sigma_row = guess(4)*sigma_range+sigma_limits(1);
    Sigma_col = guess(5)*sigma_range+sigma_limits(1);
    Theta = guess(6)*theta_range+theta_limits(1);
    Background = guess(7)*background_range+background_limits(1);

    a = ( cos(Theta)^2 / (2*Sigma_row^2)) + (sin(Theta)^2 / (2*Sigma_col^2));
    b = (-sin(2*Theta) / (4*Sigma_row^2)) + (sin(2*Theta) / (4*Sigma_col^2));
    c = ( sin(Theta)^2 / (2*Sigma_row^2)) + (cos(Theta)^2 / (2*Sigma_col^2));
    Guess_Image = Background + Amplitude*exp(-(a*(regional_indices_row-Row).^2 + 2*b*(regional_indices_row-Row).*(regional_indices_col-Col)+c*(regional_indices_col-Col).^2));
%     figure(2)
%     imagesc(Guess_Image)
    Delta = - Image_Region(:) + Guess_Image(:);
end

function ASYMMETRIC_ROTATING_LOCALIZATION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag)
    if Localized(4)*sigma_range+sigma_limits(1) < .5 || Localized(5)*sigma_range+sigma_limits(1) < .5
        if DATA.(genvarname(molecule_name)).localization_warning == 0
           % IF BAD, BUT NO WARNING ADD A WARNING
           DATA.(genvarname(molecule_name)).localization_warning = 1;
        else
            % IF BAD, WITH WARNING, REMOVE TRACK
            Unselected = cellstr(get(Unselected_molecules,'String'));
            Selected = cellstr(get(Selected_molecules, 'string'));   
            new_clicked = find(strcmp(molecule_name,Selected)==1, 1);
            Current_selection = Selected(new_clicked);
            Unselected = [Unselected; Current_selection];
            Selected(new_clicked,:) = [];
            Unselected(any(cellfun(@isempty,Unselected),2),:) = [];
            Selected(any(cellfun(@isempty,Selected),2),:) = [];
            set(Unselected_molecules,'String', sort(Unselected),'Value',1)
            set(Selected_molecules,'String', sort(Selected),'Value',1)
            if isempty(Selected) == 1
                % IF ONLY TRACK, STOP RUNNING
                toggle_autorun = [0 0];
            end
        end
    else
        % IF GOOD, REMOVE A WARNING
        theta_limits = [-pi/4 pi/4];
        theta_range =  diff(theta_limits);
        DATA.(genvarname(molecule_name)).localization_warning = 0;
        DATA.(genvarname(molecule_name)).signal_above_BG(frame) = sum(Image_Region(:)-Localized(6)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).fitted_signal(frame) = sum(fit_localization(:)-Localized(6)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).residual(frame) = resnorm;
        DATA.(genvarname(molecule_name)).row(frame) = input_pixels(2)+(Localized(1)-.5)*fitting_region;
        DATA.(genvarname(molecule_name)).col(frame) = input_pixels(1)+(Localized(2)-.5)*fitting_region;
        DATA.(genvarname(molecule_name)).amplitude(frame) = Localized(3)*amplitude_range+amplitude_limits(1);
        DATA.(genvarname(molecule_name)).background(frame) = Localized(7)*background_range+background_limits(1);
        DATA.(genvarname(molecule_name)).sigma_row(frame) = Localized(4)*sigma_range+sigma_limits(1);
        DATA.(genvarname(molecule_name)).sigma_col(frame) = Localized(5)*sigma_range+sigma_limits(1);        
        DATA.(genvarname(molecule_name)).theta(frame) =     Localized(6)*theta_range+theta_limits(1);
        DATA.(genvarname(molecule_name)).frame_enabled(lower_frame_bin:upper_frame_bin) = 1;
        DATA.(genvarname(molecule_name)).localization_info(frame) = exitflag;
        DATA.(genvarname(molecule_name)).min_frame_used = min(frame,(DATA.(genvarname(molecule_name)).min_frame_used));
        DATA.(genvarname(molecule_name)).max_frame_used = max(frame,(DATA.(genvarname(molecule_name)).max_frame_used));
    end
end

function LOCALIZE_2X_GAUSSIANS(molecule_name,input_click1,input_click2)
    input_pixels = round(mean([input_click1;input_click2]));
%     start_row1 = input_click1(2); start_row2 = input_click2(2);
%     start_col1 = input_click1(1); start_col2 = input_click2(1);
    
    offset_position1 = ((input_click1 - input_pixels) + fitting_region + .5)/(2*fitting_region+1);
    offset_position2 = ((input_click2 - input_pixels) + fitting_region + .5)/(2*fitting_region+1);
    
    pixel_limits = {[input_pixels(2)-fitting_region,input_pixels(2)+fitting_region],[input_pixels(1)-fitting_region,input_pixels(1)+fitting_region]};
%     [roi_rows, roi_cols] = ndgrid(input_pixels(2)-fitting_region:input_pixels(2)+fitting_region,input_pixels(1)-fitting_region:input_pixels(1)+fitting_region);
    Image_Region = double(imread([FILEDATA.red_channel_filelocation FILEDATA.red_channel_filename],'Index',frame,'Info',FILEDATA.red_channel_fileinfo,'PixelRegion',pixel_limits));
    BG_guess = (mean(border_region(:).*Image_Region(:),'omitnan')-background_limits(1))/background_range;
    AMP_guess = (max(Image_Region(:))-mean(border_region(:).*Image_Region(:),'omitnan')-amplitude_limits(1))/amplitude_range;
    [Localized, resnorm, residual, exitflag] = lsqnonlin(@TWO_GAUSSIAN_FIT, [offset_position1(2), offset_position1(1), offset_position2(2), offset_position2(1), AMP_guess, AMP_guess, .5, .5,BG_guess], [0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1, 1], fitting_options);
    Row_1 = Localized(1)*(2*fitting_region+1)-.5-fitting_region;
    Col_1 = Localized(2)*(2*fitting_region+1)-.5-fitting_region;
    Row_2 = Localized(3)*(2*fitting_region+1)-.5-fitting_region;
    Col_2 = Localized(4)*(2*fitting_region+1)-.5-fitting_region;
    Amplitude_1 = Localized(5)*amplitude_range+amplitude_limits(1);
    Amplitude_2 = Localized(6)*amplitude_range+amplitude_limits(1);
    Sigma_1 = Localized(7)*sigma_range+sigma_limits(1);
    Sigma_2 = Localized(8)*sigma_range+sigma_limits(1);
    fit_localization = (Localized(9)*background_range+background_limits(1))...
        + Amplitude_1*exp(-((regional_indices_row-Row_1).^2)/(2*Sigma_1^2)-((regional_indices_col-Col_1).^2)/(2*Sigma_1^2))...
        + Amplitude_2*exp(-((regional_indices_row-Row_2).^2)/(2*Sigma_2^2)-((regional_indices_col-Col_2).^2)/(2*Sigma_2^2));
    if strcmp(molecule_name,current_molecule_name)==1
        fit_region_image = repmat(mat2gray(Image_Region),[1,1,3]);
        fit_localization_image = repmat(mat2gray(fit_localization),[1,1,3]);
        residual_in_gray = mat2gray(reshape(residual,[fitting_region*2+1,fitting_region*2+1]));
        set(lp.raw_data,'CData',fit_region_image);
        set(lp.fit_data,'CData',fit_localization_image);
        set(lp.residual_data,'CData',cat(3,residual_in_gray,residual_in_gray,residual_in_gray));
    end
    TWO_GAUSSIAN_LOCALIZTION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag);
end
function [Delta] = TWO_GAUSSIAN_FIT(guess)
    Row_1 = guess(1)*(2*fitting_region+1)-.5-fitting_region;
    Col_1 = guess(2)*(2*fitting_region+1)-.5-fitting_region;
    Row_2 = guess(3)*(2*fitting_region+1)-.5-fitting_region;
    Col_2 = guess(4)*(2*fitting_region+1)-.5-fitting_region;
    Amplitude_1 = guess(5)*amplitude_range+amplitude_limits(1);
    Amplitude_2 = guess(6)*amplitude_range+amplitude_limits(1);
    Sigma_1 = guess(7)*sigma_range+sigma_limits(1);
    Sigma_2 = guess(8)*sigma_range+sigma_limits(1);
    Background = guess(9)*background_range+background_limits(1);
    Guess_Image = Background...
        + Amplitude_1*exp(-((regional_indices_row-Row_1).^2+(regional_indices_col-Col_1).^2)/(2*Sigma_1^2))...
        + Amplitude_2*exp(-((regional_indices_row-Row_2).^2+(regional_indices_col-Col_2).^2)/(2*Sigma_2^2));
    Delta = - Image_Region(:) + Guess_Image(:);
%     Delta = sum(- Image_Region(:) + Guess_Image(:)).^2;
end

function TWO_GAUSSIAN_LOCALIZTION_QUALITY_CHECK(molecule_name,Localized,fit_localization,resnorm,input_pixels,exitflag)
    if (Localized(7)*sigma_range+sigma_limits(1)) < .5 || (Localized(8)*sigma_range+sigma_limits(1)) < .5 || min(Localized(5)/Localized(6),Localized(6)/Localized(5)<.2) || min(Localized(7)/Localized(8),Localized(8)/Localized(7)<.2)
        % add in more conditions here.
        if DATA.(genvarname(molecule_name)).localization_warning == 0
           DATA.(genvarname(molecule_name)).localization_warning = 1;
           % IF BAD, BUT NO WARNING ADD A WARNING
           % Find the local peaks and try and relocalize with those.
           [rowMax, colMax] = ind2sub([2*fitting_region+1,2*fitting_region+1],find(imregionalmax(Image_Region))); % value
           sortedPeaks = sortrows([rowMax+input_pixels(2)-fitting_region,colMax+input_pixels(1)-fitting_region,Image_Region(sub2ind(size(Image_Region),rowMax,colMax))],-3);
           input_click1 = [sortedPeaks(1,2),sortedPeaks(1,1)];
           input_click2 = [sortedPeaks(2,2),sortedPeaks(2,1)];
           LOCALIZE_2X_GAUSSIANS(molecule_name,input_click1,input_click2)
        else
            % IF BAD, WITH WARNING, REMOVE TRACK
            Unselected = cellstr(get(Unselected_molecules,'String'));
            Selected = cellstr(get(Selected_molecules, 'string'));   
            new_clicked = find(strcmp(molecule_name,Selected)==1, 1);
            Current_selection = Selected(new_clicked);
            Unselected = [Unselected; Current_selection];
            Selected(new_clicked,:) = [];
            Unselected(any(cellfun(@isempty,Unselected),2),:) = [];
            Selected(any(cellfun(@isempty,Selected),2),:) = [];
            set(Unselected_molecules,'String', sort(Unselected),'Value',1)
            set(Selected_molecules,'String', sort(Selected),'Value',1)
            if isempty(Selected) == 1
                % IF ONLY TRACK, STOP RUNNING
                toggle_autorun = [0 0];
            end
        end
    else
        % IF GOOD, REMOVE A WARNING
        DATA.(genvarname(molecule_name)).localization_warning = 0;
        DATA.(genvarname(molecule_name)).signal_above_BG(frame) = sum(Image_Region(:)-Localized(9)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).fitted_signal(frame) = sum(fit_localization(:)-Localized(9)*background_range-background_limits(1));
        DATA.(genvarname(molecule_name)).residual(frame) = resnorm;
        DATA.(genvarname(molecule_name)).row_1(frame) = input_pixels(2)+Localized(1)*(2*fitting_region+1)-.5-fitting_region;
        DATA.(genvarname(molecule_name)).col_1(frame) = input_pixels(1)+Localized(2)*(2*fitting_region+1)-.5-fitting_region;
        DATA.(genvarname(molecule_name)).row_2(frame) = input_pixels(2)+Localized(3)*(2*fitting_region+1)-.5-fitting_region;
        DATA.(genvarname(molecule_name)).col_2(frame) = input_pixels(1)+Localized(4)*(2*fitting_region+1)-.5-fitting_region;
        DATA.(genvarname(molecule_name)).theta(frame) = atan2(Localized(1)-Localized(3),Localized(2)-Localized(4));  
        while DATA.(genvarname(molecule_name)).theta(frame) < 0
            DATA.(genvarname(molecule_name)).theta(frame) = DATA.(genvarname(molecule_name)).theta(frame)+pi;
        end
        DATA.(genvarname(molecule_name)).row(frame) = mean([DATA.(genvarname(molecule_name)).row_1(frame),DATA.(genvarname(molecule_name)).row_2(frame)]);
        DATA.(genvarname(molecule_name)).col(frame) = mean([DATA.(genvarname(molecule_name)).col_1(frame),DATA.(genvarname(molecule_name)).col_2(frame)]);
        DATA.(genvarname(molecule_name)).amplitude_1(frame) = Localized(5)*amplitude_range+amplitude_limits(1);
        DATA.(genvarname(molecule_name)).amplitude_2(frame) = Localized(6)*amplitude_range+amplitude_limits(1);
        DATA.(genvarname(molecule_name)).background(frame) = Localized(9)*background_range+background_limits(1);
        DATA.(genvarname(molecule_name)).sigma_1(frame) = Localized(7)*sigma_range+sigma_limits(1);
        DATA.(genvarname(molecule_name)).sigma_2(frame) = Localized(8)*sigma_range+sigma_limits(1);        
        DATA.(genvarname(molecule_name)).frame_enabled(lower_frame_bin:upper_frame_bin) = 1;
        DATA.(genvarname(molecule_name)).localization_info(frame) = exitflag;
        DATA.(genvarname(molecule_name)).min_frame_used = min(frame,(DATA.(genvarname(molecule_name)).min_frame_used));
        DATA.(genvarname(molecule_name)).max_frame_used = max(frame,(DATA.(genvarname(molecule_name)).max_frame_used));
    end
end

end