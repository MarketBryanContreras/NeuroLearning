 %% Super_BC_ phase_mod

%% initialize

% Make sure to add code from the CEH2 repo, Vadermeer lab code-shared repo,
% NeuroLearning repo before going forward
 %% Dynamic loader
[data_dir, inter_dir, mouse_group]=BC_linearTrack_dynamicLoader('experimental'); %"experimental" for archT, "control" 
%% Parameters
plot_flag = 1; % switch to 0 if you want to supress verification figures.
time_maze_start = 30; %Seconds to exclude from recording
min_trial_dur = 0.5;
%mouse_group=1; %1 for ArchT and 2 for eYFP. This just modify color of the plots
save_flag=00;
if plot_flag==0 %Making sure that if there is not plots, the save flag is off
    save_flag=00;
end
%% Generating color pallet for the mouse group
if mouse_group==1
    gcolors = [BC_color_genertor('powder_blue');  % RGB values for Group 2
        BC_color_genertor('archt_green')]; % RGB values for Group1;
else
    gcolors = [BC_color_genertor('torment_blue');  % RGB values for Group 2
        BC_color_genertor('swamp_green')]; % RGB values for Group1;
end

%% Removing non- desired recordings
inhib_dir = dir('*D4_INHIBITION*');% get all sessions with 'D4_INHIBITION'
indices_to_remove = contains({inhib_dir.name}, 'disconnected');% Find the indices of directories containing 'disconnected'
inhib_dir(indices_to_remove) = [];% Remove the directories containing 'disconnected'

%% Loop to load data from raw

for iS =1%:length(inhib_dir)
   
    %% loading
    cd([inhib_dir(iS).folder filesep inhib_dir(iS).name])
    parts = strsplit(inhib_dir(iS).name, '_');
    %Grabbing some info from the file name
    info.subject = parts{1};
    info.sess = parts{5};
    info.date = [ parts{2} '_' parts{3} '_' parts{4}];    
    cfg_csc = [];
    cfg_csc.desired_sampling_frequency = 2000;
    %Assign which csc to load to each mouse and the pattern of the ttl that represents an event in the board
    if strcmpi(info.subject, 'BC1602')
        cfg_csc.fc ={'CSC4.ncs'}; %2#'CSC4.ncs,7'
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    elseif strcmpi(info.subject, 'BC051')
        cfg_csc.fc ={'CSC5.ncs'};%5%7%2%4
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    elseif strcmpi(info.subject, 'BC1807')
        cfg_csc.fc ={'CSC6.ncs'};%3
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    elseif strcmpi(info.subject, 'BC053')
        cfg_csc.fc ={'CSC4.ncs'};%4
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).';
    elseif strcmpi(info.subject, 'BC054')
        cfg_csc.fc ={'CSC5.ncs'};%5%2%4%7
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).';
    elseif strcmpi(info.subject, 'BC011')
        cfg_csc.fc ={'CSC3.ncs'};%3%7
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    elseif strcmpi(info.subject, 'BC013')
        cfg_csc.fc ={'CSC6.ncs'};%2%4
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).';
    elseif strcmpi(info.subject, 'BC014')
        cfg_csc.fc ={'CSC7.ncs'};%4%5
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).';
    end

    [csc, evts, pos] = BC_load_NLX(cfg_csc);%Load ,csc, events and position
    fs=csc.cfg.hdr{1}.SamplingFrequency;

    %Filtering the speed to remove outliers
    % filtWinSz=60; %Size of the window in frames to  
    % pos.data(5,:)=medfilt1(pos.data(5,:),filtWinSz);
    spd=pos.data(5,:);
    %histogram(spd)
    out_idx=(spd> 35);
    spd(out_idx)=NaN;
    %histogram(spd)
    %plot(isnan(spd))
    spd=fillmissing(spd,'movmedian',60);
    %plot(isnan(spd))
    pos.data(5,:)=spd;
    % figure(2)
    % plot(pos.tvec,pos.data(5,:))
    %Restricting the csc to the values from time_maze_start onwards on the recording
    fs = csc.cfg.hdr{1}.SamplingFrequency;
    sSampCsc= time_maze_start*fs;                                          %This is the idx of the samplimg at sec 30
    restrictIvCsc=iv(csc.tvec(sSampCsc),csc.tvec(end));                    %Creates and IV from the time the mouse was placed in the maze to the end of the recording
    csc= restrict(csc,restrictIvCsc);                                      %Restrict csc to the previously created csc
    % Correction of time
    strt=csc.tvec(1);
    csc.tvec= csc.tvec-strt;                                               %correct time csc
    nn=size(evts.t,2);
    evts.t{1,nn-1}=evts.t{1,nn-1}-strt;                                    % Coreccting laser events end
    evts.t{1,nn}=evts.t{1,nn}-strt;                                        % Coreccting laser events start
    
    %Restiction of time when the mouse is placed in the maze for position,different sampling rate
    sSamp=time_maze_start*30;
    restrictIvPos=iv(pos.tvec(sSamp),pos.tvec(end));
    pos=restrict(pos,restrictIvPos);
    pos.tvec= pos.tvec-pos.tvec(1);

    %% Get period of time where the mouse is running and there is inhibition
    %This creates the intervals where the mouse received light
    iv_inhb = MS_get_evts_off(evts, pattern);
    [iv_inhb,iv_noInhb,iv_running] = BC_LT_trialfun(pos, iv_inhb, plot_flag); %This generates the graph of inhibition and velocity while it also retuns the inhibition, no inhibition and runnings epochs 
    
    if save_flag
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0001' filesep 'Fig01_' info.subject '.png' ]);
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0001' filesep 'Fig01_' info.subject '.pdf' ]);
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0001' filesep 'Fig01_' info.subject '.fig' ]);
    end
    %% Filter
    % filter the LFP in the theta band
    cfg_filt_t = [];
    cfg_filt_t.type = 'cheby1';                                            %'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [4 12];                                                % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 3;                                                  %type filter order
    cfg_filt_t.display_filter = 0;                                         % use this to see the fvtool
    theta_csc = FilterLFP(cfg_filt_t, csc);                                % filter the raw LFP using
    theta_amp_phi = BC_power_phase(theta_csc);
    theta_decimated=BC_ampSpdDcmt(theta_csc, pos);
    % Filter the LFP in the slow gamma band
    cfg_filt_t = [];
    cfg_filt_t.type = 'butter';                                            %'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [30 58];                                               % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 4;                                                  %type filter order
    cfg_filt_t.display_filter = 0;                                         % use this to see the fvtool
    SG_csc = FilterLFP(cfg_filt_t, csc);                                   % filter the raw LFP using
    SG_amp_phi = BC_power_phase(SG_csc);  
    SG_decimated=BC_ampSpdDcmt(SG_csc, pos);
    % Filter the LFP in the fast gamma band
    cfg_filt_t = [];
    cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [60 100]; % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 4; %type filter order
    cfg_filt_t.display_filter = 0; % use this to see the fvtool
    FG_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
    FG_amp_phi = BC_power_phase(FG_csc); 
    FG_decimated=BC_ampSpdDcmt(FG_csc, pos);
    %% Calculating the psd for the frequencies of interest
    %Caluclate the psd
    w = 2^12; [pxx,f]=pwelch(csc.data, hanning(w),w/2,w*2, fs);
    %Calculate the mean of the 1-50 Hz
    freqRange = [5 50];
    % Find the indices of frequencies within the desired range
    idx = f >= freqRange(1) & f <= freqRange(2);
    % Extract the frequencies and power values within this range
    frequencies_in_range = f(idx);
    power_in_range = pxx(idx);
    % Sum the power to get total power within this range
    total_power_in_range = sum(power_in_range);
    %Divide the frequencies to the mean amplitude of 1-50
    norm_power=pxx./total_power_in_range
    % Plot
    if plot_flag
        plot(f,(norm_power)); xlim([0 15])
    end
    %store the values in a strucutre
    psd=[];
    psd.frequency=f;
    psd.power=pxx;
    psd.normPower=norm_power;
    %% Plotting ampplitude and velocity according to position in the linear track
if plot_flag
    % Plot the x_pos by time and plot the int to collaborate that you got them right
    x_data= theta_decimated.data(5,:);
    fig=figure(0002);
        clf;
        x = (theta_decimated.tvec)';
        y = x_data;
        z = zeros(size(x));
        speed_data=theta_decimated.data(4,:);
        col = speed_data;  % This is the color, it varies with x in this case.
        ax1=subplot(2,1,1);
        surface([x;x],[y;y],[z;z],[col;col],... % This create a surface in the shape of a linear plot, this way it allows to mdify the color of the line accoridn to the speed of the mouse
            'facecol','no',...
            'edgecol','interp',...
            'linew',2);
        c=colorbar;
        c.Label.String='Velocity (cm/s)';
        c.Label.FontSize=12;
        c.Position= [0.910 0.5900 0.0100 0.330];
        xlim([pos.tvec(1) pos.tvec(end)]);
        set(gca, 'TickDir', 'out');
        ax1=subplot(2,1,1);
        %Amplitude
        subplot(2,1,2)
        x_data= theta_decimated.data(5,:);
       
        
        x = (theta_decimated.tvec)';
        y = x_data;
        z = zeros(size(x));
        amp_data=theta_decimated.data(1,:);
        col = amp_data;  % This is the color, it varies with x in this case.
        ax2=subplot(2,1,2);
        surface([x;x],[y;y],[z;z],[col;col],... % This create a surface in the shape of a linear plot, this way it allows to mdify the color of the line accoridn to the speed of the mouse
            'facecol','no',...
            'edgecol','interp',...
            'linew',2);
        c2=colorbar;
        c2.Label.String='Tta Amplitude Normalized';
        c2.Label.FontSize=12;
        c2.Position= [0.910 0.1200 0.0100 0.330];
        xlim([pos.tvec(1) pos.tvec(end)]);
        set(gca, 'TickDir', 'out');
        %  (*) Asthetics (*)
        %Give common xlabel, ylabel, and title
        han=axes(fig,'visible','off');
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        han.YLabel.Position=[-0.0300 0.5000 0];
        ylabel(han,'Position in the x axis (cm)','FontSize', 16);
        han.XLabel.Position=[0.5000 -0.060 0];
        xlabel(han,'Time(s)','FontSize', 16);
        title(han,'Example of mouse displacement during linear track','FontSize', 20);
        fig = gcf;                   % Get current figure handle
        fig.Color = [1 1 1];
        fig.Color = [1 1 1];         % Set background color to white
        fig.Position = [100, 100, 1600, 700];  % [x, y, width, height]
        hold off;
        % figure(1921)
        % scatter(theta_decimated.data(4,:),theta_decimated.data(1,:)) %This plots the speed vs the theta amplitude
        % R=corr2((theta_decimated.data(4,:)),(theta_decimated.data(1,:))) %This calculates the correlation between speed and amplitude but I have to correct for nan values I think
        % lsline
end


if save_flag
    saveas(figure(2),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0002' filesep 'Fig02_' info.subject '.png' ]);
    saveas(figure(2),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0002' filesep 'Fig02_' info.subject '.pdf' ]);
    saveas(figure(2),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0002' filesep 'Fig02_' info.subject '.fig' ]);
end
    %% Restricting data to intervals of inhb, noInhb
    %inhb
    csc_inhb = restrict(csc, iv_inhb);
    theta_inhb=restrict(theta_csc, iv_inhb);
    SG_inhb = restrict(SG_csc, iv_inhb);
    FG_inhb = restrict(FG_csc, iv_inhb);
    pos_inhb= restrict(pos, iv_inhb);
    thetaD_inhb=restrict(theta_decimated, iv_inhb);
    %no-inhb
    csc_noinhb = restrict(csc, iv_noInhb);
    theta_noinhb = restrict(theta_csc, iv_noInhb);
    SG_noinhb = restrict(SG_csc, iv_noInhb);
    FG_noinhb = restrict(FG_csc, iv_noInhb);
    pos_noinhb= restrict(pos, iv_noInhb);
    thetaD_noinhb=restrict(theta_decimated, iv_noInhb);
    %running
    csc_running = restrict(csc, iv_running);
    theta_running = restrict(theta_csc, iv_running);
    SG_running = restrict(SG_csc, iv_running);
    FG_running = restrict(FG_csc, iv_running);
    pos_running= restrict(pos, iv_running);
   %% ploting the speed vs amplitude in inhibition and no inhibition
   % tresh=0;
   % SpdInhb=thetaD_inhb.data(4,:);
   % AmpInhb=thetaD_inhb.data(1,:);
   % outIdxInhb=SpdInhb>tresh;
   % SpdInhb=SpdInhb(outIdxInhb);
   % AmpInhb=AmpInhb(outIdxInhb);
   % 
   % SpdNoInhb=thetaD_noinhb.data(4,:);
   % outIdxNoInhb=SpdNoInhb>tresh;
   % AmpNoInhb=thetaD_noinhb.data(1,:);
   % SpdNoInhb=SpdNoInhb(outIdxNoInhb);
   % AmpNoInhb=AmpNoInhb(outIdxNoInhb);
   % 
   % if plot_flag
   %     figure(1001)
   %     s1=subplot(1,2,1);
   %     scatter(SpdInhb,AmpInhb, 'MarkerEdgeColor',BC_color_genertor('Archt_green'),'SizeData',11)
   %     lsq1=lsline(s1);lsq1.Color='g',lsq1.LineWidth=2;
   %     %xlim([0 50]); ylim([0 1]);
   %     xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
   %     RI=corr2(SpdInhb,AmpInhb)
   %     title(sprintf('Silencing R=%f',RI))
   % 
   %     s2=subplot(1,2,2);
   %     scatter(SpdNoInhb,AmpNoInhb,'MarkerEdgeColor',BC_color_genertor('Oxford_blue'),'SizeData',11)
   %     lsq2=lsline(s2);lsq2.Color='g',lsq2.LineWidth=2;
   %     %xlim([0 50]); ylim([0 1]);
   %     xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
   %     RNI=corr2((thetaD_noinhb.data(4,:)),(thetaD_noinhb.data(1,:)))
   %     title(sprintf('No Silencing R=%f',RNI))
   % 
   % end
    %% Calculating and plotting modulation index 
    %For slow gamma
    %SG-Inhb
    SGphiAmpNormInhb=BC_phase_amp_norm_bins(theta_inhb,SG_inhb);
    SGInhb_modidx=MS_ModIdx(SGphiAmpNormInhb);
    if plot_flag
        if mouse_group ==1;
            BC_plot_modidx(SGphiAmpNormInhb,BC_color_genertor('Archt_green'),SGInhb_modidx,'SG', 'Light' )
        else
            BC_plot_modidx(SGphiAmpNormInhb,BC_color_genertor('Swamp_green'),SGInhb_modidx,'SG', 'Light' )
        end
    end
    if save_flag
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0003' filesep 'Fig03_' info.subject '.png' ]);
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0003' filesep 'Fig03_' info.subject '.pdf' ]);
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0003' filesep 'Fig03_' info.subject '.fig' ]);
    end
%SG-NoInhb
    SGphiAmpNormNoInhb = BC_phase_amp_norm_bins(theta_noinhb,SG_noinhb);
    SGNoInhb_modidx=MS_ModIdx(SGphiAmpNormNoInhb);
    if plot_flag
        if mouse_group ==1;
        BC_plot_modidx(SGphiAmpNormNoInhb,BC_color_genertor('Powder_blue'),SGNoInhb_modidx, 'SG', 'No light')
        else
        BC_plot_modidx(SGphiAmpNormNoInhb,BC_color_genertor('Torment_blue'),SGNoInhb_modidx, 'SG', 'No light')
        end 
    end
    if save_flag
        saveas(figure(2),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0004' filesep 'Fig04_' info.subject '.png' ]);
        saveas(figure(2),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0004' filesep 'Fig04_' info.subject '.pdf' ]);
        saveas(figure(2),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0004' filesep 'Fig04_' info.subject '.fig' ]);
    end
    [SGshift_mean, SGshift_std]=LTshifted_meanModIdx(theta_running,SG_running);
%Z scores for SG. This calculates the z score with the mean and the std from the random shift
    z_SGInhb_modidx = (SGInhb_modidx - SGshift_mean) / SGshift_std;
    z_SGNoInhb_modidx = (SGNoInhb_modidx - SGshift_mean) / SGshift_std;
% %Prototype of a bar graph for the mod idx comparison between light and no light
% if plot_flag
%     figure(19)
%     bSG=bar(["No Light" "Light"],[0.3*10^-4  ;0.35*10^-4 ],'FaceColor','flat');
%     % Set the face color of each bar individually
%     for k = 1:2 % Loop through the number of columns in data
%         bSG.CData(k,:) = gcolors(k,:);
%     end
%     ylabel('Modulation Index');
%     set(gca,'fontsize', 14);
%     bSG.EdgeColor = 'none';
%     set(gca, 'TickDir', 'out');  % Move ticks outside the plot
%     set(gca,'box','off');
%     fig = gcf;                   % Get current figure handle
%     fig.Color = [1 1 1];         % Set background color to white
%     fig.Position = [100, 100, 600, 800];  % [x, y, width, height]
% end

%For Fast Gamma
%FG-Inhb
    FGphiAmpNormInhb=BC_phase_amp_norm_bins(theta_inhb,FG_inhb);
    FGInhb_modidx=MS_ModIdx(FGphiAmpNormInhb);
    if plot_flag
        if mouse_group ==1;
        BC_plot_modidx(FGphiAmpNormInhb,BC_color_genertor('Archt_green'),FGInhb_modidx,'FG', 'Light' )
        else
        BC_plot_modidx(FGphiAmpNormInhb,BC_color_genertor('Swamp_green'),FGInhb_modidx,'FG', 'Light' )
        end
    end
      if save_flag
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0005' filesep 'Fig05_' info.subject '.png' ]);
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0005' filesep 'Fig05_' info.subject '.pdf' ]);
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0005' filesep 'Fig05_' info.subject '.fig' ]);
    end
%%FG-NoInhb
    FGphiAmpNormNoInhb = BC_phase_amp_norm_bins(theta_noinhb,FG_noinhb);
    FGNoInhb_modidx=MS_ModIdx(FGphiAmpNormNoInhb);
    if plot_flag
        if mouse_group ==1;
        BC_plot_modidx(FGphiAmpNormNoInhb,BC_color_genertor('Powder_blue'),FGNoInhb_modidx, 'FG', 'No light')
        else
        BC_plot_modidx(FGphiAmpNormNoInhb,BC_color_genertor('Torment_blue'),FGNoInhb_modidx, 'FG', 'No light')
        end
    end
    if save_flag
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0006' filesep 'Fig06_' info.subject '.png' ]);
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0006' filesep 'Fig06_' info.subject '.pdf' ]);
        saveas(figure(1),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0006' filesep 'Fig06_' info.subject '.fig' ]);
    end
    %Shifted
    [FGshift_mean,FGshift_std]=LTshifted_meanModIdx(theta_running,FG_running);
    %Z scores for FG

    z_FGInhb_modidx = (FGInhb_modidx - FGshift_mean) / FGshift_std; 
    z_FGNoInhb_modidx = (FGNoInhb_modidx - FGshift_mean) / FGshift_std;  
    %% 
    % Identify the csc corresponding to silencing and no silencing

    %Apply psd and reduce the xlim to 0 to 15


    %Calculate the band power of delta, theta, sg and fg


    %Normalize to the band power to 140
    %% Epochs by epochs spectral analysis
    
    % inhibition
    win_s = 256;
    r_inhib = [];
    r_noinhib = [];
    r_run = [];
    spd_ihb=[];
    tta_amp_ihb= [];
    tta_pwr_ihb= [];
    tta_phi_ihb= [];
    sg_amp_ihb= [];
    sg_pwr_ihb= [];
    sg_phi_ihb= [];
    fg_amp_ihb= [];
    fg_pwr_ihb= [];
    fg_phi_ihb= [];
    epoch_ihb= [];

    tta_amp_ihbTotal= [];
    tta_pwr_ihbTotal= [];
    tta_phi_ihbTotal= [];
    sg_amp_ihbTotal= [];
    sg_pwr_ihbTotal= [];
    sg_phi_ihbTotal= [];
    fg_amp_ihbTotal= [];
    fg_pwr_ihbTotal= [];
    fg_phi_ihbTotal= [];
    epoch_ihbTotal= [];
    inhbTbl=[];

    psd_inhb=[];
    psd_inhb_norm=[];

    for ii = length(iv_inhb.tstart):-1:1
        if (iv_inhb.tend(ii) - iv_inhb.tstart(ii)) < min_trial_dur %Testing that the epoch last more than the treshold
            t_bp_inhib(ii)= NaN;
            sg_bp_inhib(ii) = NaN;
            fg_bp_inhib(ii) = NaN;
            modidx_SG_inhib(ii)= NaN;
            modidx_FG_inhib(ii)= NaN;
            r_inhib = NaN(length(1:0.1:160), length(1:0.1:160));
            spd_ihb=[spd_ihb NaN];
            tta_amp_ihb= [tta_amp_ihb NaN];
            tta_pwr_ihb= [tta_pwr_ihb NaN];
            tta_phi_ihb= [tta_phi_ihb NaN];
            sg_amp_ihb= [sg_amp_ihb NaN];
            sg_pwr_ihb= [sg_pwr_ihb NaN];
            sg_phi_ihb= [sg_phi_ihb NaN];
            fg_amp_ihb= [fg_amp_ihb NaN];
            fg_pwr_ihb= [fg_pwr_ihb NaN];
            fg_phi_ihb= [fg_phi_ihb NaN];
            epoch_ihb= [epoch_ihb ii];

            tta_amp_ihbTotal= [tta_amp_ihbTotal NaN];
            tta_pwr_ihbTotal= [tta_pwr_ihbTotal NaN];
            tta_phi_ihbTotal= [tta_phi_ihbTotal NaN];
            sg_amp_ihbTotal= [sg_amp_ihbTotal NaN];
            sg_pwr_ihbTotal= [sg_pwr_ihbTotal NaN];
            sg_phi_ihbTotal= [sg_phi_ihbTotal NaN];
            fg_amp_ihbTotal= [fg_amp_ihbTotal NaN];
            fg_pwr_ihbTotal= [fg_pwr_ihbTotal NaN];
            fg_phi_ihbTotal= [fg_phi_ihbTotal NaN];
            epoch_ihbTotal= [epoch_ihbTotal ii];

            psd_inhb=[psd_inhb; NaN];

            outTbl=BC_EpochSpdAmpDcmtdTbl(NaN, NaN, NaN,iv_inhb,ii);
        else
            this_csc = restrict(csc, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            
            % get the power
            t_bp =  bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [4 12]); %Calculates the theta band power
            sg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [30 58]);%Calculates the sg band power
            fg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [60 100]);%Calculates the fg band power
            
            ref_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 50]);%Calculates the ref band power
            
            %Store the theta, sg, and fg power in their varibales
            t_bp_inhib(ii)= t_bp;
            sg_bp_inhib(ii) = sg_bp;
            fg_bp_inhib(ii) = fg_bp;
            
            t_bp_inhib_norm(ii)= t_bp/ ref_bp;
            sg_bp_inhib_norm(ii) = sg_bp/ ref_bp;
            fg_bp_inhib_norm(ii) = fg_bp/ ref_bp;

            %Get the PSD
            %Caluclate the psd
            w = 2^12; [pxx,f]=pwelch(csc.data, hanning(w),w/2,w*2, fs);
            psd_inhb=[psd_inhb; pxx'];
            %Calculate the mean of the 1-50 Hz
            freqRange = [5 50];
            % Find the indices of frequencies within the desired range
            idx = f >= freqRange(1) & f <= freqRange(2);
            % Extract the frequencies and power values within this range
            frequencies_in_range = f(idx);
            power_in_range = pxx(idx);
            % Sum the power to get total power within this range
            total_power_in_range = sum(power_in_range);
            %Divide the frequencies to the mean amplitude of 1-50
            norm_power=pxx./total_power_in_range;
            psd_inhb_norm=[psd_inhb_norm; norm_power'];

            %Create a table with epoch#, condition, amplitude, velocity
            outTbl=BC_EpochSpdAmpDcmtdTbl(theta_decimated, SG_decimated,FG_decimated,iv_inhb,ii);       
            %Create a table with epoch#, condition, amplitude, velocity
            %The following code create the veactor for the creation of the table
            this_theta_amp=restrict(theta_amp_phi, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            this_Sg_amp=restrict(SG_amp_phi, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            this_Fg_amp=restrict(FG_amp_phi, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            this_pos= restrict(pos, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            
            %Creating a copy of the amplitude vectors so I can use all the values in a second analysis
            this_theta_amp2=this_theta_amp;
            this_Sg_amp2=this_Sg_amp;
            this_Fg_amp2=this_Fg_amp;
            %Combining to a single structure
            this_theta_amp2.data(4:6,:)=this_Sg_amp2.data;
            this_theta_amp2.data(7:9,:)=this_Fg_amp2.data;
            this_theta_amp2.label={'tta_amp','tta_pwr','tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi'};
            %Extracting into vectors for setting a table later
            tta_amp_ihbTotal= [tta_amp_ihbTotal this_theta_amp2.data(1,:)];
            tta_pwr_ihbTotal= [tta_pwr_ihbTotal this_theta_amp2.data(2,:)];
            tta_phi_ihbTotal= [tta_phi_ihbTotal this_theta_amp2.data(3,:)];
            sg_amp_ihbTotal= [sg_amp_ihbTotal this_theta_amp2.data(4,:)];
            sg_pwr_ihbTotal= [sg_pwr_ihbTotal this_theta_amp2.data(5,:)];
            sg_phi_ihbTotal= [sg_phi_ihbTotal this_theta_amp2.data(6,:)];
            fg_amp_ihbTotal= [fg_amp_ihbTotal this_theta_amp2.data(7,:)];
            fg_pwr_ihbTotal= [fg_pwr_ihbTotal this_theta_amp2.data(8,:)];
            fg_phi_ihbTotal= [fg_phi_ihbTotal this_theta_amp2.data(9,:)];
            epoch_ihbTotal= [epoch_ihbTotal repmat(ii,1,length(this_theta_amp2.data))];

            %Continuation of the speed-amplitude analysis
            this_speed=this_pos;
            this_speed.data=[];this_speed.data=this_pos.data(5,:);this_speed.label={'speed'};this_speed.units={'cm/s'};
            extraction_idx=nearest_idx(this_speed.tvec, this_theta_amp.tvec); %Get the idx of the amplitude which tvec is more similar to those of the speed tvec
            %Prunning the amplitude to only the time points shared with the velocity vector
            this_theta_amp.data=this_theta_amp.data(1:3,extraction_idx);
            this_Sg_amp.data=this_Sg_amp.data(1:3,extraction_idx);
            this_Fg_amp.data=this_Fg_amp.data(1:3,extraction_idx);
            %Combining everyhing into a single structure
            epoch_spd_amp_phi=[];epoch_spd_amp_phi.type='tsd';epoch_spd_amp_phi.tvec=this_speed.tvec;epoch_spd_amp_phi.cfg=this_theta_amp.cfg;
            epoch_spd_amp_phi.data(1,:)=this_speed.data;
            epoch_spd_amp_phi.data(2:4,:)=this_theta_amp.data;
            epoch_spd_amp_phi.data(5:7,:)=this_Sg_amp.data;
            epoch_spd_amp_phi.data(8:10,:)=this_Fg_amp.data;
            epoch_spd_amp_phi.data(11,:)=repmat(ii,1,length(this_speed.tvec));
            epoch_spd_amp_phi.label={'spd', 'tta_amp','tta_pwr','tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi', '#epoch'};
            %Extracting the data from the structure and combining vector for all the epochs
            spd_ihb=[spd_ihb epoch_spd_amp_phi.data(1,:)];
            tta_amp_ihb= [tta_amp_ihb epoch_spd_amp_phi.data(2,:)];
            tta_pwr_ihb= [tta_pwr_ihb epoch_spd_amp_phi.data(3,:)];
            tta_phi_ihb= [tta_phi_ihb epoch_spd_amp_phi.data(4,:)];
            sg_amp_ihb= [sg_amp_ihb epoch_spd_amp_phi.data(5,:)];
            sg_pwr_ihb= [sg_pwr_ihb epoch_spd_amp_phi.data(6,:)];
            sg_phi_ihb= [sg_phi_ihb epoch_spd_amp_phi.data(7,:)];
            fg_amp_ihb= [fg_amp_ihb epoch_spd_amp_phi.data(8,:)];
            fg_pwr_ihb= [fg_pwr_ihb epoch_spd_amp_phi.data(9,:)];
            fg_phi_ihb= [fg_phi_ihb epoch_spd_amp_phi.data(10,:)]; 
            epoch_ihb= [epoch_ihb epoch_spd_amp_phi.data(11,:)];


            % Phase mod SG and FG
            this_th = restrict(theta_csc, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            this_sg = restrict(SG_csc, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            this_fg = restrict(FG_csc, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            %SG
            this_SG_phi_amp=BC_phase_amp_norm_bins(this_th,this_sg);
            modidx_SG_inhib(ii) =MS_ModIdx(this_SG_phi_amp);
            %FG
            this_FG_phi_amp=BC_phase_amp_norm_bins(this_th,this_fg);
            modidx_FG_inhib(ii) =MS_ModIdx(this_FG_phi_amp);
            % cross freq coupling
            [~, F, ~,P] = spectrogram(this_csc.data,hanning(win_s),win_s/2,1:0.1:160,this_csc.cfg.hdr{1}.SamplingFrequency); % spectrogram -- will take a while to compute!, add 'yaxis for visualization porpuses
            %F: Frequency axis of the spectrogram. P: Power spectral density (or magnitude squared of the Fourier transform) of the signal.
            [r_inhib(:,:,ii),~] = corrcoef(10*log10(P')); % correlation matrix (across frequ`encies) of spectrogram
        end
        inhbTbl=[inhbTbl;outTbl];
    end
    %putting the speed-pwr into a table so its easier to extract later
    spdPwrTbl_inhb=table( spd_ihb',tta_amp_ihb',tta_pwr_ihb', tta_phi_ihb',sg_amp_ihb',sg_pwr_ihb',sg_phi_ihb',fg_amp_ihb',fg_pwr_ihb',fg_phi_ihb',epoch_ihb', 'VariableNames',{'spd','tta_amp','tta_pwr', 'tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi','epoch'});
    %creating a table for all the variables of the amplitude        
    spdPwrTbl_inhbTotal= table( tta_amp_ihbTotal',tta_pwr_ihbTotal', tta_phi_ihbTotal',sg_amp_ihbTotal',sg_pwr_ihbTotal',sg_phi_ihbTotal',fg_amp_ihbTotal',fg_pwr_ihbTotal',fg_phi_ihbTotal',epoch_ihbTotal', 'VariableNames',{'tta_amp','tta_pwr', 'tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi','epoch'});      
%Plotting the psd with stdandar error
std_error_psd_norm=std(psd_inhb_norm)./sqrt(size(psd_inhb_norm,2));



    % NO inhibition
            spd_noIhb=[]; 
            tta_amp_noIhb= [];
            tta_pwr_noIhb= [];
            tta_phi_noIhb= [];
            sg_amp_noIhb= [];
            sg_pwr_noIhb= [];
            sg_phi_noIhb= [];
            fg_amp_noIhb= [];
            fg_pwr_noIhb= [];
            fg_phi_noIhb= []; 
            epoch_noIhb= [];
            tta_amp_noIhbTotal= [];
            tta_pwr_noIhbTotal= [];
            tta_phi_noIhbTotal= [];
            sg_amp_noIhbTotal= [];
            sg_pwr_noIhbTotal= [];
            sg_phi_noIhbTotal= [];
            fg_amp_noIhbTotal= [];
            fg_pwr_noIhbTotal= [];
            fg_phi_noIhbTotal= []; 
            epoch_noIhbTotal= [];

            NoInhbTbl=[];
            
    for ii = length(iv_noInhb.tstart):-1:1
        if (iv_noInhb.tend(ii) - iv_noInhb.tstart(ii)) < min_trial_dur
            t_bp_noinhib(ii)= NaN;
            sg_bp_noinhib(ii) = NaN;
            fg_bp_noinhib(ii) = NaN;
            modidx_SG_noinhib(ii)= NaN;
            modidx_FG_noinhib(ii)= NaN;
            r_noinhib = NaN(length(1:0.1:160), length(1:0.1:160));
            spd_noIhb=[spd_noIhb NaN];
            tta_amp_noIhb= [tta_amp_noIhb NaN];
            tta_pwr_noIhb= [tta_pwr_noIhb NaN];
            tta_phi_noIhb= [tta_phi_noIhb NaN];
            sg_amp_noIhb= [sg_amp_noIhb NaN];
            sg_pwr_noIhb= [sg_pwr_noIhb NaN];
            sg_phi_noIhb= [sg_phi_noIhb NaN];
            fg_amp_noIhb= [fg_amp_noIhb NaN];
            fg_pwr_noIhb= [fg_pwr_noIhb NaN];
            fg_phi_noIhb= [fg_phi_noIhb NaN]; 
            epoch_noIhb= [epoch_noIhb ii];
            
            tta_amp_noIhbTotal= [tta_amp_noIhbTotal NaN];
            tta_pwr_noIhbTotal= [tta_pwr_noIhbTotal NaN];
            tta_phi_noIhbTotal= [tta_phi_noIhbTotal NaN];
            sg_amp_noIhbTotal= [sg_amp_noIhbTotal NaN];
            sg_pwr_noIhbTotal= [sg_pwr_noIhbTotal NaN];
            sg_phi_noIhbTotal= [sg_phi_noIhbTotal NaN];
            fg_amp_noIhbTotal= [fg_amp_noIhbTotal NaN];
            fg_pwr_noIhbTotal= [fg_pwr_noIhbTotal NaN];
            fg_phi_noIhbTotal= [fg_phi_noIhbTotal NaN];
            epoch_noIhbTotal= [epoch_noIhbTotal ii];
            NIoutTbl=BC_EpochSpdAmpDcmtdTbl(NaN, NaN, NaN,iv_noInhb,ii);
        else
            
            this_csc = restrict(csc, iv_noInhb.tstart(ii), iv_noInhb.tend(ii));
            
            % get the power
            t_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [5 12]);
            sg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [30 58]);
            fg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [60 100]);
            
            ref_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 50]);
            
            t_bp_noinhib(ii)= t_bp;
            sg_bp_noinhib(ii) = sg_bp;
            fg_bp_noinhib(ii) = fg_bp;
            
            t_bp_noinhib_norm(ii)= t_bp/ ref_bp;
            sg_bp_noinhib_norm(ii) = sg_bp/ ref_bp;
            fg_bp_noinhib_norm(ii) = fg_bp/ ref_bp;
            
            NIoutTbl=BC_EpochSpdAmpDcmtdTbl(theta_decimated, SG_decimated,FG_decimated,iv_noInhb,ii);
            
            %Instantaneus speed_pwr_analaysis
            %Create a table with epoch#, condition, amplitude, velocity
            this_theta_amp=restrict(theta_amp_phi, iv_noInhb.tstart(ii), iv_noInhb.tend(ii));
            this_Sg_amp=restrict(SG_amp_phi, iv_noInhb.tstart(ii), iv_noInhb.tend(ii));
            this_Fg_amp=restrict(FG_amp_phi, iv_noInhb.tstart(ii), iv_noInhb.tend(ii));
            
            %Creating a copy of the amplitude vectors so I can use all the values
            this_theta_amp2=this_theta_amp;
            this_Sg_amp2=this_Sg_amp;
            this_Fg_amp2=this_Fg_amp;
            %Combining to a single structure
            this_theta_amp2.data(4:6,:)=this_Sg_amp2.data;
            this_theta_amp2.data(7:9,:)=this_Fg_amp2.data;
            this_theta_amp2.label={'tta_amp','tta_pwr','tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi'};
            %Extracting into vectors for setting a table later
            tta_amp_noIhbTotal= [tta_amp_noIhbTotal this_theta_amp2.data(1,:)];
            tta_pwr_noIhbTotal= [tta_pwr_noIhbTotal this_theta_amp2.data(2,:)];
            tta_phi_noIhbTotal= [tta_phi_noIhbTotal this_theta_amp2.data(3,:)];
            sg_amp_noIhbTotal= [sg_amp_noIhbTotal this_theta_amp2.data(4,:)];
            sg_pwr_noIhbTotal= [sg_pwr_noIhbTotal this_theta_amp2.data(5,:)];
            sg_phi_noIhbTotal= [sg_phi_noIhbTotal this_theta_amp2.data(6,:)];
            fg_amp_noIhbTotal= [fg_amp_noIhbTotal this_theta_amp2.data(7,:)];
            fg_pwr_noIhbTotal= [fg_pwr_noIhbTotal this_theta_amp2.data(8,:)];
            fg_phi_noIhbTotal= [fg_phi_noIhbTotal this_theta_amp2.data(9,:)];
            epoch_noIhbTotal= [epoch_noIhbTotal repmat(ii,1,length(this_theta_amp2.data))];
            
            this_pos= restrict(pos, iv_noInhb.tstart(ii), iv_noInhb.tend(ii));
            this_speed=this_pos;
            this_speed.data=[];this_speed.data=this_pos.data(5,:);this_speed.label={'speed'};this_speed.units={'cm/s'};
            extraction_idx=nearest_idx(this_speed.tvec, this_theta_amp.tvec); %Get the idx of the amplitude more similar to those of the speed
            %Prunning the amplitude to only the time points shared with the velocity vector
            this_theta_amp.data=this_theta_amp.data(1:3,extraction_idx);
            this_Sg_amp.data=this_Sg_amp.data(1:3,extraction_idx);
            this_Fg_amp.data=this_Fg_amp.data(1:3,extraction_idx);
            
            %Combining everything into a single structure
            
            epoch_spd_amp_phi=[];epoch_spd_amp_phi.type='tsd';epoch_spd_amp_phi.tvec=this_speed.tvec;epoch_spd_amp_phi.cfg=this_theta_amp.cfg;
            epoch_spd_amp_phi.data(1,:)=this_speed.data;
            epoch_spd_amp_phi.data(2:4,:)=this_theta_amp.data;
            epoch_spd_amp_phi.data(5:7,:)=this_Sg_amp.data;
            epoch_spd_amp_phi.data(8:10,:)=this_Fg_amp.data;
            epoch_spd_amp_phi.data(11,:)=repmat(ii,1,length(this_speed.tvec));
            epoch_spd_amp_phi.label={'spd', 'tta_amp','tta_pwr','tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi', '#epoch'};
            %Extracting the data from the structure and combining vector for all the epcohs
            spd_noIhb=[spd_noIhb epoch_spd_amp_phi.data(1,:)];
            tta_amp_noIhb= [tta_amp_noIhb epoch_spd_amp_phi.data(2,:)];
            tta_pwr_noIhb= [tta_pwr_noIhb epoch_spd_amp_phi.data(3,:)];
            tta_phi_noIhb= [tta_phi_noIhb epoch_spd_amp_phi.data(4,:)];
            sg_amp_noIhb= [sg_amp_noIhb epoch_spd_amp_phi.data(5,:)];
            sg_pwr_noIhb= [sg_pwr_noIhb epoch_spd_amp_phi.data(6,:)];
            sg_phi_noIhb= [sg_phi_noIhb epoch_spd_amp_phi.data(7,:)];
            fg_amp_noIhb= [fg_amp_noIhb epoch_spd_amp_phi.data(8,:)];
            fg_pwr_noIhb= [fg_pwr_noIhb epoch_spd_amp_phi.data(9,:)];
            fg_phi_noIhb= [fg_phi_noIhb epoch_spd_amp_phi.data(10,:)]; 
            epoch_noIhb= [epoch_noIhb epoch_spd_amp_phi.data(11,:)];
            % Phase mod SG and FG
            this_th = restrict(theta_csc, iv_noInhb.tstart(ii), iv_noInhb.tend(ii));
            this_sg = restrict(SG_csc, iv_noInhb.tstart(ii), iv_noInhb.tend(ii));
            this_fg = restrict(FG_csc, iv_noInhb.tstart(ii), iv_noInhb.tend(ii));
            %SG
            this_SG_phi_amp=BC_phase_amp_norm_bins(this_th,this_sg);
            modidx_SG_noinhib(ii) =MS_ModIdx(this_SG_phi_amp);
            %FG
            this_FG_phi_amp=BC_phase_amp_norm_bins(this_th,this_fg);
            modidx_FG_noinhib(ii) =MS_ModIdx(this_FG_phi_amp);

            % cross freq coupling
            [~, F, ~,P] = spectrogram(this_csc.data,hanning(win_s),win_s/2,1:0.1:160,this_csc.cfg.hdr{1}.SamplingFrequency); % spectrogram -- will take a while to compute!
            [r_noinhib(:,:,ii),~] = corrcoef(10*log10(P')); % correlation matrix (across frequencies) of spectrogram
                        
        end
        NoInhbTbl=[NoInhbTbl; NIoutTbl];
    end
 %putting the speed-pwr into a table so its easier to extract later
    spdPwrTbl_noInhb=table( spd_noIhb',tta_amp_noIhb',tta_pwr_noIhb', tta_phi_noIhb',sg_amp_noIhb',sg_pwr_noIhb',sg_phi_noIhb',fg_amp_noIhb',fg_pwr_noIhb',fg_phi_noIhb',epoch_noIhb', 'VariableNames',{'spd','tta_amp','tta_pwr', 'tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi','epoch'});
 %creating a table for all the variables of the amplitude        
    spdPwrTbl_noInhbTotal= table( tta_amp_noIhbTotal',tta_pwr_noIhbTotal', tta_phi_noIhbTotal',sg_amp_noIhbTotal',sg_pwr_noIhbTotal',sg_phi_noIhbTotal',fg_amp_noIhbTotal',fg_pwr_noIhbTotal',fg_phi_noIhbTotal',epoch_noIhbTotal', 'VariableNames',{'tta_amp','tta_pwr', 'tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi','epoch'});      
                      
    % Running
    spd_running=[];
    tta_amp_running= [];
    tta_pwr_running= [];
    tta_phi_running= [];
    sg_amp_running= [];
    sg_pwr_running= [];
    sg_phi_running= [];
    fg_amp_running= [];
    fg_pwr_running= [];
    fg_phi_running= [];
    epoch_running= [];
    
    tta_amp_runningTotal= [];
    tta_pwr_runningTotal= [];
    tta_phi_runningTotal= [];
    sg_amp_runningTotal= [];
    sg_pwr_runningTotal= [];
    sg_phi_runningTotal= [];
    fg_amp_runningTotal= [];
    fg_pwr_runningTotal= [];
    fg_phi_runningTotal= [];
    epoch_runningTotal= [];

    RunTbl=[];
    for ii = length(iv_running.tstart):-1:1
        if (iv_running.tend(ii) - iv_running.tstart(ii)) < min_trial_dur
            t_bp_run(ii)= NaN;
            sg_bp_run(ii) = NaN;
            fg_bp_run(ii) = NaN;

            modidx_SG_run(ii)= NaN;
            modidx_FG_run(ii)= NaN;
            spd_running=[spd_running NaN];
            tta_amp_running= [tta_amp_running NaN];
            tta_pwr_running= [tta_pwr_running NaN];
            tta_phi_running= [tta_phi_running NaN];
            sg_amp_running= [sg_amp_running NaN];
            sg_pwr_running= [sg_pwr_running NaN];
            sg_phi_running= [sg_phi_running NaN];
            fg_amp_running= [fg_amp_running NaN];
            fg_pwr_running= [fg_pwr_running NaN];
            fg_phi_running= [fg_phi_running NaN]; 
            epoch_running= [epoch_running ii];
            
            tta_amp_runningTotal= [tta_amp_runningTotal NaN];
            tta_pwr_runningTotal= [tta_pwr_runningTotal NaN];
            tta_phi_runningTotal= [tta_phi_runningTotal NaN];
            sg_amp_runningTotal= [sg_amp_runningTotal NaN];
            sg_pwr_runningTotal= [sg_pwr_runningTotal NaN];
            sg_phi_runningTotal= [sg_phi_runningTotal NaN];
            fg_amp_runningTotal= [fg_amp_runningTotal NaN];
            fg_pwr_runningTotal= [fg_pwr_runningTotal NaN];
            fg_phi_runningTotal= [fg_phi_runningTotal NaN]; 
            epoch_runningTotal= [epoch_runningTotal ii];
            RunOutTbl=BC_EpochSpdAmpDcmtdTbl(NaN, NaN, NaN,iv_running,ii);

          

        else
            
            this_csc = restrict(csc, iv_running.tstart(ii), iv_running.tend(ii));
            
            % get the power
            t_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [5 12]);
            sg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [30 58]);
            fg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [60 100]);
            
            ref_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 50]);
            
            t_bp_run_norm(ii)= t_bp/ ref_bp;
            sg_bp_run_norm(ii) = sg_bp/ ref_bp;
            fg_bp_run_norm(ii) = fg_bp/ ref_bp;

            RunOutTbl=BC_EpochSpdAmpDcmtdTbl(theta_decimated, SG_decimated,FG_decimated,iv_running,ii);
            %Instantaneus speed_pwr_analaysis
            %Create a table with epoch#, condition, amplitude, velocity
            this_theta_amp=restrict(theta_amp_phi, iv_running.tstart(ii), iv_running.tend(ii));
            this_Sg_amp=restrict(SG_amp_phi, iv_running.tstart(ii), iv_running.tend(ii));
            this_Fg_amp=restrict(FG_amp_phi, iv_running.tstart(ii), iv_running.tend(ii));
            
            %Creating a copy of the amplitude vectors so I can use all the values
            this_theta_amp2=this_theta_amp;
            this_Sg_amp2=this_Sg_amp;
            this_Fg_amp2=this_Fg_amp;
            %Combining to a single structure
            this_theta_amp2.data(4:6,:)=this_Sg_amp2.data;
            this_theta_amp2.data(7:9,:)=this_Fg_amp2.data;
            this_theta_amp2.label={'tta_amp','tta_pwr','tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi'};
            %Extracting into vectors for setting a table later
            tta_amp_runningTotal= [tta_amp_runningTotal this_theta_amp2.data(1,:)];
            tta_pwr_runningTotal= [tta_pwr_runningTotal this_theta_amp2.data(2,:)];
            tta_phi_runningTotal= [tta_phi_runningTotal this_theta_amp2.data(3,:)];
            sg_amp_runningTotal= [sg_amp_runningTotal this_theta_amp2.data(4,:)];
            sg_pwr_runningTotal= [sg_pwr_runningTotal this_theta_amp2.data(5,:)];
            sg_phi_runningTotal= [sg_phi_runningTotal this_theta_amp2.data(6,:)];
            fg_amp_runningTotal= [fg_amp_runningTotal this_theta_amp2.data(7,:)];
            fg_pwr_runningTotal= [fg_pwr_runningTotal this_theta_amp2.data(8,:)];
            fg_phi_runningTotal= [fg_phi_runningTotal this_theta_amp2.data(9,:)];
            epoch_runningTotal= [epoch_runningTotal repmat(ii,1,length(this_theta_amp2.data))];
            
            this_pos= restrict(pos, iv_running.tstart(ii), iv_running.tend(ii));
            this_speed=this_pos;
            this_speed.data=[];this_speed.data=this_pos.data(5,:);this_speed.label={'speed'};this_speed.units={'cm/s'};
            extraction_idx=nearest_idx(this_speed.tvec, this_theta_amp.tvec); %Get the idx of the amplitude more similar to those of the speed
            %Prunning the amplitude to only the time points shared with the velocity vector
            this_theta_amp.data=this_theta_amp.data(1:3,extraction_idx);
            this_Sg_amp.data=this_Sg_amp.data(1:3,extraction_idx);
            this_Fg_amp.data=this_Fg_amp.data(1:3,extraction_idx);
            
            %Combining everything into a single structure
            
            epoch_spd_amp_phi=[];epoch_spd_amp_phi.type='tsd';epoch_spd_amp_phi.tvec=this_speed.tvec;epoch_spd_amp_phi.cfg=this_theta_amp.cfg;
            epoch_spd_amp_phi.data(1,:)=this_speed.data;
            epoch_spd_amp_phi.data(2:4,:)=this_theta_amp.data;
            epoch_spd_amp_phi.data(5:7,:)=this_Sg_amp.data;
            epoch_spd_amp_phi.data(8:10,:)=this_Fg_amp.data;
            epoch_spd_amp_phi.data(11,:)=repmat(ii,1,length(this_speed.tvec));
            epoch_spd_amp_phi.label={'spd', 'tta_amp','tta_pwr','tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi', '#epoch'};
            %Extracting the data from the structure and combining vector for all the epcohs
            spd_running=[spd_running epoch_spd_amp_phi.data(1,:)];
            tta_amp_running= [tta_amp_running epoch_spd_amp_phi.data(2,:)];
            tta_pwr_running= [tta_pwr_running epoch_spd_amp_phi.data(3,:)];
            tta_phi_running= [tta_phi_running epoch_spd_amp_phi.data(4,:)];
            sg_amp_running= [sg_amp_running epoch_spd_amp_phi.data(5,:)];
            sg_pwr_running= [sg_pwr_running epoch_spd_amp_phi.data(6,:)];
            sg_phi_running= [sg_phi_running epoch_spd_amp_phi.data(7,:)];
            fg_amp_running= [fg_amp_running epoch_spd_amp_phi.data(8,:)];
            fg_pwr_running= [fg_pwr_running epoch_spd_amp_phi.data(9,:)];
            fg_phi_running= [fg_phi_running epoch_spd_amp_phi.data(10,:)]; 
            epoch_running= [epoch_running epoch_spd_amp_phi.data(11,:)];
            % Phase mod SG and FG
            this_th = restrict(theta_csc,iv_running.tstart(ii), iv_running.tend(ii));
            this_sg = restrict(SG_csc, iv_running.tstart(ii), iv_running.tend(ii));
            this_fg = restrict(FG_csc, iv_running.tstart(ii), iv_running.tend(ii));
            %SG
            this_SG_phi_amp=BC_phase_amp_norm_bins(this_th,this_sg);
            modidx_SG_run(ii) =MS_ModIdx(this_SG_phi_amp);
            %FG
            this_FG_phi_amp=BC_phase_amp_norm_bins(this_th,this_fg);
            modidx_FG_run(ii) =MS_ModIdx(this_FG_phi_amp);
          % cross freq coupling
                [~, F, ~,P] = spectrogram(this_csc.data,hanning(win_s),win_s/2,1:0.1:160,this_csc.cfg.hdr{1}.SamplingFrequency); % spectrogram -- will take a while to compute!
       
                [r_run(:,:,ii),~] = corrcoef(10*log10(P')); % correlation matrix (across frequencies) of spectrogram's power
                
        end
        RunTbl=[RunTbl; RunOutTbl];
    end
 %putting the speed-pwr into a table so its easier to extract later
    spdPwrTbl_running=table( spd_running',tta_amp_running',tta_pwr_running', tta_phi_running',sg_amp_running',sg_pwr_running',sg_phi_running',fg_amp_running',fg_pwr_running',fg_phi_running',epoch_running', 'VariableNames',{'spd','tta_amp','tta_pwr', 'tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi','epoch'});
 %creating a table for all the variables of the amplitude        
    spdPwrTbl_runningTotal= table( tta_amp_runningTotal',tta_pwr_runningTotal', tta_phi_runningTotal',sg_amp_runningTotal',sg_pwr_runningTotal',sg_phi_runningTotal',fg_amp_runningTotal',fg_pwr_runningTotal',fg_phi_runningTotal',epoch_runningTotal', 'VariableNames',{'tta_amp','tta_pwr', 'tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi','epoch'});      
                   
    %% plot power and cross freq coupling
    
    if plot_flag
        figure(0007)
        clf
        subplot(3,3,[1 4])
        boxplot([t_bp_inhib, t_bp_noinhib],[zeros(size(t_bp_inhib)), ones(size(t_bp_noinhib))])
        title('Normalized SG ModIdx ArchT');
        ylabel('Z-Scores');
        ax = gca;
        %        ax.
        title('Theta power');
        set(gca, 'XTickLabel', {'Silencing', 'No silencing'})
        
        colors = [BC_color_genertor('powder_blue');  % RGB values for Group 2
                  BC_color_genertor('archt_green')]; % RGB values for Group1;
        h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
        for aa = 1:numel(h)
            patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
        end
        % Customize line and whisker colors
        set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
        % Adjust plot aesthetics
        set(gca, 'TickDir', 'out');  % Move ticks outside the plot
        box off;                     % Turn off the box around the plot
        
        % Slow Gamma
        subplot(3,3,[2 5])
        boxplot([sg_bp_inhib, sg_bp_noinhib],[zeros(size(sg_bp_inhib)), ones(size(sg_bp_noinhib))])
        title('SG power')
        set(gca, 'XTickLabel', {'Silencing', 'No silencing'})
        h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
        for aa = 1:numel(h)
            patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
        end
        % Customize line and whisker colors
        set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
        % Adjust plot aesthetics
        set(gca, 'TickDir', 'out');  % Move ticks outside the plot
        box off;       
        
        % Fast Gamma
        subplot(3,3,[3 6])
        boxplot([fg_bp_inhib, fg_bp_noinhib],[zeros(size(fg_bp_inhib)), ones(size(fg_bp_noinhib))])
        title('FG power')
        set(gca, 'XTickLabel', {'Silencing', 'No silencing'})
        h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
        for aa = 1:numel(h)
            patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
        end
        % Customize line and whisker colors
        set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
        % Adjust plot aesthetics
        set(gca, 'TickDir', 'out');  % Move ticks outside the plot
        box off;       
        subplot(3,3,7)
        imagesc(F,F,nanmean(r_inhib,3));
        caxis([-0.1 1]); axis xy; colorbar; grid on;
        set(gca,'XLim',[0 100],'YLim',[0 100],'FontSize',14,'XTick',0:20:140,'YTick',0:20:140);
        title('Inhibition')
        
        subplot(3,3,8)
        imagesc(F,F,nanmean(r_noinhib,3));
        caxis([-0.1 1]); axis xy; colorbar; grid on;
        set(gca,'XLim',[0 100],'YLim',[0 100],'FontSize',14,'XTick',0:20:140,'YTick',0:20:140);
        title('No Inhibition')
        
        subplot(3,3,9)
        %        this_r = ;
        %        this_r(logical(eye(size(this_r)))) = NaN;
        imagesc(F,F,nanmean(r_run,3));
        caxis([-0.1 .5]); axis xy; colorbar; grid on;
        set(gca,'XLim',[0 100],'YLim',[0 100],'FontSize',14,'XTick',0:20:140,'YTick',0:20:140);
        title('Running')
        
        % Adjust figure properties
        fig = gcf;                   % Get current figure handle
        fig.Color = [1 1 1];         % Set background color to white
        % Resize the figure (optional)
        fig.Position = [200, 200, 1000, 700];  % [x, y, width, height]
        %saveas(gcf, 'FG_NormModIdx.png');%%%% fill in nice plotting %%%%%%%

    end
    
    if save_flag
        saveas(figure(7),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0007' filesep 'Fig07_' info.subject '.png' ]);
        saveas(figure(7),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0007' filesep 'Fig07_' info.subject '.pdf' ]);
        saveas(figure(7),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0007' filesep 'Fig07_' info.subject '.fig' ]);
    end
    %% Plot the CoModulation matrix
    
    % plots are duplicates for some reason. WIP
    
   

    cfg_como.A_step = .5;
    cfg_como.P_step = .5;
    cfg_como.phi_bins = 18;

    [CoMoI, phi_f, amp_f] = MS_phase_freq(cfg_como, csc_inhb, [4 12], [30 100]);
    [CoMoNI, phi_f, amp_f] = MS_phase_freq(cfg_como, csc_noinhb, [4 12], [30 100]);
    [CoMoR, phi_f, amp_f] = MS_phase_freq(cfg_como, csc_running, [4 12], [30 100]);
    if plot_flag
        figure(0008)
        clf
        ax1=subplot(1,2,1); cla;
        %surf(phi_f,amp_f,CoMoI')
        %hold on
        imagesc(phi_f, amp_f, CoMoI');
        
         set(gca, 'ydir', 'normal') %axis('xy')
         clim([0 10^-3.5])
         title('Silencing')
         xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)');
         colorbar('Location', 'southoutside')

        ax2=subplot(1,2,2); cla;
        imagesc(phi_f, amp_f, CoMoNI');
        %surf(CoMoNI)
        set(gca, 'ydir', 'normal')
        caxis([0 10^-3])
        title('No Silencing')
        xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)');
        colorbar('Location', 'southoutside')

        % subplot(1,3,3); cla;
        % imagesc(phi_f, amp_f, CoMoR');
        % set(gca, 'ydir', 'normal')
        % caxis([0 10^-3])
        % title('Running')
        %xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)');
        %colorbar('Location', 'southoutside')

        SetFigure([], gcf)
        maximize

if save_flag
        saveas(figure(8),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0008' filesep 'Fig08_' info.subject '.png' ]);
        saveas(figure(8),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0008' filesep 'Fig08_' info.subject '.pdf' ]);
        saveas(figure(8),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0008' filesep 'Fig08_' info.subject '.fig' ]);
end

    end
       
    %% Wavelet
    if plot_flag
        figure(1009)
        cwt(csc.data, fs);
        [cfs, frq]=cwt(csc.data, fs);
        AX = gca;
        [minf,maxf] = cwtfreqbounds(numel(csc.data),fs);
        
        freq = 2.^(round(log2(minf)):round(log2(maxf)));
        AX.YTickLabelMode = 'auto';
        AX.YTick = freq;
        ylim([1 140])
        iv_inhb_min=BC_iv2min(iv_inhb);
        iv_noInhb_min=BC_iv2min(iv_noInhb);
        
        h01=LTplotIvBarsScalogram(iv_inhb_min,BC_color_genertor('Archt_green'),0.4)
        h02=LTplotIvBarsScalogram(iv_noInhb_min,BC_color_genertor('Web_orange'),0.2)
        
        %Aesthetics
        leg=legend([h01;h02],{'Light','Running No Light'});
        legendFontSize = 12;         % Adjust the font size as needed
        leg.Position=[0.732 0.935 0.09047 0.0487];
        set(leg, 'FontSize', legendFontSize);
        box off;
        legend boxon ;
        set(gca, 'TickDir', 'out');
        fig = gcf;                   % Get current figure handle
        fig.Color = [1 1 1];         % Set background color to white
        fig.Position = [100, 100, 1200, 700];  % [x, y, width, height]
    end
if save_flag
        saveas(figure(9),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0009' filesep 'Fig09_' info.subject '.png' ]);
        %saveas(figure(9),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0009' filesep 'Fig09_' info.subject '.pdf' ]);
        saveas(figure(9),[inter_dir filesep 'AutomaticFigures' filesep 'Fig0009' filesep 'Fig09_' info.subject '.fig' ]);
end
    
%% Figure example LFP

if plot_flag

    fig190=figure(190)
    clf
    
    
    hold on

    %Silencing
    %Theta filtered
    ax1=subplot(5,2,1,'position', [0.05, 0.82, 0.43, 0.16])
    plot(csc.tvec, (theta_csc.data), 'color',BC_color_genertor('Swamp_green') , 'linewidth', 1);
    ylabel({'Theta filtered signal','(mV)'})
    h01=LTplotIvBars(iv_inhb,theta_csc.data,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,theta_csc.data,BC_color_genertor('Burnt_orange'),0.1);
    ylim([min(theta_csc.data) max(theta_csc.data)]);
    xlim([0 max(csc.tvec)]);
    title('Silencing')
    
    %Theta phase
    ax3=subplot(5,2,3,'position', [0.05, 0.63, 0.43, 0.16])
    hold on
    plot(csc.tvec, theta_amp_phi.data(3,:), 'color',BC_color_genertor('Swamp_green') , 'linewidth', 1);
    ylabel({'Theta phase','(rad)'})
    ylim([min(theta_amp_phi.data(3,:))-0.1 max(theta_amp_phi.data(3,:))+0.1]);
    xlim([0 max(csc.tvec)]);
    
    % SG Amplitude 
    ax5=subplot(5,2,5,'position', [0.05, 0.44, 0.43, 0.16]);
    hold on
    plot(csc.tvec, (SG_csc.data), 'color',BC_color_genertor('Web_orange') , 'linewidth', 1);
    plot(csc.tvec, abs(hilbert(SG_csc.data)), 'color', BC_color_genertor('Web_orange'), 'linewidth', .5);
    ylabel({'SG amplitude envelope','(mV)'})
    h01=LTplotIvBars(iv_inhb,SG_csc.data,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,SG_csc.data,BC_color_genertor('Burnt_orange'),0.1);
    ylim([min(SG_csc.data) max(SG_csc.data)]);
    xlim([0 max(csc.tvec)]);

% FG amplitude 
    ax7=subplot(5,2,7,'position', [0.05, 0.25, 0.43, 0.16]);
    hold on
    plot(csc.tvec, (FG_csc.data), 'color', BC_color_genertor('Red_crayola'), 'linewidth', 1);
    plot(csc.tvec, abs(hilbert(FG_csc.data)), 'color', BC_color_genertor('Red_crayola'), 'linewidth', .5);
    ylabel({'FG amplitude envelope','(mV)'})    
    h01=LTplotIvBars(iv_inhb,FG_csc.data,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,FG_csc.data,BC_color_genertor('Burnt_orange'),0.1);
     ylim([min(FG_csc.data) max(FG_csc.data)]);
    xlim([0 max(csc.tvec)]);

% Raw signal 
    ax9=subplot(5,2,9,'position', [0.05, 0.06, 0.43, 0.16]);
    plot(csc.tvec, csc.data+0.0005, 'color', BC_color_genertor('Oxford_blue'), 'linewidth', 1);
    ylabel({'Raw LFP signal', '(mV)'})
    h01=LTplotIvBars(iv_inhb,csc.data+0.0005,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,csc.data+0.0005,BC_color_genertor('Burnt_orange'),0.1);
    ylim([min(csc.data+0.0005) max(csc.data+0.0005)]);
    xlim([0 max(csc.tvec)]);

   
    
    %No silencing
    %Theta filtered
    ax2=subplot(5,2,2,'position', [0.52,0.82, 0.43, 0.16]);
    plot(csc.tvec, (theta_csc.data), 'color',BC_color_genertor('Swamp_green') , 'linewidth', 1);
    h01=LTplotIvBars(iv_inhb,theta_csc.data,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,theta_csc.data,BC_color_genertor('Burnt_orange'),0.1);
    ylim([min(theta_csc.data) max(theta_csc.data)]);
    xlim([0 max(csc.tvec)]);
    title('No silencing')

     %Theta phase
    ax4=subplot(5,2,4,'position', [0.52, 0.63, 0.43, 0.16])
    plot(csc.tvec, theta_amp_phi.data(3,:), 'color',BC_color_genertor('Swamp_green') , 'linewidth', 1);
    
    h01=LTplotIvBars(iv_inhb,theta_csc.data,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,theta_csc.data,BC_color_genertor('Burnt_orange'),0.1);
    ylim([min(theta_amp_phi.data(3,:)) max(theta_amp_phi.data(3,:))]);
    xlim([0 max(csc.tvec)]);

    %SG amplitude
    ax6=subplot(5,2,6,'position', [0.52,0.44, 0.43, 0.16]);
    hold on
    plot(csc.tvec, (SG_csc.data), 'color',BC_color_genertor('Web_orange') , 'linewidth', 1);
    plot(csc.tvec, abs(hilbert(SG_csc.data)), 'color', BC_color_genertor('Web_orange'), 'linewidth', .5);
    h01=LTplotIvBars(iv_inhb,SG_csc.data,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,SG_csc.data,BC_color_genertor('Burnt_orange'),0.1);
    ylim([min(SG_csc.data) max(SG_csc.data)]);
    xlim([0 max(csc.tvec)]);

    %FG amplitude
    ax8=subplot(5,2,8,'position', [0.52,0.25, 0.43, 0.16]);
    hold on
    plot(csc.tvec, (FG_csc.data), 'color', BC_color_genertor('Red_crayola'), 'linewidth', 1);
    plot(csc.tvec, abs(hilbert(FG_csc.data)), 'color', BC_color_genertor('Red_crayola'), 'linewidth', .5);
    h01=LTplotIvBars(iv_inhb,FG_csc.data,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,FG_csc.data,BC_color_genertor('Burnt_orange'),0.1);
    ylim([min(FG_csc.data) max(FG_csc.data)]);
    xlim([0 max(csc.tvec)]);

    %Raw data
    ax10=subplot(5,2,10,'position', [0.52, 0.06, 0.43, 0.16]);
    plot(csc.tvec, csc.data+0.0005, 'color', BC_color_genertor('Oxford_blue'), 'linewidth', 1);
    h01=LTplotIvBars(iv_inhb,csc.data+0.0005,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,csc.data+0.0005,BC_color_genertor('Burnt_orange'),0.1);
    ylim([min(csc.data+0.0005) max(csc.data+0.0005)]);
    xlim([0 max(csc.tvec)]);

    link1=linkprop([ax1,ax2],'Ylim');
    link2=linkprop([ax3,ax4],'Ylim');
    link3=linkprop([ax5,ax6],'Ylim');
    link4=linkprop([ax7,ax8],'Ylim');
    link5=linkprop([ax9,ax10],'Ylim');
    linkaxes([ax1,ax3,ax5,ax7,ax9],'x')
    linkaxes([ax2,ax4,ax6,ax8,ax10],'x')

    % Adjust plot aesthetics
    set(ax1, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax1, 'box','off');
    set(ax2, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax2, 'box','off');
    set(ax3, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax3, 'box','off');
    set(ax4, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax4, 'box','off');
    set(ax5, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax5, 'box','off');
    set(ax6, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax6, 'box','off');
    set(ax7, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax7, 'box','off');
    set(ax8, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax8, 'box','off');
    set(ax9, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax9, 'box','off');
    set(ax10, 'TickDir', 'out');  % Move ticks outside the plot
    set(ax10, 'box','off');
   
end

    %% save outputs
    %inhb
    out.(info.subject).(info.sess).t_bp_inhib = t_bp_inhib; 
    out.(info.subject).(info.sess).sg_bp_inhib = sg_bp_inhib;
    out.(info.subject).(info.sess).fg_bp_inhib = fg_bp_inhib;


    out.(info.subject).(info.sess).t_bp_inhib_norm = t_bp_inhib_norm; 
    out.(info.subject).(info.sess).sg_bp_inhib_norm = sg_bp_inhib_norm;
    out.(info.subject).(info.sess).fg_bp_inhib_norm = fg_bp_inhib_norm;

    out.(info.subject).(info.sess).modidx_SG_inhib = modidx_SG_inhib;
    out.(info.subject).(info.sess).modidx_FG_inhib = modidx_FG_inhib;
    out.(info.subject).(info.sess).z_SGInhb_modidx=z_SGInhb_modidx;
    out.(info.subject).(info.sess).z_FGInhb_modidx=z_FGInhb_modidx;
    out.(info.subject).(info.sess).speedPwrTbl_inhb=spdPwrTbl_inhb;
    out.(info.subject).(info.sess).speedPwrTbl_inhbTotal=spdPwrTbl_inhbTotal;

    out.(info.subject).(info.sess).InhbDcmtSpeedAmpTbl=inhbTbl;

    %noinb
    out.(info.subject).(info.sess).t_bp_noinhib = t_bp_noinhib; 
    out.(info.subject).(info.sess).sg_bp_noinhib = sg_bp_noinhib;
    out.(info.subject).(info.sess).fg_bp_noinhib = fg_bp_noinhib;

    out.(info.subject).(info.sess).t_bp_noinhib_norm = t_bp_noinhib_norm; 
    out.(info.subject).(info.sess).sg_bp_noinhib_norm = sg_bp_noinhib_norm;
    out.(info.subject).(info.sess).fg_bp_noinhib_norm = fg_bp_noinhib_norm;

    out.(info.subject).(info.sess).modidx_SG_noinhib = modidx_SG_noinhib;
    out.(info.subject).(info.sess).modidx_FG_noinhib = modidx_FG_noinhib;
    out.(info.subject).(info.sess).z_SGNoInhb_modidx=z_SGNoInhb_modidx;
    out.(info.subject).(info.sess).z_FGNoInhb_modidx=z_FGNoInhb_modidx;
    out.(info.subject).(info.sess).speedPwrTbl_noInhb=spdPwrTbl_noInhb
    out.(info.subject).(info.sess).speedPwrTbl_noInhbTotal=spdPwrTbl_noInhbTotal;

    out.(info.subject).(info.sess).NoInhbDcmtSpeedAmpTbl=NoInhbTbl;

    %running
    out.(info.subject).(info.sess).t_bp_noinhib = t_bp_noinhib; 
    out.(info.subject).(info.sess).sg_bp_noinhib = sg_bp_noinhib;
    out.(info.subject).(info.sess).fg_bp_noinhib = fg_bp_noinhib;
    out.(info.subject).(info.sess).modidx_SG_noinhib = modidx_SG_noinhib;
    out.(info.subject).(info.sess).modidx_FG_noinhib = modidx_FG_noinhib;
    out.(info.subject).(info.sess).speedPwrTbl_running=spdPwrTbl_running;
    out.(info.subject).(info.sess).speedPwrTbl_runningTotal=spdPwrTbl_runningTotal;

    out.(info.subject).(info.sess).RunningDcmtSpeedAmpTbl=RunTbl;


end
% Formating for saving output

cd(inter_dir)
save('out_eyfp_12_Aug_24.mat','out')
%out_eyfp_07_nov_23=out;save('out_eyfp_07_nov_23.mat','out_eyfp_07_nov_23')