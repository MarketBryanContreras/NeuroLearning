 %% Super_BC_ phase_mod

%% initialize

% Make sure to add code from the CEH2 repo, Vadermeer lab code-shared repo,
% NeuroLearning repo before going forward
%% Dynamic loader
sys=computer;
if contains(sys,'PCWIN')
    % Dynamic based on computer user windows.
    data_dir = [getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\Archt'];
    %data_dir = [getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\eyfp'];
    inter_dir=[getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter']
elseif contains(sys,'MAC')
    % Dynamic based on computer user for mac
    data_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'raw' filesep 'Archt'];
    %data_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'raw' filesep 'eyfp'];
    inter_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'inter'];
else disp("This system is not supported on the dynamic loader yet, pleae add it or contact BC for it")
end
cd(data_dir)
%% Parameters
plot_flag = 1; % switch to 0 if you want to supress verification figures.
time_maze_start = 30;
min_trial_dur = 0.5;
mouse_group=1; %1 for ArchT and 2 for eYFP. This just modify color of the plots

%% Generating color pallet for the mouse group
if mouse_group==1
    gcolors = [BC_color_genertor('powder_blue');  % RGB values for Group 2
        BC_color_genertor('archt_green')]; % RGB values for Group1;
else
    gcolors = [BC_color_genertor('torment_blue');  % RGB values for Group 2
        BC_color_genertor('swamp_green')]; % RGB values for Group1;
end
%% Loop over sessions

% get all sessions with 'D4_INHIBITION'
inhib_dir = dir('*D4_INHIBITION*');
% Find the indices of directories containing 'disconnected'
indices_to_remove = contains({inhib_dir.name}, 'disconnected');

% Remove the directories containing 'disconnected'
inhib_dir(indices_to_remove) = [];

% Loop to load data from raw

for iS =1%:length(inhib_dir)
   
    %% loading
    cd([inhib_dir(iS).folder filesep inhib_dir(iS).name])
    
    parts = strsplit(inhib_dir(iS).name, '_');
    
    info.subject = parts{1};
    info.sess = parts{5};
    info.date = [ parts{2} '_' parts{3} '_' parts{4}];
    
    cfg_csc = [];
    cfg_csc.desired_sampling_frequency = 2000;
    
    if strcmpi(info.subject, 'BC1602')
        cfg_csc.fc ={'CSC2.ncs'}; %'CSC4.ncs'
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    elseif strcmpi(info.subject, 'BC051')
        cfg_csc.fc ={'CSC5.ncs'};
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    elseif strcmpi(info.subject, 'BC1807')
        cfg_csc.fc ={'CSC3.ncs'};
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    elseif strcmpi(info.subject, 'BC053')
        cfg_csc.fc ={'CSC2.ncs'};
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).';
    elseif strcmpi(info.subject, 'BC054')
        cfg_csc.fc ={'CSC5.ncs'};
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).';
    elseif strcmpi(info.subject, 'BC011')
        cfg_csc.fc ={'CSC2.ncs'};
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    elseif strcmpi(info.subject, 'BC013')
        cfg_csc.fc ={'CSC2.ncs'};
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).';
    elseif strcmpi(info.subject, 'BC014')
        cfg_csc.fc ={'CSC2.ncs'};
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).';
    end
    
    [csc, evts, pos] = BC_load_NLX(cfg_csc);
    
    fs = csc.cfg.hdr{1}.SamplingFrequency;
    sSampCsc= time_maze_start*fs;                                          %This is the idx of the samplimg at sec 30
    restrictIvCsc=iv(csc.tvec(sSampCsc),csc.tvec(end));                    %Creates and IV from the time the mouse was placed in the maze to the end of the recording
    csc= restrict(csc,restrictIvCsc);
    % Correction of time
    strt=csc.tvec(1);
    csc.tvec= csc.tvec-strt;                                               %correct time csc
    nn=size(evts.t,2);
    evts.t{1,nn-1}=evts.t{1,nn-1}-strt;                                    % Coreccting laser events start
    evts.t{1,nn}=evts.t{1,nn}-strt;                                        % Coreccting laser events end
    
    %Restiction of time when the mouse is placed in the maze
    sSamp=time_maze_start*30;
    restrictIvPos=iv(pos.tvec(sSamp),pos.tvec(end));
    pos=restrict(pos,restrictIvPos);
    pos.tvec= pos.tvec-pos.tvec(1);
    
    
    iv_inhb = MS_get_evts_off(evts, pattern);
    
    
    %% Trial split
    
    [iv_inhb,iv_noInhb, iv_running] = BC_LT_trialfun(pos, iv_inhb, plot_flag); %This generates the graph of inhibition and velocity while it also retuns the inhibition, no inhibition and runnings epochs 
    
    % control for movement
    pos_spd = pos;
    pos_spd.data = [];
    %Get the x and y position data into a structure
    pos_spd.data(1,:) = pos.data(3,:);
    pos_spd.data(2,:) = pos.data(4,:);
    pos_spd.label = [];
    pos_spd.label = {'x', 'y'}; 

    spd = getLinSpd([],pos_spd);
    cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold =5.0 ;% speed limit in cm/sec
    iv_fast = TSDtoIV(cfg,spd);                                            % only keep intervals with speed above thresh
 
    % Restrict data so it includes fast intervals only
       iv_running = IntersectIV([], iv_fast, iv_running); 
    %% Filter
    % filter the LFP in the theta band
    cfg_filt_t = [];
    cfg_filt_t.type = 'cheby1';                                            %'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [5 12];                                                % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 3;                                                  %type filter order
    cfg_filt_t.display_filter = 0;                                         % use this to see the fvtool
    theta_csc = FilterLFP(cfg_filt_t, csc);                                % filter the raw LFP using
    theta_pow = BC_power(theta_csc);
    theta_phi = angle(hilbert(theta_csc.data));
    
    % Filter the LFP in the slow gamma band
    cfg_filt_t = [];
    cfg_filt_t.type = 'butter';                                            %'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [30 58];                                               % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 4;                                                  %type filter order
    cfg_filt_t.display_filter = 0;                                         % use this to see the fvtool
    SG_csc = FilterLFP(cfg_filt_t, csc);                                   % filter the raw LFP using
    SG_pow = BC_power(SG_csc);
    SG_phi = angle(hilbert(SG_csc.data));
    
    % Filter the LFP in the fast gamma band
    cfg_filt_t = [];
    cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [60 100]; % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 4; %type filter order
    cfg_filt_t.display_filter = 0; % use this to see the fvtool
    FG_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
    FG_pow = BC_power(FG_csc);
    FG_phi = angle(hilbert(FG_csc.data));
    
  
    %% Restricting data to intervals of inhb, noInhb
    %inhb
    csc_inhb = restrict(csc, iv_inhb);
    theta_inhb=restrict(theta_csc, iv_inhb);
    SG_inhb = restrict(SG_csc, iv_inhb);
    FG_inhb = restrict(FG_csc, iv_inhb);
    
    %no-inhb
    csc_noinhb = restrict(csc, iv_noInhb);
    theta_noinhb = restrict(theta_csc, iv_noInhb);
    SG_noinhb = restrict(SG_csc, iv_noInhb);
    FG_noinhb = restrict(FG_csc, iv_noInhb);
    
    %running
    csc_running = restrict(csc, iv_running);
    theta_running = restrict(theta_csc, iv_running);
    SG_running = restrict(SG_csc, iv_running);
    FG_running = restrict(FG_csc, iv_running);
    
    %% Get phase mod over trial

    SGphiAmpNormInhb=BC_phase_amp_norm_bins(theta_inhb,SG_inhb);
    SGInhb_modidx=MS_ModIdx(SGphiAmpNormInhb);
    if plot_flag
        if mouse_group ==1;
        BC_plot_modidx(SGphiAmpNormInhb,BC_color_genertor('Archt_green'),SGInhb_modidx,'SG', 'Light' )
        else
        BC_plot_modidx(SGphiAmpNormInhb,BC_color_genertor('Swamp_green'),SGInhb_modidx,'SG', 'Light' )
        end 
    end
    SGphiAmpNormNoInhb = BC_phase_amp_norm_bins(theta_noinhb,SG_noinhb);
    SGNoInhb_modidx=MS_ModIdx(SGphiAmpNormNoInhb);
    if plot_flag
        if mouse_group ==1;
        BC_plot_modidx(SGphiAmpNormNoInhb,BC_color_genertor('Powder_blue'),SGNoInhb_modidx, 'SG', 'No light')
        else
        BC_plot_modidx(SGphiAmpNormNoInhb,BC_color_genertor('Torment_blue'),SGNoInhb_modidx, 'SG', 'No light')
        end 
    end
    [SGshift_mean, SGshift_std]=LTshifted_meanModIdx(theta_running,SG_running);
%Z scores for SG
    z_SGInhb_modidx = (SGInhb_modidx - SGshift_mean) / SGshift_std;
    z_SGNoInhb_modidx = (SGNoInhb_modidx - SGshift_mean) / SGshift_std;
%Prototype of a bar graph for the mod idx comparison between light and no light
if plot_flag
    figure(19)
    bSG=bar(["No Light" "Light"],[0.3*10^-4  ;0.35*10^-4 ],'FaceColor','flat');
    % Set the face color of each bar individually
    for k = 1:2 % Loop through the number of columns in data
        bSG.CData(k,:) = gcolors(k,:);
    end
    ylabel('Modulation Index');
    set(gca,'fontsize', 14);
    bSG.EdgeColor = 'none';
    set(gca, 'TickDir', 'out');  % Move ticks outside the plot
    set(gca,'box','off');
    fig = gcf;                   % Get current figure handle
    fig.Color = [1 1 1];         % Set background color to white
    fig.Position = [100, 100, 600, 800];  % [x, y, width, height]
end

    %FG
    FGphiAmpNormInhb=BC_phase_amp_norm_bins(theta_inhb,FG_inhb);
    FGInhb_modidx=MS_ModIdx(FGphiAmpNormInhb);
    if plot_flag
        if mouse_group ==1;
        BC_plot_modidx(FGphiAmpNormInhb,BC_color_genertor('Archt_green'),FGInhb_modidx,'FG', 'Light' )
        else
        BC_plot_modidx(FGphiAmpNormInhb,BC_color_genertor('Swamp_green'),FGInhb_modidx,'FG', 'Light' )
        end
    end
    %NoInhb
    FGphiAmpNormNoInhb = BC_phase_amp_norm_bins(theta_noinhb,FG_noinhb);
    FGNoInhb_modidx=MS_ModIdx(FGphiAmpNormNoInhb);
    if plot_flag
        if mouse_group ==1;
        BC_plot_modidx(FGphiAmpNormNoInhb,BC_color_genertor('Powder_blue'),FGNoInhb_modidx, 'FG', 'No light')
        else
        BC_plot_modidx(FGphiAmpNormNoInhb,BC_color_genertor('Torment_blue'),FGNoInhb_modidx, 'FG', 'No light')
        end
    end
    %Shifted
    [FGshift_mean,FGshift_std]=LTshifted_meanModIdx(theta_running,FG_running);
    %Z scores for FG

    z_FGInhb_modidx = (FGInhb_modidx - FGshift_mean) / FGshift_std; 
    z_FGNoInhb_modidx = (FGNoInhb_modidx - FGshift_mean) / FGshift_std;  
    %% Epochs by epochs analysis
    
    % inhibition
    win_s = 256;
            r_inhib = []; 
            r_noinhib = []; 
            r_run = []; 

    for ii = length(iv_inhb.tstart):-1:1
        if (iv_inhb.tend(ii) - iv_inhb.tstart(ii)) < min_trial_dur
            t_bp_inhib(ii)= NaN;
            sg_bp_inhib(ii) = NaN;
            fg_bp_inhib(ii) = NaN;
            modidx_SG_inhib(ii)= NaN;
            modidx_FG_inhib(ii)= NaN;
            r_inhib = NaN(length(1:0.1:160), length(1:0.1:160)); 
        else
            this_csc = restrict(csc, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            
            % get the power
            t_bp =  bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [5 12]);
            sg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [30 58]);
            fg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [60 100]);
            
            ref_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 50]);
            
            t_bp_inhib(ii)= t_bp/ ref_bp;
            sg_bp_inhib(ii) = sg_bp/ ref_bp;
            fg_bp_inhib(ii) = fg_bp/ ref_bp;
            
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
    end

    % NO inhibition
    for ii = length(iv_noInhb.tstart):-1:1
        if (iv_noInhb.tend(ii) - iv_noInhb.tstart(ii)) < min_trial_dur
            t_bp_noinhib(ii)= NaN;
            sg_bp_noinhib(ii) = NaN;
            fg_bp_noinhib(ii) = NaN;
            modidx_SG_noinhib(ii)= NaN;
            modidx_FG_noinhib(ii)= NaN;
            r_noinhib = NaN(length(1:0.1:160), length(1:0.1:160));
        else
            
            this_csc = restrict(csc, iv_noInhb.tstart(ii), iv_noInhb.tend(ii));
            
            % get the power
            t_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [5 12]);
            sg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [30 58]);
            fg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [60 100]);
            
            ref_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 50]);
            
            t_bp_noinhib(ii)= t_bp/ ref_bp;
            sg_bp_noinhib(ii) = sg_bp/ ref_bp;
            fg_bp_noinhib(ii) = fg_bp/ ref_bp;
            
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
    end
 
    % Running
    for ii = length(iv_running.tstart):-1:1
        if (iv_running.tend(ii) - iv_running.tstart(ii)) < min_trial_dur
            t_bp_run(ii)= NaN;
            sg_bp_run(ii) = NaN;
            fg_bp_run(ii) = NaN;

            modidx_SG_run(ii)= NaN;
            modidx_FG_run(ii)= NaN;

        else
            
            this_csc = restrict(csc, iv_running.tstart(ii), iv_running.tend(ii));
            
            % get the power
            t_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [5 12]);
            sg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [30 58]);
            fg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [60 100]);
            
            ref_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 50]);
            
            t_bp_run(ii)= t_bp/ ref_bp;
            sg_bp_run(ii) = sg_bp/ ref_bp;
            fg_bp_run(ii) = fg_bp/ ref_bp;
            
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
    end
   
    %% plot power and cross freq coupling
    
    if plot_flag
        figure(101)
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
    
    %% Plot the CoModulation matrix
    
    % plots are duplicates for some reason. WIP
    
    if plot_flag
        
        cfg_como.A_step = .5;
        cfg_como.P_step = .5;
        cfg_como.phi_bins = 18;
        
        [CoMoI, phi_f, amp_f] = MS_phase_freq(cfg_como, csc_inhb, [4 12], [30 100]);
        [CoMoNI, phi_f, amp_f] = MS_phase_freq(cfg_como, csc_noinhb, [4 12], [30 100]);
        [CoMoR, phi_f, amp_f] = MS_phase_freq(cfg_como, csc_running, [4 12], [30 100]);
        
%%
        figure(202)
        clf
        subplot(1,2,1); cla;
        surf(phi_f,amp_f,CoMoI')
        hold on
        imagesc(phi_f, amp_f, CoMoI');
        zlim([-0.5e-3 2e-3])
        % 
        % set(gca, 'ydir', 'normal') %axis('xy')
        % clim([0 10^-3])
        % title('Silencing')
        % xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)');
        % colorbar('Location', 'southoutside')

        subplot(1,2,2); cla; 
        %imagesc(phi_f, amp_f, CoMoNI');
        surf(CoMoNI)
        % set(gca, 'ydir', 'normal')
        % caxis([0 10^-3])
        % title('No Silencing')
        % xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)'); 
        %                 colorbar('Location', 'southoutside')

        % subplot(1,3,3); cla; 
        % imagesc(phi_f, amp_f, CoMoR');
        % set(gca, 'ydir', 'normal')
        % caxis([0 10^-3])
        % title('Running')
        %         xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)'); 
        %         colorbar('Location', 'southoutside')

                SetFigure([], gcf)
                maximize
        
    end
       
    %% Wavelet
    if plot_flag
        figure(109)
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
        
        h01=LTplotIvBarsScalogram(iv_inhb_min,BC_color_genertor('Archt_green'),0.2)
        h02=LTplotIvBarsScalogram(iv_noInhb_min,BC_color_genertor('Web_orange'),0.2)
        
        %Aesthetics
        leg=legend([h01;h02],{'Silencing','Running no silencing'});
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
    
    %% trial averaged spectra
    
% 
%     Triggered_Spec_FT(csc, iv_inhb.tstart, 'Opto On', [1:0.2:120], [-2 0], [-5 10])
% 
%     figure
%         Triggered_Spec_FT(csc, iv_noInhb.tstart, 'Opto Off', [1:0.2:120], [-2 0], [-5 10])
%         

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
    plot(csc.tvec, theta_phi, 'color',BC_color_genertor('Swamp_green') , 'linewidth', 1);
    ylabel({'Theta phase','(rad)'})
    ylim([min(theta_phi)-0.1 max(theta_phi)+0.1]);
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
    plot(csc.tvec, theta_phi, 'color',BC_color_genertor('Swamp_green') , 'linewidth', 1);
    
    h01=LTplotIvBars(iv_inhb,theta_csc.data,BC_color_genertor('Archt_green'),0.1);
    h02=LTplotIvBars(iv_noInhb,theta_csc.data,BC_color_genertor('Burnt_orange'),0.1);
    ylim([min(theta_phi) max(theta_phi)]);
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
    out.(info.subject).(info.sess).modidx_SG_inhib = modidx_SG_inhib;
    out.(info.subject).(info.sess).modidx_FG_inhib = modidx_FG_inhib;
    out.(info.subject).(info.sess).z_SGInhb_modidx=z_SGInhb_modidx;
    out.(info.subject).(info.sess).z_FGInhb_modidx=z_FGInhb_modidx;
    %noinb
    out.(info.subject).(info.sess).t_bp_noinhib = t_bp_noinhib; 
    out.(info.subject).(info.sess).sg_bp_noinhib = sg_bp_noinhib;
    out.(info.subject).(info.sess).fg_bp_noinhib = fg_bp_noinhib;
    out.(info.subject).(info.sess).modidx_SG_noinhib = modidx_SG_noinhib;
    out.(info.subject).(info.sess).modidx_FG_noinhib = modidx_FG_noinhib;
    out.(info.subject).(info.sess).z_SGNoInhb_modidx=z_SGNoInhb_modidx;
    out.(info.subject).(info.sess).z_FGNoInhb_modidx=z_FGNoInhb_modidx;
    %running
    out.(info.subject).(info.sess).t_bp_noinhib = t_bp_noinhib; 
    out.(info.subject).(info.sess).sg_bp_noinhib = sg_bp_noinhib;
    out.(info.subject).(info.sess).fg_bp_noinhib = fg_bp_noinhib;
    out.(info.subject).(info.sess).modidx_SG_noinhib = modidx_SG_noinhib;
    out.(info.subject).(info.sess).modidx_FG_noinhib = modidx_FG_noinhib;

end
% Formating for saving output

%cd(inter_dir)
%out_Archt_07_nov_23=out;save('out_Archt_07_nov_23.mat','out_Archt_07_nov_23')
%out_eyfp_07_nov_23=out;save('out_eyfp_07_nov_23.mat','out_eyfp_07_nov_23')