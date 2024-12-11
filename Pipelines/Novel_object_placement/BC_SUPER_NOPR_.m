%% Sandbox for NOPR analsys

%% Initialize

%% I have to create an automatic loader 
% Go to the directory of some data
%BC053 D1
%data_dir= '/Users/bryancontrerasmercado/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/CHRNA2_NOVEL_OBJECT/raw_data/NOPR/BC053_2023_11_16_D1_HAB_T2';
%BC053 D2
%data_dir= '/Users/bryancontrerasmercado/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/CHRNA2_NOVEL_OBJECT/raw_data/NOPR/BC053_2023_11_16_D1_HAB_T2';
% data_dir='/Users/bryancontrerasmercado/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/CHRNA2_NOVEL_OBJECT/raw_data/Behavior/BC1807_2023_07_07_D2_NOPR';
%mac
%data_dir= 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_NOVEL_OBJECT\raw_data\NOPR\BC053_2023_11_16_D1_HAB_T2';
%data_dir= 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_NOVEL_OBJECT\raw_data\Behavior\BC1807_2023_07_07_D2_NOPR'
%'windows
     
% cd(data_dir)
%% General parameters
plot_flag = 00;
video_flag=00;
save_output=01;
save_figures=00;
%% Rolling through the folders with dynamic loader
sys=computer;
if contains(sys,'PCWIN')
    % Dynamic based on computer user windows.
    data_dir = [getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_NOVEL_OBJECT\raw_data\Behavior'];
    %data_dir = [getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\eyfp'];
    inter_dir=[getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_NOVEL_OBJECT\inter'];
elseif contains(sys,'MAC')
    % Dynamic based on computer user for mac
    data_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_NOVEL_OBJECT' filesep 'raw_data' filesep 'Behavior'];
    %data_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'raw' filesep 'eyfp'];
    inter_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_NOVEL_OBJECT' filesep 'inter'];
else disp("This system is not supported on the dynamic loader yet, pleae add it or contact BC for it")
end
cd(data_dir)

% get the folders from D1 and D2
D1_dir = dir('*D1_hab*');
D2_dir=dir('*D2_nopr*');
%% Initializing outputs
    out=[];
    sleep_mod_idx=[];

%% Loop to load data from raw

for iS=1:length(D1_dir)
    %% Colloecting subject info
    %subjectFolder = fullfile(data_dir, inhib_dir(iS).name);

    cd([D1_dir(iS).folder filesep D1_dir(iS).name])
    parts = pwd;
    infoD1= infoSession(parts);
    %% Individual parameters
    if infoD1.subject=="BC011";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC3.ncs';lfp_chanD2='CSC3.ncs';%6...5,4
    elseif infoD1.subject=="BC1807";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC5.ncs';lfp_chanD2='CSC5.ncs';%6...7...4.3
    elseif infoD1.subject=="BC054";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC7.ncs';lfp_chanD2='CSC7.ncs';%7...6,5,4,3,
    elseif infoD1.subject=="BC053";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC7.ncs';lfp_chanD2='CSC7.ncs';
    elseif infoD1.subject=="BC051" ; %Adjudted to match cable number, due to mouse bitting cable during end of recording in D1 and changing cables in D2 and tracked the corresponding cable to match corresponding cable
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC6.ncs';lfp_chanD2='CSC2.ncs';%2...3,4,5,6,7
    elseif infoD1.subject=="BC014";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC4.ncs';lfp_chanD2='CSC4.ncs';
    elseif infoD1.subject=="BC013";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC7.ncs';lfp_chanD2='CSC7.ncs';
    end
    %% Loading data

    % Load the CSC guide
    cfg = [];
    cfg.fc = {lfp_chan};
    cscD1 = MS_LoadCSC(cfg);

    % Load the EMG
    cfg_emg = [];
    cfg_emg.fc = {emg_chan};
    emgD1= MS_LoadCSC(cfg_emg);

    % load the events
    evtsD1 = LoadEvents([]);

    %% Restrict the data only to the sleep phase
    [csc_sD1, emg_sD1]=Sleep_restrict(evtsD1,cscD1,emgD1,infoD1);
 
    %% plot a bit of data for quality check end verify the awake states
    if plot_flag

        figure(4)
        clf

        ax(1) = subplot(2,1,1);
        plot(csc_sD1.tvec, csc_sD1.data);
        legend('HC LFP')

        ax(2) = subplot(2,1,2);
        plot(emg_sD1.tvec, emg_sD1.data - .001, 'r') % offset it a bit
        legend('emg')

        linkaxes(ax, 'x'); % locks the x axes of the subplots so if you zoom on one the other zooms in.

        % fix the annoying wide limits on the x axis
        xlim([csc_sD1.tvec(1) csc_sD1.tvec(end)])
    end

    %% sleep state *adapt

    % specify known wake times or periods to ignore. Note; The perios are in
    % the interval format, so you want to write when do they start and when do
    % they end in pairs
    
        if strcmpi(infoD1.subject,"BC1807")
            wake_tD1 = [0 3292 3657 5083 8332 10042 12118 13092];
        elseif strcmpi(infoD1.subject,"BC053") %Specify the wake period according to subject and session
            wake_tD1 = [0 3093 5015 5631 8824 9340 10605 10880 11531 11803 12450 13020 13618 13743 14188 14397];
        elseif strcmpi(infoD1.subject,"BC011")
            wake_tD1 = [0 1493 3455 4407 10405 13405 14654 14710];
        elseif strcmpi(infoD1.subject,"BC013")
            wake_tD1 = [0 4296 4735 5514 9823 10377 12838 13281];
        elseif strcmpi(infoD1.subject,"BC014")
            wake_tD1 = [0 4869 6952 9715 11557 11861 12753 13105 13301 13531];
        elseif strcmpi(infoD1.subject,"BC051") % check this animal, apparenlty there was only 2 hrs of sleep 
            wake_tD1 = [0 954 1542 1603 2770 2822 3081 3144 4588 4661 6238 8179];
        elseif strcmpi(infoD1.subject,"BC054") 
            wake_tD1 = [0 2054 2710 3859 4111 4621 4971 5234 6148 6328 6695 7153 7498 8198 8410 8725 9564 9999 11770 13106 13683 14402];
        end
        %To add to part2
        % if strcmpi(infoD2.subject,"BC1807")
        %     wake_t = [0 4735 7179 9298 10043 11642 13919 14403];
        % elseif strcmpi(infoD2.subject,"BC011")
        %     wake_t = [0 9199 10295  11070];
        % elseif strcmpi(infoD2.subject,"BC013")
        %     wake_t = [0 688 1817 3130 4335 4783 6820 7511 10772 11917];
        % elseif strcmpi(infoD2.subject,"BC014")
        %     wake_t = [0 3730 ];
        % elseif strcmpi(infoD2.subject,"BC051") % check this animal, apparenlty there was only 2 hrs of sleep
        %     wake_t = [0 2453 2944 4320 4721 5119 6653 7569 10306 11233 13043 13639];
        % elseif strcmpi(infoD2.subject,"BC053") %Specify the wake period according to subject and session
        %     %wake_t = [0 3093 5015 5631 8824 9340 10605 10880 11531 11803 12450 13020 13618 13743 14188 14397];
        % elseif strcmpi(infoD2.subject,"BC054") % check this animal, apparenlty there was only 2 hrs of sleep
        %     wake_t = [0 2518 4597 5331 6521 8947 9211 7153 9607 10094 12110 12942];
        % end
        %end

    wake_idxD1 = nearest_idx(wake_tD1, csc_sD1.tvec); %Converts time to correspondant sample idx
    wake_idxD1 = reshape(wake_idxD1,2, length(wake_idxD1)/2)'; %reshape(columns, rows)
    %% score the sleep
    [hypnoD1, csc_outD1, emg_outD1] = dSub_Sleep_screener(0, csc_sD1, emg_sD1, wake_idxD1);  % can add in 'wake_idx' as the last input.
    %% Obtaining the time stamps for these sleep states
    [iv_awakeD1, iv_swsD1, iv_remD1]= BC_sleep_iv_extractor(hypnoD1);
    %% Obtaining IV for inhibition in case this is session D2
 
     %% Filter CSC in the theta and gamma bands
    % filter the LFP in the theta band
    cfg_filt_t = [];
    cfg_filt_t.type = 'cheby1';                                            %'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [5 12];                                                % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 3;                                                  %type filter order
    cfg_filt_t.display_filter = 0;                                         % use this to see the fvtool
    theta_cscD1 = FilterLFP(cfg_filt_t, csc_outD1);                                % filter the raw LFP using
    
    % Filter the LFP in the slow gamma band
    cfg_filt_sg = [];
    cfg_filt_sg.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_sg.f  = [30 58]; % freq range to match Mizuseki et al. 2011
    cfg_filt_sg.order = 4; %type filter orderf
    cfg_filt_sg.display_filter = 0; % use this to see the fvtool
    SG_cscD1 = FilterLFP(cfg_filt_sg, csc_outD1); % filter the raw LFP using

    % Filter the LFP in the fast gamma band
    cfg_filt_fg = [];
    cfg_filt_fg.type = 'butter';                                            %'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_fg.f  = [60 100];                                               % freq range to match Mizuseki et al. 2011
    cfg_filt_fg.order = 4;                                                  %type filter order
    cfg_filt_fg.display_filter = 0;                                         % use this to see the fvtool
    FG_cscD1 = FilterLFP(cfg_filt_fg, csc_outD1);                                   % filter the raw LFP using

    %% Band power extraction per REM episode
    min_trial_dur=10;
    win_s = 256;
    %loop over the REM IV
    for ii =length(iv_remD1.tstart):-1:1
        %if (iv_rem.tend(thisIV) - iv_rem.tstart(thisIV)) < min_trial_dur %Testing that the epoch last more than the treshold

            %Extract the band power and mod idx for each rem episode
            this_csc = restrict(csc_sD1, iv_remD1.tstart(ii), iv_remD1.tend(ii));

            t_bpD1 =  bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [4 12]); %Calculates the theta band power
            sg_bpD1 = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [30 58]);%Calculates the sg band power
            fg_bpD1 = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [60 100]);%Calculates the fg band power

            ref_bpD1 = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 50]);%Calculates the ref band power
            %Store the theta, sg, and fg power in their varibales
            t_bp_inhibD1(ii)= t_bpD1;
            sg_bp_inhibD1(ii) = sg_bpD1;
            fg_bp_inhibD1(ii) = fg_bpD1;
            
            t_bp_inhib_normD1(ii)= t_bpD1/ ref_bpD1;
            sg_bp_inhib_normD1(ii) = sg_bpD1/ ref_bpD1;
            fg_bp_inhib_normD1(ii) = fg_bpD1/ ref_bpD1;

            %Mod-idx analysis

            this_thD1 = restrict(theta_cscD1, iv_remD1.tstart(ii), iv_remD1.tend(ii));
            this_sgD1 = restrict(SG_cscD1, iv_remD1.tstart(ii), iv_remD1.tend(ii));
            this_fgD1 = restrict(FG_cscD1, iv_remD1.tstart(ii), iv_remD1.tend(ii));
            %SG
            this_SG_phi_ampD1=BC_phase_amp_norm_bins(this_thD1,this_sgD1);
            modidx_SG_inhibD1(ii) =MS_ModIdx(this_SG_phi_ampD1);
            %FG
            this_FG_phi_ampD1=BC_phase_amp_norm_bins(this_thD1,this_fgD1);
            modidx_FG_inhibD1(ii) =MS_ModIdx(this_FG_phi_ampD1);
            % cross freq coupling
            [~, F, ~,P] = spectrogram(this_csc.data,hanning(win_s),win_s/2,1:0.1:160,this_csc.cfg.hdr{1}.SamplingFrequency); % spectrogram -- will take a while to compute!, add 'yaxis for visualization porpuses
            %F: Frequency axis of the spectrogram. P: Power spectral density (or magnitude squared of the Fourier transform) of the signal.
            [r_remD1(:,:,ii),~] = corrcoef(10*log10(P')); % correlation matrix (across frequ`encies) of spectrogram
        %end
        
        % Putting these values in a structure 
        out.(infoD1.subject).(infoD1.session).t_bp_inhib = t_bp_inhibD1;
        out.(infoD1.subject).(infoD1.session).sg_bp_inhib = sg_bp_inhibD1;
        out.(infoD1.subject).(infoD1.session).fg_bp_inhib = fg_bp_inhibD1;

        out.(infoD1.subject).(infoD1.session).t_bp_inhib_norm = t_bp_inhib_normD1;
        out.(infoD1.subject).(infoD1.session).sg_bp_inhib_norm = sg_bp_inhib_normD1;
        out.(infoD1.subject).(infoD1.session).fg_bp_inhib_norm = fg_bp_inhib_normD1;

        out.(infoD1.subject).(infoD1.session).modidx_SG_inhib = modidx_SG_inhibD1;
        out.(infoD1.subject).(infoD1.session).modidx_FG_inhib = modidx_FG_inhibD1;
        % out.(infoD1.subject).(info.sess).z_SGInhb_modidx=z_SGInhb_modidxD1;
        % out.(infoD1.subject).(info.sess).z_FGInhb_modidx=z_FGInhb_modidxD1;
    end
        %% Going to D2 session
    cd([D2_dir(iS).folder filesep D2_dir(iS).name])
    parts = pwd;
    infoD2=infoSession(parts);
    
    %% Loading data

    % Load the CSC guide
    cfg = [];
    cfg.fc = {lfp_chanD2};
    cscD2 = MS_LoadCSC(cfg);

    % Load the EMG
    cfg_emg = [];
    cfg_emg.fc = {emg_chan};
    emgD2= MS_LoadCSC(cfg_emg);

    % load the events
    evtsD2 = LoadEvents([]);

    %% Restrict the data only to the sleep phase
    [csc_sD2, emg_sD2]=Sleep_restrict(evtsD2,cscD2,emgD2,infoD2);
 
    %% plot a bit of data for quality check end verify the awake states
    if plot_flag

        figure(4)
        clf

        ax(1) = subplot(2,1,1);
        plot(csc_sD2.tvec, csc_sD2.data);
        legend('HC LFP')

        ax(2) = subplot(2,1,2);
        plot(emg_sD2.tvec, emg_sD2.data - .001, 'r') % offset it a bit
        legend('emg')

        linkaxes(ax, 'x'); % locks the x axes of the subplots so if you zoom on one the other zooms in.

        % fix the annoying wide limits on the x axis
        xlim([csc_sD2.tvec(1) csc_sD2.tvec(end)])
    end

    %% sleep state
    if strcmpi(infoD2.subject,"BC1807")
        wake_tD2 = [0 4735 7179 9298 10043 11642 13919 14403];
    elseif strcmpi(infoD2.subject,"BC011")
        wake_tD2 = [0 9199 10295  11070];
    elseif strcmpi(infoD2.subject,"BC013")
        wake_tD2 = [0 688 1817 3130 4335 4783 6820 7511 10772 11917];
    elseif strcmpi(infoD2.subject,"BC014")
        wake_tD2 = [0 3730 ];
    elseif strcmpi(infoD2.subject,"BC051") % check this animal, apparenlty there was only 2 hrs of sleep
        wake_tD2 = [0 2453 2944 4320 4721 5119 6653 7569 10306 11233 13043 13639];
    elseif strcmpi(infoD2.subject,"BC053") %Specify the wake period according to subject and session
        %wake_tD2 = [0 3093 5015 5631 8824 9340 10605 10880 11531 11803 12450 13020 13618 13743 14188 14397];
    elseif strcmpi(infoD2.subject,"BC054") % check this animal, apparenlty there was only 2 hrs of sleep
        wake_tD2 = [0 2518 4597 5331 6521 8947 9211 7153 9607 10094 12110 12942];
    end

    wake_idxD2 = nearest_idx(wake_tD2, csc_sD2.tvec); %Converts time to correspondant sample idx
    wake_idxD2 = reshape(wake_idxD2,2, length(wake_idxD2)/2)'; %reshape(columns, rows)
    %% score the sleep
    [hypnoD2, csc_outD2, emg_outD2] = dSub_Sleep_screener(0, csc_sD2, emg_sD2, wake_idxD2);  % can add in 'wake_idx' as the last input.
    %% Obtaining the time stamps for these sleep states
    [iv_awakeD2, iv_swsD2, iv_remD2]= BC_sleep_iv_extractor(hypnoD2);
    %% Obtaining IV for inhibition in case this is session D2
 
     %% Filter CSC in the theta and gamma bands
    % filter the LFP in the theta band
    cfg_filt_t = [];
    cfg_filt_t.type = 'cheby1';                                            %'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [5 12];                                                % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 3;                                                  %type filter order
    cfg_filt_t.display_filter = 0;                                         % use this to see the fvtool
    theta_cscD2 = FilterLFP(cfg_filt_t, csc_outD2);                                % filter the raw LFP using
    
    % Filter the LFP in the slow gamma band
    cfg_filt_sg = [];
    cfg_filt_sg.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_sg.f  = [30 58]; % freq range to match Mizuseki et al. 2011
    cfg_filt_sg.order = 4; %type filter orderf
    cfg_filt_sg.display_filter = 0; % use this to see the fvtool
    SG_cscD2 = FilterLFP(cfg_filt_sg, csc_outD2); % filter the raw LFP using

    % Filter the LFP in the fast gamma band
    cfg_filt_fg = [];
    cfg_filt_fg.type = 'butter';                                            %'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_fg.f  = [60 100];                                               % freq range to match Mizuseki et al. 2011
    cfg_filt_fg.order = 4;                                                  %type filter order
    cfg_filt_fg.display_filter = 0;                                         % use this to see the fvtool
    FG_cscD2 = FilterLFP(cfg_filt_fg, csc_outD2);                                   % filter the raw LFP using

    %% Band power extraction per REM episode
    min_trial_dur=10;
    win_s = 256;
    %loop over the REM IV
    for ii =length(iv_remD2.tstart):-1:1
        %if (iv_rem.tend(thisIV) - iv_rem.tstart(thisIV)) < min_trial_dur %Testing that the epoch last more than the treshold

            %Extract the band power and mod idx for each rem episode
            this_csc = restrict(csc_sD2, iv_remD2.tstart(ii), iv_remD2.tend(ii));

            t_bpD2 =  bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [4 12]); %Calculates the theta band power
            sg_bpD2 = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [30 58]);%Calculates the sg band power
            fg_bpD2 = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [60 100]);%Calculates the fg band power

            ref_bpD2 = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 50]);%Calculates the ref band power
            %Store the theta, sg, and fg power in their varibales
            t_bp_inhibD2(ii)= t_bpD2;
            sg_bp_inhibD2(ii) = sg_bpD2;
            fg_bp_inhibD2(ii) = fg_bpD2;
            
            t_bp_inhib_normD2(ii)= t_bpD2/ ref_bpD2;
            sg_bp_inhib_normD2(ii) = sg_bpD2/ ref_bpD2;
            fg_bp_inhib_normD2(ii) = fg_bpD2/ ref_bpD2;

            %Mod-idx analysis

            this_thD2 = restrict(theta_cscD2, iv_remD2.tstart(ii), iv_remD2.tend(ii));
            this_sgD2 = restrict(SG_cscD2, iv_remD2.tstart(ii), iv_remD2.tend(ii));
            this_fgD2 = restrict(FG_cscD2, iv_remD2.tstart(ii), iv_remD2.tend(ii));
            %SG
            this_SG_phi_ampD2=BC_phase_amp_norm_bins(this_thD2,this_sgD2);
            modidx_SG_inhibD2(ii) =MS_ModIdx(this_SG_phi_ampD2);
            %FG
            this_FG_phi_ampD2=BC_phase_amp_norm_bins(this_thD2,this_fgD2);
            modidx_FG_inhibD2(ii) =MS_ModIdx(this_FG_phi_ampD2);
            % cross freq coupling
            [~, F, ~,P] = spectrogram(this_csc.data,hanning(win_s),win_s/2,1:0.1:160,this_csc.cfg.hdr{1}.SamplingFrequency); % spectrogram -- will take a while to compute!, add 'yaxis for visualization porpuses
            %F: Frequency axis of the spectrogram. P: Power spectral density (or magnitude squared of the Fourier transform) of the signal.
            [r_remD2(:,:,ii),~] = corrcoef(10*log10(P')); % correlation matrix (across frequ`encies) of spectrogram
        %end
        
        % Putting these values in a structure 
        out.(infoD2.subject).(infoD2.session).t_bp_inhib = t_bp_inhibD2;
        out.(infoD2.subject).(infoD2.session).sg_bp_inhib = sg_bp_inhibD2;
        out.(infoD2.subject).(infoD2.session).fg_bp_inhib = fg_bp_inhibD2;

        out.(infoD2.subject).(infoD2.session).t_bp_inhib_norm = t_bp_inhib_normD2;
        out.(infoD2.subject).(infoD2.session).sg_bp_inhib_norm = sg_bp_inhib_normD2;
        out.(infoD2.subject).(infoD2.session).fg_bp_inhib_norm = fg_bp_inhib_normD2;

        out.(infoD2.subject).(infoD2.session).modidx_SG_inhib = modidx_SG_inhibD2;
        out.(infoD2.subject).(infoD2.session).modidx_FG_inhib = modidx_FG_inhibD2;
        % out.(infoD1.subject).(info.sess).z_SGInhb_modidx=z_SGInhb_modidxD1;
        % out.(infoD1.subject).(info.sess).z_FGInhb_modidx=z_FGInhb_modidxD1;
    end

    %% Making a collection fx to plot 
    %Creating a folder for the subject and
    MainFldr='/Users/bryancontrerasmercado/Documents/01_Projects/02_OLM_silencing/[03]Figures/07_REM_distribution';
    fldrnames=dir(MainFldr);
    fldrnames=fldrnames(~ismember({fldrnames.name},{'.','..','.DS_Store'}));
    SfolderPath=[];
    % MIdxFldr='/Users/bryancontrerasmercado/Documents/01_Projects/02_OLM_silencing/[03]Figures/07_REM_distribution/Mod_idx';
    % BPNormFldr= '/Users/bryancontrerasmercado/Documents/01_Projects/02_OLM_silencing/[03]Figures/07_REM_distribution/Normalized_Band_powers';
    % BPRawFldr='/Users/bryancontrerasmercado/Documents/01_Projects/02_OLM_silencing/[03]Figures/07_REM_distribution/Raw_Band_Powers';
    % folders={'MIdxFldr','BPNormFldr','BPRawFldr'};
    for iD=3:-1:1
        folderpath=fullfile(MainFldr,(fldrnames(iD).name));
        SfolderPath{iD}= fullfile(MainFldr,(fldrnames(iD).name),(infoD1.subject));
        cd(folderpath);
        if ~exist(SfolderPath{iD}, 'dir')
            % If the folder doesn't exist, create it
            mkdir(SfolderPath{iD});
        end
    end
    %Assigning positionsto the figures
    fig_position1=[200 200 800 600];
    fig_position2=[1000 200 800 600];
    fig_position3=[1800 200 800 600];
    stamp=sprintf([(infoD1.subject),' ',(lfp_chan),'D1', newline,(lfp_chanD2),'D2']);

    %%% Plot Bp (with Theta SG and FG)
    fig1=figure(1);
    fig1.Name='Raw Band Power';fig1.Position= fig_position1;
    % plot Theta
    subplot(1, 3, 1)
    y_t_bp_inhib={out.(infoD1.subject).D1.t_bp_inhib',out.(infoD2.subject).D2.t_bp_inhib'};
    plot_t_bp_inhib=rm_raincloud(y_t_bp_inhib,[BC_color_genertor('Oxford_blue'); BC_color_genertor('Archt_green')]);
    title(['Band Power' newline 'Theta']);
    legend({'D1','', 'D2'})
    % plot SG
    subplot(1, 3, 2)
    y_sg_bp_inhib={out.(infoD1.subject).D1.sg_bp_inhib',out.(infoD2.subject).D2.sg_bp_inhib'};
    plot_sg_bp_inhib=rm_raincloud(y_sg_bp_inhib,[BC_color_genertor('Oxford_blue'); BC_color_genertor('Archt_green')]);
    title(['Band Power' newline 'SG']);
    legend({'D1','', 'D2'})
    % plot FG
    subplot(1, 3, 3)
    y_fg_bp_inhib={out.(infoD1.subject).D1.fg_bp_inhib',out.(infoD2.subject).D2.fg_bp_inhib'};
    plot_fg_bp_inhib=rm_raincloud(y_fg_bp_inhib,[BC_color_genertor('Oxford_blue'); BC_color_genertor('Archt_green')]);
    title(['Band Power' newline 'FG']);
    legend({'D1','', 'D2'})
    annotation('textbox',[.85 .4 .1 .15], 'String', stamp);
    %Save it
    fileName=[(infoD1.subject),'_',lfp_chan(1:4),'.png'];
    fig1path=fullfile(SfolderPath{3},fileName);
    saveas(fig1,fig1path)

    %%% Plot Bp Normalized(with Theta SG and FG)
    fig2=figure(2);
    fig2.Name='Normalized Band Power';fig2.Position= fig_position2;
    % plot Theta
    subplot(1, 3, 1)
    y_t_bp_inhib_norm={out.(infoD1.subject).D1.t_bp_inhib_norm',out.(infoD2.subject).D2.t_bp_inhib_norm'};
    plot_t_bp_inhib_norm=rm_raincloud(y_t_bp_inhib_norm,[BC_color_genertor('Oxford_blue'); BC_color_genertor('Archt_green')]);
    title(['Normalized Band Power' newline 'Theta']);
    legend({'D1','', 'D2'})
    % plot SG
    subplot(1, 3, 2)
    y_sg_bp_inhib_norm={out.(infoD1.subject).D1.sg_bp_inhib_norm',out.(infoD2.subject).D2.sg_bp_inhib_norm'};
    plot_sg_bp_inhib_norm=rm_raincloud(y_sg_bp_inhib_norm,[BC_color_genertor('Oxford_blue'); BC_color_genertor('Archt_green')]);
    title(['Normalized Band Power' newline 'SG']);
    legend({'D1','', 'D2'})
    % plot FG
    subplot(1, 3, 3)
    y_fg_bp_inhib_norm={out.(infoD1.subject).D1.fg_bp_inhib_norm',out.(infoD2.subject).D2.fg_bp_inhib_norm'};
    plot_fg_bp_inhib_norm=rm_raincloud(y_fg_bp_inhib_norm,[BC_color_genertor('Oxford_blue'); BC_color_genertor('Archt_green')]);
    title(['Normalized Band Power' newline 'FG']);
    legend({'D1','', 'D2'})
    annotation('textbox',[.85 .4 .1 .15], 'String', stamp);
    %Save it
    fileName=[(infoD1.subject),'_',lfp_chan(1:4),'.png'];
    fig2path=fullfile(SfolderPath{2},fileName);
    saveas(fig2,fig2path)

    %%% Plot ModIdx (with SG and FG)
    fig3=figure(3);
    fig3.Name='Modulation Index';fig3.Position= fig_position3;
    % plot SG
    subplot(1, 2, 1)
    y_modidx_SG_inhib={out.(infoD1.subject).D1.modidx_SG_inhib',out.(infoD2.subject).D2.modidx_SG_inhib'};
    plot_modidx_SG_inhib=rm_raincloud(y_modidx_SG_inhib,[BC_color_genertor('Oxford_blue'); BC_color_genertor('Archt_green')]);
    title(['Modulation Index' newline 'SG']);
    legend({'D1','', 'D2'})
    % plot FG
    subplot(1, 2, 2)
    y_modidx_FG_inhib={out.(infoD1.subject).D1.modidx_FG_inhib',out.(infoD2.subject).D2.modidx_FG_inhib'};
    plot_modidx_FG_inhib=rm_raincloud(y_modidx_FG_inhib,[BC_color_genertor('Oxford_blue'); BC_color_genertor('Archt_green')]);
    title(['Modulation Index' newline 'FG']);
    legend({'D1','', 'D2'})
    annotation('textbox',[.85 .4 .1 .15], 'String', stamp);
     %Save it
    fileName=[(infoD1.subject),'_',lfp_chan(1:4),'.png'];
    fig3path=fullfile(SfolderPath{1},fileName);
    saveas(fig3,fig3path)
    % for iF=1:3
    % figure(iF);clf;
    % end
end

    %% Obtaining the Mod IDX for the sleep periods
    % Restrict the filtered csc to the IV of SWS and REM
    theta_awake= restrict(theta_cscD1,iv_awake);
    theta_sws= restrict(theta_cscD1,iv_sws);
    theta_rem= restrict(theta_cscD1,iv_rem);
    
    SG_awake= restrict(SG_cscD1,iv_awake);
    SG_sws= restrict(SG_cscD1,iv_sws);
    SG_rem= restrict(SG_cscD1,iv_rem);

    FG_awake= restrict(FG_cscD1,iv_awake);
    FG_sws= restrict(FG_cscD1,iv_sws);
    FG_rem= restrict(FG_cscD1,iv_rem);

    % Restricting the non filtered CSC to the sleep periods
    CSC_Awk= restrict(csc_outD1, iv_awake);
    CSC_Sws= restrict(csc_outD1, iv_sws);
    CSC_Rem= restrict(csc_outD1, iv_rem);

    % Calculating the modIDX for the three sleep phases
    SGAwakeBins=BC_phase_amp_norm_bins(theta_awake,SG_awake);
    SGawakemodidx=MS_ModIdx(SGAwakeBins);
    FGAwakeBins=BC_phase_amp_norm_bins(theta_awake,FG_awake);
    FGawakemodidx=MS_ModIdx(FGAwakeBins);

    SGswsBins=BC_phase_amp_norm_bins(theta_sws,SG_sws);
    SGswsmodidx=MS_ModIdx(SGswsBins);
    FGswsBins=BC_phase_amp_norm_bins(theta_sws,FG_sws);
    FGswsmodidx=MS_ModIdx(FGswsBins);

    SGremBins=BC_phase_amp_norm_bins(theta_rem,SG_rem);
    SGremmodidx=MS_ModIdx(SGremBins);
    FGremBins=BC_phase_amp_norm_bins(theta_rem,FG_rem);
    FGremmodidx=MS_ModIdx(FGAwakeBins);

    % Storing the Mod idx in a matrix
    this_mouse_modidx=[SGawakemodidx FGawakemodidx; SGswsmodidx FGswsmodidx; SGremmodidx FGremmodidx];

    % if plot_flag
    %     figure(2001)
    %     subplot(2,1,1) %plot the modulation index for FG 
    %     plot([1 2 3],this_mouse_modidx(:,1));xticks([1:1:3]);xticklabels([{'Awake' 'SWS' 'REM'}]); ylabel('SG Mod Idx');
    %     subplot(2,1,2) %plot the modulation index for FG
    %     plot([1 2 3],this_mouse_modidx(:,2));xticks([1:1:3]);xticklabels([{'Awake' 'SWS' 'REM'}]); ylabel('FG Mod Idx');
    % end
    %sleep_mod_idx(:,:,iS)=this_mouse_modidx;
%% Computing the CoMo analysis fo  this mouse
    % cfg_como.A_step = 2; %I am using 2 for time-processing reasons
    % cfg_como.P_step = .5; %0.5
    % cfg_como.phi_bins = 18; %18
    % 
    % This_CoMo=[];
    % [This_CoMo.CoMoAwk, This_CoMo.phi_f, This_CoMo.amp_f] = MS_phase_freq(cfg_como, CSC_Awk, [4 12], [30 100]);
    % [This_CoMo.CoMoSws, phi_f, amp_f] = MS_phase_freq(cfg_como, CSC_Sws, [4 12], [30 100]);
    % [This_CoMo.CoMoRem, ~, ~] = MS_phase_freq(cfg_como, CSC_Rem, [4 12], [30 100]);
    % 
    % 
    % if plot_flag
    %     figure(203)
    %     clf
    %     subplot(1,3,1); cla;
    %     %surf(phi_f,amp_f,CoMoAwk')
    %     %zlim([-0.5e-3 2e-3])
    %     %hold on
    %     imagesc(phi_f, amp_f, CoMoAwk');
    %     set(gca, 'ydir', 'normal') %axis('xy')
    %     clim([0 10^-3.4])
    %     title('Awake')
    %     xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)');
    %     colorbar('Location', 'southoutside')
    % 
    %     subplot(1,3,2); cla;
    %     %surf(phi_f,amp_f,CoMoAwk')
    %     %zlim([-0.5e-3 2e-3])
    %     %hold on
    %     imagesc(phi_f, amp_f, CoMoSws');
    %     set(gca, 'ydir', 'normal') %axis('xy')
    %     %clim([0 10^-3.4])
    %     title('SWS')
    %     xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)');
    %     colorbar('Location', 'southoutside')
    % 
    %     subplot(1,3,3); cla;
    %     imagesc(phi_f, amp_f, CoMoRem');
    %     set(gca, 'ydir', 'normal') %axis('xy')
    %     %clim([0 10^-3.4])
    %     title('REM')
    %     xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)');
    %     colorbar('Location', 'southoutside')
    %     SetFigure([], gcf)
    %     maximize
        
    %end
%% Analizying the powerband of the sleep periods
bpower_slp=[];
bpower_slp.awk.Tta = bandpower(CSC_Awk.data,CSC_Awk.cfg.hdr{1}.SamplingFrequency, [5 12]);
bpower_slp.awk.SG = bandpower(CSC_Awk.data,CSC_Awk.cfg.hdr{1}.SamplingFrequency, [30 58]);
bpower_slp.awk.FG = bandpower(CSC_Awk.data,CSC_Awk.cfg.hdr{1}.SamplingFrequency, [60 100]);

bpower_slp.sws.Tta = bandpower(CSC_Sws.data,CSC_Sws.cfg.hdr{1}.SamplingFrequency, [5 12]);
bpower_slp.sws.SG = bandpower(CSC_Sws.data,CSC_Sws.cfg.hdr{1}.SamplingFrequency, [30 58]);
bpower_slp.sws.FG = bandpower(CSC_Sws.data,CSC_Sws.cfg.hdr{1}.SamplingFrequency, [60 100]);

bpower_slp.rem.Tta = bandpower(CSC_Rem.data,CSC_Rem.cfg.hdr{1}.SamplingFrequency, [5 12]);
bpower_slp.rem.SG = bandpower(CSC_Rem.data,CSC_Rem.cfg.hdr{1}.SamplingFrequency, [30 58]);
bpower_slp.rem.FG = bandpower(CSC_Rem.data,CSC_Rem.cfg.hdr{1}.SamplingFrequency, [60 100]);
sleep_colors= [...
            0.3467    0.5360    0.6907;
            0.9153    0.2816    0.2878;
            0.4416    0.7490    0.4322];

% Creating  a structure for the bar graph to compare powers
bp_bar=[];
lclBP=[];
slpStates=fieldnames(bpower_slp);
bpCat=fieldnames(bpower_slp.(slpStates{1}));
for iS=1:length(slpStates)
    for iC=1:length(bpCat)
        lclBP=[lclBP; bpower_slp.(slpStates{iS}).(bpCat{iC})];
    end
    bp_bar=[bp_bar lclBP];
    lclBP=[];
end
    %% Ploting the mod idx over the sleep phases for comparison
    if plot_flag
        figure(3)
        clf
        x=[1 2 3];
        step=0.15;
        width=.2;
        labels=[{"Awake" "SWS" "REM"}];
        step=[step*-1 step];
        hold on
        for ii=2:-1:1
            bar((x+step(ii)),this_mouse_modidx(:,ii),'BarWidth', .2,'EdgeColor','none' );
        end
        xticks([1 2 3]);xticklabels(labels);legend({'Fast gamma' 'Slow Gamma'})
        colororder(sleep_colors);
        ylabel("Modulation index (AU)"); 
        set(gca,'fontsize', 20);
        cYLim=gca().YLim;cYLim=cYLim(2);
        new_axis=[0:(cYLim/4):cYLim];
        yticks(new_axis);
        tText=[sprintf('ModIdx of mouse %s on %s with %s' , info.subject, info.session, lfp_chan )];
        title(tText); colororder(sleep_colors);
  
    end
    if save_figures
        cd('/Users/bryancontrerasmercado/Documents/01_Projects/02_OLM_silencing/11_sleep_guide/ModIdx')

        modIdxFileName=[sprintf('%s_SleepModIDX_%s_%s.',info.subject, info.session,(lfp_chan(1:4)))];
        saveas(figure(3),[modIdxFileName 'png']);
        saveas(figure(3),[modIdxFileName 'fig']);
    end
    % FG_pow = BC_power(FG_csc);
    % FG_phi = angle(hilbert(FG_csc.data));
%% Plotting the band powers
if plot_flag
    figure(2)
    bar(bp_bar,'EdgeColor','none');
    legend({"Awake" "SWS" "REM"});
    xticklabels({'Theta' 'SG' 'FG'});
    ylabel('Band power (V)');
    tText=[sprintf('Band power of mouse %s on %s with %s' , info.subject, info.session, lfp_chan )];
    title(tText); colororder(sleep_colors);
end
if save_figures
    cd('/Users/bryancontrerasmercado/Documents/01_Projects/02_OLM_silencing/11_sleep_guide/BandPower')
    BPFileName=[sprintf('%s_SleepModIDX_%s_%s.',info.subject, info.session,(lfp_chan(1:4)))];
    if save_figures
        cd('/Users/bryancontrerasmercado/Documents/01_Projects/02_OLM_silencing/11_sleep_guide/BandPower');
        BpFileName=[sprintf('%s_SleepBandBower_%s_%s.',info.subject, info.session,(lfp_chan(1:4)))];
        saveas(figure(2),[BpFileName 'png']);
        saveas(figure(2),[BpFileName 'fig']);
    end
end
    %% Getting the percentage of sleep sates
    [y,x]=histcounts(hypnoD1.data,[0.5:1:3.5]); %Y=count, x=rnage values
    y_per=(y/sum(y))*100; %Percentage of Wake, SWS and REM
    sleep_time_sec=(y./(csc_s.cfg.hdr{1,1}.SamplingFrequency))';

    if plot_flag
        figure(222)
        clf
        subplot(1,2,1)
        b=bar([1:1:3],y_per);
        b.FaceColor = 'flat';
        cord = linspecer(5);
        %cord=orderedcolors('reef');
        %Color assignment
        for ci=1:3
            b.CData(ci,:) = cord(ci,:);
        end
        b.LineStyle= "none";
        xticklabels({'Wake','SWS','REM'})
        ylim([0 80]);
        % Add percentage labels on top of each bar
        for ii = 1:numel(y_per)
            text(ii, y_per(ii), sprintf('%.1f%%', y_per(ii)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
        ylabel('Time spent per state [%]')

        %Donut
        subplot(1,2,2)
        d=donutchart(y_per, {'Wake','SWS','REM'})
        sleep_colors= [...
            0.3467    0.5360    0.6907;
            0.9153    0.2816    0.2878;
            0.4416    0.7490    0.4322];
        colororder(sleep_colors);
        sgtitle(sprintf('Sleep phases %s on  session %s', info.subject, info.session), 'fontweight', 'bold', 'fontsize', 16);
    end
    %-----TO DO----- add a label of the total time that the mice spent sleeping

    %% Extracting the  sleep data
    out.(info.subject).(info.session).sleep.percetages=y_per;
    out.(info.subject).(info.session).sleep.times_sec=sleep_time_sec;
    out.(info.subject).(info.session).sleep.mod_idx=this_mouse_modidx;
    %out.(info.subject).(info.session).sleep.CoMo=This_CoMo;
    out.(info.subject).sleep_labels=hypnoD1.labels;
    out.(info.subject).(info.session).sleep.powers=bpower_slp;


    %% END OF SLEEP ANALYSIS
    %% Tracking movement and setting up intervals when it explored the objects
    if info.session == "D2"
        % Load the DLC data
        distances=[];

        pos = MS_DLC2TSD_divide(cd, [], [9.3 9.3]); %Need to ask Eric how to solve for the size of the cage. For now we will assume its 4.5
        files=fieldnames(pos);
        files_idx=(find((contains((fieldnames(pos)),'File'))));
        nfiles=length(files(files_idx));
        files=files(files_idx);

        %Correct time and remove time where the mouse was not in OF
        for iF=1:nfiles;
            pos.(files{iF}).tvec=pos.(files{iF}).tvec-pos.(files{iF}).tvec(1);
            %Correct positions of objects to the mean of the location for both objects
            for ii= 5:8
                pos.(files{iF}).data(ii,:)=mean(pos.(files{iF}).data(ii,:));
            end
            %Calculate the distance between the objects and the mouse nose
            %A. Calculate the distance between the nose and object A (Object A is the
            %one on the left)
            for iframe= length(pos.(files{iF}).data(1,:)):-1:1
                x1=pos.(files{iF}).data(1,iframe);
                y1=pos.(files{iF}).data(2,iframe);
                x2=pos.(files{iF}).data(5,iframe);
                y2=pos.(files{iF}).data(6,iframe);
                d=sqrt(((x2-x1)^2)+((y2-y1)^2));
                distances.(files{iF})(1,iframe)=d;
                if iframe==1
                    clear iframe;clear x1;clear x2;clear y1;clear y2;clear d;
                end
            end
            %B.Same for object B (Object B is the one on the right)
            for iframe= length(pos.(files{iF}).data(1,:)):-1:1
                x1=pos.(files{iF}).data(1,iframe);
                y1=pos.(files{iF}).data(2,iframe);
                x2=pos.(files{iF}).data(7,iframe);
                y2=pos.(files{iF}).data(8,iframe);
                d=sqrt(((x2-x1)^2)+((y2-y1)^2));
                distances.(files{iF})(2,iframe)=d;
                if iframe==1
                    clear iframe;clear x1;clear x2;clear y1;clear y2;clear d;
                end
            end       
        end
        distances.labels=(pos.File1.label(3:4))';


        %%% You are here in this function
%% Creating a heatmap of the position of the mouse
% Number of bins for the heatmap (resolution of the grid)
if plot_flag
    numBins = 40;
    sec2delete=40;
    GaussianS=1.2;
    for iF=1:nfiles
        X=pos.(files{iF}).data(1,sec2delete*30:end);
        Y=pos.(files{iF}).data(2,sec2delete*30:end);
        [heatmapData,C] = hist3([X', Y'], [numBins, numBins]);
        heatmapData=(heatmapData./length(X))*100;
        heatmapData_smoothed = imgaussfilt(heatmapData, GaussianS);
        % Plot the heatmap
        figure (1600+iF);
        imagesc(C{1}', C{2}',heatmapData_smoothed'); % Transpose because hist3's output is transposed
        hold on
        plot(X,Y,'w')
        clim([0 0.2])
        %set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
        colorbar; % Show color scale
        title('Heatmap of Mouse Position in Open Field');
        %xlabel('Position (cm) ');
        %ylabel('Position (cm)');
    end
end
        %% Plot the position of the mouse
        %%---To do--- Adapt this cell to the new structures
       
        radious=5;
        minFrames=5; % The minimum number of frames where the mouse is in radious
        
        for iF=1:nfiles
            %Initialize figure
            if plot_flag
            figure (20+iF)
            end
            %Actual position
            %subplot(8,1,1:4);
            %plot the actual position
            time_frames=[1:length(pos.("File"+iF).data)];
            [aX,aY]=BC_Circle_plot(plot_flag,radious,pos.("File"+iF).data(5,1),pos.("File"+iF).data(6,1), BC_color_genertor('Swamp_green'));
            hold on
            [bX,bY]=BC_Circle_plot(plot_flag,radious,pos.("File"+iF).data(7,1),pos.("File"+iF).data(8,1), BC_color_genertor('Torment_blue'));
            inA= inpolygon(pos.("File"+iF).data(1,:),pos.("File"+iF).data(2,:),aX,aY);
            aidx=find(inA);
            aInt=BC_object_interaction_intervals(aidx,minFrames);%Get periods where the mouse spends more than 5 frames (0.33 sec)
            timeA=(sum(diff(aInt,1,2))/30);
            inB= inpolygon(pos.("File"+iF).data(1,:),pos.("File"+iF).data(2,:),bX,bY);
            bidx=find(inB);
            bInt=BC_object_interaction_intervals(bidx,minFrames);
            timeB=(sum(diff(bInt,1,2))/30);
            fprintf('The mouse <strong>%s</strong> interacted <strong>%d</strong> times with the object 1 for a total of <strong>%.2f</strong> seconds and <strong>%d</strong> times with the object 2 for a total of <strong>%.2f</strong> seconds in the <strong>%d</strong> experiment \n ', info.subject ,size(aInt,1),timeA,size(bInt,1), timeB, iF)

            %Lets verify that these idx correspond to times the mouse was with the objects

            if plot_flag==1
                if length(aInt) > length(bInt)
                    n = length(aInt); % Number of colors

                else
                    n = length(bInt); % Number of colors

                end
                cmap = autumn(n); % You can replace 'parula' with any other colormap name 'parula', 'jet', 'hsv', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'lines', and 'colorcube'.

                for jj=1:1:size(aInt,1)
                    scatter(pos.("File"+iF).data(1,aInt(jj,1):aInt(jj,2)), pos.("File"+iF).data(2,aInt(jj,1):aInt(jj,2)), 'o', 'MarkerFaceColor', cmap(jj,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
                end
                for jj=1:1:size(bInt,1)
                    scatter(pos.("File"+iF).data(1,bInt(jj,1):bInt(jj,2)), pos.("File"+iF).data(2,bInt(jj,1):bInt(jj,2)), 'o', 'MarkerFaceColor', cmap(jj,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
                end

                %Position of the two objects
                scatter(pos.("File"+iF).data(5,1),pos.("File"+iF).data(6,1) , 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'none')
                scatter(pos.("File"+iF).data(7,1),pos.("File"+iF).data(8,1) , 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
                xlim([min(pos.("File"+iF).data(1,:)) max(pos.("File"+iF).data(1,:))]);
                ylim([min(pos.("File"+iF).data(2,:)) max(pos.("File"+iF).data(2,:))]);

                if video_flag==0;
                    scatter(pos.("File"+iF).data(1,:), pos.("File"+iF).data(2,:), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.01);%Position of the mouse during the whole experiment in blue
                else
                    %%Trying the video
                    for ii=1400:-1:900
                        scatter(pos.("File"+iF).data(1,ii), pos.("File"+iF).data(2,ii), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.02);
                        F(ii) = getframe(gcf);
                        drawnow
                    end
                    writerObj = VideoWriter('myVideo.avi');
                    writerObj.FrameRate = 10;
                    open(writerObj);
                    % write the frames to the video
                    for i=1:length(F)
                        % convert the image to a frame
                        frame = F(i) ;
                        writeVideo(writerObj, frame);
                    end
                    % close the writer object
                    close(writerObj);
                end
            end
            % %%% Rest of the plot
            % %plot X position of the mouse
            % subplot(8,1,5);
            % plot(time_frames,pos.("File"+iF).data(1,:));
            % xlim([0,time_frames(end)])
            % %xlabel('Time(s)');
            % ylabel('X pos (cm)');
            % %plot Y position of the mouse
            % subplot(8,1,6);
            % plot(time_frames,pos.("File"+iF).data(2,:), 'Color',BC_color_genertor('Oxford_blue'));
            % xlim([0,time_frames(end)])
            % %xlabel('Time(s)');
            % ylabel('Y pos (cm)');
            % %Plot the distances to object 1
            % subplot(8,1,7);
            % plot(time_frames, distances.("File"+iF)(1,:), 'm');
            % yline(5, 'Color', 'r')
            % xlim([0,time_frames(end)])
            % %xlabel('Time(s)');
            % ylabel('Distance to object A');
            % %Plot the distnaces to object 2
            % subplot(8,1,8);
            % plot(time_frames, distances.("File"+iF)(2,:),'r');
            % yline(5, 'Color', 'r')
            % xlim([0,time_frames(end)])
            % xlabel('Time (s)');
            % ylabel('Distance to object B');
            % 
            % hold off

            %Extracting the data
            out.(info.subject).(info.session).("Behavior"+iF).nObjAInteractions=length(aInt);
            out.(info.subject).(info.session).("Behavior"+iF).nObjBInteractions=length(bInt);
            out.(info.subject).(info.session).("Behavior"+iF).IntervalsObjAInteractions=aInt;
            out.(info.subject).(info.session).("Behavior"+iF).IntervalsObjBInteractions=bInt;
            out.(info.subject).(info.session).("Behavior"+iF).ObjATime= timeA;
            out.(info.subject).(info.session).("Behavior"+iF).ObjBTime= timeB;
            out.(info.subject).(info.session).("Behavior"+iF).pos=pos.("File"+iF);

            clear wake_t;


        end
    end


% Formating for saving output
if save_output
    cd(inter_dir)
    save(["out-" + date + ".mat"],'out')   
end



%% Lets assign the intervals where te mice is inisde the radious
%
% %disp(cmap);
% %Calculate the angle of the mouse head and if it is pointing towards the object
%
% %Figure out intervals where the mouse nose was close to the object and the
% %head direction is pointing towards the closest object
% %Figure out whih data points are inside the radious
%
% for iF=1:nfiles
%     inA= inpolygon(pos.("File"+iF).data(1,:),pos.("File"+iF).data(2,:),aX,aY);
%     aidx=find(inA);
%     aInt=BC_object_interaction_intervals(aidx,minFrames);%Get periods where the mouse spends more than 5 frames (0.33 sec)
%     timeA=(sum(diff(aInt,1,2))/30);
%     inB= inpolygon(pos.("File"+iF).data(1,:),pos.("File"+iF).data(2,:),bX,bY);
%     bidx=find(inB);
%     bInt=BC_object_interaction_intervals(bidx,minFrames);
%     timeB=(sum(diff(bInt,1,2))/30);
%     fprintf('The mouse <strong>%s</strong> interacted <strong>%d</strong> times with the object 1 for a total of <strong>%.2f</strong> seconds and <strong>%d</strong> times with the object 2 for a total of <strong>%.2f</strong> seconds in the <strong>%d</strong> experiment \n ', info.subject ,size(aInt,1),timeA,size(bInt,1), timeB, iF)
%
%     %Lets verify that these idx correspond to times the mouse was with the objects
%     for jj=1:1:size(aInt,1)
%         scatter(pos.("File"+iF).data(1,aInt(jj,1):aInt(jj,2)), pos.("File"+iF).data(2,aInt(jj,1):aInt(jj,2)), 'o', 'MarkerFaceColor', cmap(jj,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
%     end
%     for jj=1:1:size(bInt,1)
%         scatter(pos.("File"+iF).data(1,bInt(jj,1):bInt(jj,2)), pos.("File"+iF).data(2,bInt(jj,1):bInt(jj,2)), 'o', 'MarkerFaceColor', cmap(jj,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
%     end
% end
%% Automate collection of data for more individuals
% Is 5 mc normally the case? how can I define the perimeter of my object ?
% Would this have to be done based on subject?

%Remove periods where the mouse is on top of the object (if any)