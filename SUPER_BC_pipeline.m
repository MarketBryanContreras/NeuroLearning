%% Super_BC_ phase_mod

%% initialize

% CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2';
% mvdm_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared';
% BC_dir = 'C:\Users\ecarm\Documents\GitHub\NeuroLearning';
% 
% addpath(genpath(mvdm_dir));
% addpath(genpath(CEH2_dir));
% addpath(genpath(BC_dir));

% data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw';
% inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter';

data_dir = 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw';
inter_dir = 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter';

cd(data_dir)

%% parameters
plot_flag = 1; % switch to 0 if you want to supress verification figures.
time_maze_start = 30;
min_trial_dur = 2;

%% loop over sessions

cd(data_dir)

% get all sessions with 'INHIBITION'
inhib_dir = dir('*INHIBITION*');


% load data from raw

for iS = 13%1:length(inhib_dir)
    
    cd([inhib_dir(iS).folder filesep inhib_dir(iS).name])
    
    parts = strsplit(inhib_dir(iS).name, '_');
    
    info.subject = parts{1};
    info.sess = parts{5};
    info.date = [ parts{2} '_' parts{3} '_' parts{4}];
    
    cfg_csc = [];
    cfg_csc.desired_sampling_frequency = 2000;
    
    if strcmpi(info.subject, 'BC1602')
        cfg_csc.fc ={'CSC4.ncs'};
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    elseif strcmpi(info.subject, 'BC1807')
        cfg_csc.fc ={'CSC4.ncs'};
        pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).';
    end
    
    [csc, evts, pos] = BC_load_NLX(cfg_csc);
    
    fs = csc.cfg.hdr{1}.SamplingFrequency;
    sSampCsc= time_maze_start*fs; %This is the idx of the samplimg at sec 30
    restrictIvCsc=iv(csc.tvec(sSampCsc),csc.tvec(end)); %Creates and IV from the time the mouse was placed in the maze to the end of the recording
    csc= restrict(csc,restrictIvCsc);
    % Correction of time
    strt=csc.tvec(1);
    csc.tvec= csc.tvec-strt; %correct time csc
    nn=size(evts.t,2);
    evts.t{1,nn-1}=evts.t{1,nn-1}-strt; % Coreccting laser events start
    evts.t{1,nn}=evts.t{1,nn}-strt; % Coreccting laser events end
    
    %Restiction of time when the mouse is plces in the maze
    sSamp=time_maze_start*30;
    restrictIvPos=iv(pos.tvec(sSamp),pos.tvec(end));
    pos=restrict(pos,restrictIvPos);
    pos.tvec= pos.tvec-pos.tvec(1);
    
    
    iv_inhb = MS_get_evts_off(evts, pattern);
    
    
    %% trial split
    
    [iv_inhb,iv_noInhb, iv_running] = BC_LT_trialfun(pos, iv_inhb, 1);
    
    %% Filter
    % filter the LFP in the theta band
    cfg_filt_t = [];
    cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [5 12]; % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 3; %type filter order
    cfg_filt_t.display_filter = 0; % use this to see the fvtool
    theta_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
    theta_pow = BC_power(theta_csc);
    theta_phi = angle(hilbert(theta_csc.data));
    
    % Filter the LFP in the slow gamma band
    cfg_filt_t = [];
    cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [30 58]; % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 4; %type filter order
    cfg_filt_t.display_filter = 0; % use this to see the fvtool
    SG_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
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
    
    %% get phase mod over trial
    
    SGphiAmpNormInhb=BC_phase_amp_norm_bins(theta_inhb,SG_inhb);
    SGInhb_modidx=MS_ModIdx(SGphiAmpNormInhb);
    BC_plot_modidx(SGphiAmpNormInhb,BC_color_genertor('Archt_green'),SGInhb_modidx,'SG', 'silencing' )
    
    SGphiAmpNormNoInhb = BC_phase_amp_norm_bins(theta_noinhb,SG_noinhb);
    SGNoInhb_modidx=MS_ModIdx(SGphiAmpNormNoInhb);
    BC_plot_modidx(SGphiAmpNormNoInhb,BC_color_genertor('Powder_blue'),SGNoInhb_modidx, 'SG', 'no silencing')
    
    [SGshift_mean, SGshift_std]=LTshifted_meanModIdx(theta_running,SG_running);
    
    
    FGphiAmpNormInhb=BC_phase_amp_norm_bins(theta_inhb,FG_inhb);
    FGInhb_modidx=MS_ModIdx(FGphiAmpNormInhb);
    BC_plot_modidx(FGphiAmpNormInhb,BC_color_genertor('Archt_green'),FGInhb_modidx,'FG', 'silencing' )
    
    %NoInhb
    FGphiAmpNormNoInhb = BC_phase_amp_norm_bins(theta_noinhb,FG_noinhb);
    FGNoInhb_modidx=MS_ModIdx(FGphiAmpNormNoInhb);
    BC_plot_modidx(FGphiAmpNormNoInhb,BC_color_genertor('Powder_blue'),FGNoInhb_modidx, 'FG', 'no silencing')
    
    %Shifted
    [FGshift_mean,FGshift_std]=LTshifted_meanModIdx(theta_running,FG_running);
    
    
    %% epochs by epochs analysis
    
    % inhibition
    for ii = length(iv_inhb.tstart):-1:1
        if (iv_inhb.tend(ii) - iv_inhb.tstart(ii)) < min_trial_dur
            t_bp_inhib(ii)= NaN;
            sg_bp_inhib(ii) = NaN;
            fg_bp_inhib(ii) = NaN;
            this_SG_inhib(ii)= NaN;
            this_FG_inhib(ii)= NaN;
        else
            this_csc = restrict(csc, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            
            % get the power
            t_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [5 12]);
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
            this_SG_inhib(ii) =MS_ModIdx(this_SG_phi_amp);
            %FG
            this_FG_phi_amp=BC_phase_amp_norm_bins(this_th,this_fg);
            this_FG_inhib(ii) =MS_ModIdx(this_FG_phi_amp);
                        
        end
        
    end
    
    % NO inhibition
    for ii = length(iv_noInhb.tstart):-1:1
        if (iv_noInhb.tend(ii) - iv_noInhb.tstart(ii)) < min_trial_dur
            t_bp_noinhib(ii)= NaN;
            sg_bp_noinhib(ii) = NaN;
            fg_bp_noinhib(ii) = NaN;
            this_SG_noinhib(ii)= NaN;
            this_FG_noinhib(ii)= NaN;
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
            this_SG_noinhib(ii) =MS_ModIdx(this_SG_phi_amp);
            %FG
            this_FG_phi_amp=BC_phase_amp_norm_bins(this_th,this_fg);
            this_FG_noinhib(ii) =MS_ModIdx(this_FG_phi_amp);
        end
        
    end
 
    % Running
    for ii = length(iv_running.tstart):-1:1
        if (iv_running.tend(ii) - iv_running.tstart(ii)) < min_trial_dur
            t_bp_run(ii)= NaN;
            sg_bp_run(ii) = NaN;
            fg_bp_run(ii) = NaN;
            this_SG_run(ii)= NaN;
            this_FG_run(ii)= NaN;
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
            this_SG_run(ii) =MS_ModIdx(this_SG_phi_amp);
            %FG
            this_FG_phi_amp=BC_phase_amp_norm_bins(this_th,this_fg);
            this_FG_run(ii) =MS_ModIdx(this_FG_phi_amp);
        end
        
    end
    
    %% save outputs
    %inhb
    out.(info.subject).(info.sess).t_bp_inhib = t_bp_inhib; 
    out.(info.subject).(info.sess).sg_bp_inhib = sg_bp_inhib;
    out.(info.subject).(info.sess).fg_bp_inhib = fg_bp_inhib;
    out.(info.subject).(info.sess).this_SG_inhib = this_SG_inhib;
    out.(info.subject).(info.sess).this_FG_inhib = this_FG_inhib;
    
    %noinb
    out.(info.subject).(info.sess).t_bp_noinhib = t_bp_noinhib; 
    out.(info.subject).(info.sess).sg_bp_noinhib = sg_bp_noinhib;
    out.(info.subject).(info.sess).fg_bp_noinhib = fg_bp_noinhib;
    out.(info.subject).(info.sess).this_SG_noinhib = this_SG_noinhib;
    out.(info.subject).(info.sess).this_FG_noinhib = this_FG_noinhib;
    
    %running
    out.(info.subject).(info.sess).t_bp_noinhib = t_bp_noinhib; 
    out.(info.subject).(info.sess).sg_bp_noinhib = sg_bp_noinhib;
    out.(info.subject).(info.sess).fg_bp_noinhib = fg_bp_noinhib;
    out.(info.subject).(info.sess).this_SG_noinhib = this_SG_noinhib;
    out.(info.subject).(info.sess).this_FG_noinhib = this_FG_noinhib;
    %% zscore the values
    
    
end





