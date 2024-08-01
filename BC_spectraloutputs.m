function out=BC_EpochsSpectralOutputsCondition(inputs, )
%% The function BC_EpochsSpectralOutputsCondition:
%            [Provide a general description of the function here]
%   Inputs:
%            -input1[type]: description 
% 
% 
%   Outputs:
%             -output1[type]: description 
% 
% 
%  First version BC 30-Jul-2024 
%% Analysis part
    
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

    for ii = length(iv_inhb.tstart):-1:1
        if (iv_inhb.tend(ii) - iv_inhb.tstart(ii)) < min_trial_dur
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
        else
            this_csc = restrict(csc, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            
            % get the power
            t_bp =  bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [5 12]);
            sg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [30 58]);
            fg_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [60 100]);
            
            ref_bp = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 50]);
            
            t_bp_inhib(ii)= t_bp;
            sg_bp_inhib(ii) = sg_bp;
            fg_bp_inhib(ii) = fg_bp;
            
            t_bp_inhib_norm(ii)= t_bp/ ref_bp;
            sg_bp_inhib_norm(ii) = sg_bp/ ref_bp;
            fg_bp_inhib_norm(ii) = fg_bp/ ref_bp;
            
            %Create a table with epoch#, condition, amplitude, velocity
            %The following code create the veactor for the creation of the table
            this_theta_amp=restrict(theta_amp_phi, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            this_Sg_amp=restrict(SG_amp_phi, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            this_Fg_amp=restrict(FG_amp_phi, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            this_pos= restrict(pos, iv_inhb.tstart(ii), iv_inhb.tend(ii));
            
            %Creating a copy of the amplitude vectors so I can use all the values ion a second analysis
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
            %Extracting the data from the structure and combining vector for all the epcohs
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
    end
    %putting the speed-pwr into a table so its easier to extract later
    spdPwrTbl_inhb=table( spd_ihb',tta_amp_ihb',tta_pwr_ihb', tta_phi_ihb',sg_amp_ihb',sg_pwr_ihb',sg_phi_ihb',fg_amp_ihb',fg_pwr_ihb',fg_phi_ihb',epoch_ihb', 'VariableNames',{'spd','tta_amp','tta_pwr', 'tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi','epoch'});
    %creating a table for all the variables of the amplitude        
    spdPwrTbl_inhbTotal= table( tta_amp_ihbTotal',tta_pwr_ihbTotal', tta_phi_ihbTotal',sg_amp_ihbTotal',sg_pwr_ihbTotal',sg_phi_ihbTotal',fg_amp_ihbTotal',fg_pwr_ihbTotal',fg_phi_ihbTotal',epoch_ihbTotal', 'VariableNames',{'tta_amp','tta_pwr', 'tta_phi','sg_amp','sg_pwr','sg_phi','fg_amp','fg_pwr','fg_phi','epoch'});      
          

%% outputs
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