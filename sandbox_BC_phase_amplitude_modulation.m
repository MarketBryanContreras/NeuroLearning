%% sandbox  BC laser pulse start-end
% Go to some data
%BC1602
%data_dir = "C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC1602_07_12_2023_D4_INHIBITION";
%BC1807
data_dir = 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC1807_07_12_2023_D4_INHIBITION';
cd (data_dir)

%% Assigning color pallete
Archt_green = [0.530 0.820 0.645];
Oxford_blue = [0.039 0.137 0.259];
Powder_blue = [0.682 0.772 0.921];
Burt_orange= [0.758 0.348 0.249];
Red_crayola= [0.937 0.176 0.337];
Web_orange= [1.00 0.678 0.020];
Archt_green_alpha = [0.530 0.820 0.645 0.3];
Oxford_blue_alpha= [0.039, 0.137, 0.259, 0.300];
Powder_blue_alpha= [0.682 0.772 0.921 0.3];
Burt_orange_alpha= [0.933 0.423 0.302 0.3];
Web_orange_alpha= [1.00 0.678 0.020 0.3];
%% load some data
evts = LoadEvents([]); 
cfg =[];
cfg.fc = {'CSC6.ncs'};

csc = MS_LoadCSC(cfg);
corrected_time=csc.tvec-min(csc.tvec); % creates n array of time corrected for the time that neuralynx started recording
%% get some events of interest
%on_ts2 = evts.t{(strcmp('TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).',evts.label))};
on_idx = find(contains(evts.label, '(0x0002)'));  % find the correct label. 
on_ts = evts.t{on_idx}; % grab those times
off_idx = find(contains(evts.label, '(0x0000)'));  % find the corresponding off label
off_ts = evts.t{off_idx}; 

% convert to 'iv' (inverval) format for simplicity
laser_iv = iv(on_ts, off_ts); 

%% alternative. Use a little function based on the cell above. 
pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).'; 
iv_out = MS_get_evts_off(evts, pattern); 

%% plot some examples
figure(101)
clf

hold on

max_csc = max(csc.data); 
min_csc = min(csc.data); 
offset = (max_csc - min_csc);%I use these lines to calculate the y limits (+- 10% of these value) 
range_final= [min_csc-offset*0.1 max_csc+offset*0.1];
height=(max_csc+offset*0.1)-(min_csc-offset*0.1);

% c_ord = linspecer(length(on_ts));
% 
% for ii = 1:length(laser_iv.tstart)
%     
%     rectangle('position', [laser_iv.tstart(ii)-min(csc.tvec), (min_csc-offset*0.1), (laser_iv.tend(ii)-min(csc.tvec)) - (laser_iv.tstart(ii)-min(csc.tvec)), height], 'facecolor', Archt_green_alpha, 'edgecolor', Archt_green_alpha)
%     
% end
for ii = 1:length(laser_iv.tstart)
    
    h=fill([(laser_iv.tstart(ii)-min(csc.tvec)) (laser_iv.tend(ii)-min(csc.tvec)) (laser_iv.tend(ii)-min(csc.tvec)) (laser_iv.tstart(ii)-min(csc.tvec))], [range_final(1) range_final(1) range_final(2) range_final(2)],Archt_green, 'FaceAlpha', 0.3, 'LineStyle',"none");
   
end
% The function fill to create rentangles use the arguments in the following fashion:fill([xmin xmax xmax xmin], [ymin ymin ymax ymax],Burt_orange, 'FaceAlpha', 0.3, 'LineStyle',"none")
%Adventage of fill over rentanlge is that it can be used in the label of the figure 
% plot the data on top
corrected_time=csc.tvec-min(csc.tvec); % creates n array of time corrected for the time that neuralynx started recording
plot(corrected_time, csc.data, 'Color', Oxford_blue );
xlim([corrected_time(1) corrected_time(end)]);
ylim(range_final);
xlabel('Time (s)');
ylabel('Potential (µV)');
title('Periods of inhibition in linear track', 'Fontsize', 20);
lasleg=legend(h, 'Laser 520nm');
legendFontSize = 14; % Adjust the font size as needed
set(lasleg, 'FontSize', legendFontSize);
legend boxoff ;
%print('-dpng', 'Inhibition_epochs.png', '-r300'); %code to save the figure

%with 300 dpi
%% gettting the power and phase

% filter the LFP in the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [5 12]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 3; %type filter order
cfg_filt_t.display_filter = 0; % use this to see the fvtool

theta_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
theta_amp = abs(hilbert(theta_csc.data)); % get the amplitude
theta_phi  = angle(hilbert(theta_csc.data(1,:))); %get the phase
%theta_csc.data = theta_csc.data(1,:); 



%% Get the amplitud of slow gamma

% filter the LFP in the slow gamma band
cfg_filt_t = [];
cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [30 58]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 4; %type filter order
cfg_filt_t.display_filter = 0; % use this to see the fvtool

SG_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using

SG_amp = abs(hilbert(SG_csc.data)); % get the amplitud


%% Get ampllitud of fast gamma

% filter the LFP in the fast gamma band
cfg_filt_t = [];
cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [70 90]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 4; %type filter order
cfg_filt_t.display_filter = 0; % use this to see the fvtool

FG_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
FG_amp = abs(hilbert(FG_csc.data)); % get the amplitude for Fast Gamm


%% Try some Phase-amp coupling

mod_th_g = MS_ModIdx_win(theta_csc, SG_csc, 5*theta_csc.cfg.hdr{1}.SamplingFrequency, iv_out);
title('Theta - SG Mod Index')
%print('-dpng', 'Thetha_slow_gamma_modulation_index.png', '-r300'); %code to save the figure

mod_th_fg = MS_ModIdx_win(theta_csc, FG_csc, 5*theta_csc.cfg.hdr{1}.SamplingFrequency, iv_out);
title('Theta - FG Mod Index')
%print('-dpng', 'Thetha_fast_gamma_modulation_index.png', '-r300'); %code to save the figure
%% Plot only the modulation indexes
figure(101)
clf
plot(corrected_time, mod_th_g, 'Color', Powder_blue); 
hold on
plot(corrected_time, mod_th_fg, 'Color', Red_crayola); 
maxMod = max(mod_th_g); 
minMod = min(mod_th_g); 
ModOffset = (maxMod - minMod);
added_height= 0.2;
Mod_range_final= [minMod-ModOffset*added_height maxMod+ModOffset*added_height];
height=(maxMod+ModOffset*added_height)-(minMod-ModOffset*added_height);
for ii = 1:length(iv_out.tstart)
    
    h=fill([(iv_out.tstart(ii)-min(theta_csc.tvec)) (iv_out.tend(ii)-min(theta_csc.tvec)) (iv_out.tend(ii)-min(theta_csc.tvec)) (iv_out.tstart(ii)-min(theta_csc.tvec))], [Mod_range_final(1) Mod_range_final(1) Mod_range_final(2) Mod_range_final(2)],Archt_green, 'FaceAlpha', 0.3, 'LineStyle',"none");
   hold on
end
%fill([xmin xmax xmax xmin], [ymin ymin ymax ymax],Burt_orange, 'FaceAlpha', 0.3, 'LineStyle',"none")
ylabel('Mod Idx' )
xlabel('Time(s)')
legend({'Theta - SG', 'Theta - FG', 'Laser 520nm'})
xlim([corrected_time(1) corrected_time(end)])
ylim([Mod_range_final(1) Mod_range_final(2)])
graph_width= 1000;
graph_height=500;
set (gcf, 'position',[ 100 100 graph_width graph_height]);
box off;
print('-dpng', 'Thetha_fast_gamma_slow_gamma_modulation_index_only.png', '-r300'); %code to save the figure

%% Temp Mod IDX for laser on blocks. 

csc_laser = restrict(csc, iv_out);

theta_laser = restrict(theta_csc, iv_out);
SG_laser = restrict(SG_csc, iv_out);
FG_laser = restrict(FG_csc, iv_out);

mod_th_sg_l = MS_ModIdx_win(theta_laser, SG_laser, 5*theta_csc.cfg.hdr{1}.SamplingFrequency, iv_out);
title('Theta - SG Mod Index (no laser)')
%print('-dpng', 'Thetha_slow_gamma_modulation_index.png', '-r300'); %code to save the figure

mod_th_fg_l = MS_ModIdx_win(theta_laser, FG_laser, 5*theta_csc.cfg.hdr{1}.SamplingFrequency, iv_out);
title('Theta - FG Mod Index (no laser)')
%print('-dpng', 'Thetha_fast_gamma_modulation_index.png', '-r300'); %code to save the figure
% Problem with this is that is not corresoinding to the original
% modulatiuon index values, better solution would be to restrict the data
% in the modulation indeces and plot it instead of plotting feeding it into
% MS_ModIDx_win again
%% temporary comparison against the non-laser periods
iv_no_laser = InvertIV(iv_out, evts.t{1}, evts.t{2}); % grab all the non-laser intervals. 

theta_no_laser = restrict(theta_csc, iv_no_laser);
SG_no_laser = restrict(SG_csc, iv_no_laser);
FG_no_laser = restrict(FG_csc, iv_no_laser);

mod_th_sg_no_l = MS_ModIdx_win(theta_no_laser, SG_no_laser, 10*theta_csc.cfg.hdr{1}.SamplingFrequency, 0);
title('Theta - SG Mod Index (no laser)')
% c_ord = linspecer(length(on_ts));
Artch_green = [0.530 0.820 0.645];
for ii = 1:length(laser_iv.tstart)
    
    rectangle('position', [laser_iv.tstart(ii), min_csc, laser_iv.tend(ii) - laser_iv.tstart(ii), offset], 'facecolor', [Artch_green .2], 'edgecolor', [Artch_green .2])
    
end

mod_th_fg_no_l = MS_ModIdx_win(theta_no_laser, FG_no_laser, 10*theta_csc.cfg.hdr{1}.SamplingFrequency, 0);
title('Theta - FG Mod Index (no laser)')




%% get a baseline mod_value for anytime the animal is moving but not during the lasers;
% save this for when you have 

% % get the linear speed
% speed = getLinSpd([],pos); % linear speed
% 
% % Threshold speed
% cfg = []; cfg.method = 'raw'; cfg.operation = 'range'; cfg.threshold = [3 20]; % speed limit in cm/sec
% iv_move = TSDtoIV(cfg,speed); % only keep intervals with speed above thresh
% 
% 
% % restrict to movement
%    cfg_spd = [];
%     cfg_spd.method = 'raw'; 
%     cfg_spd.operation = '>'; 
%     cfg_spd.threshold = cfg.spd_thresh;
%     iv_move = TSDtoIV(cfg_spd,spd); % only keep intervals with speed above thresh
% 
%     
% 
% csc_move = restrict(csc, iv_move); 


%% slow / fast ratio


SG_FG_ratio = smooth(SG_amp, 2*theta_csc.cfg.hdr{1}.SamplingFrequency)./smooth(FG_amp, 2*theta_csc.cfg.hdr{1}.SamplingFrequency);

figure(10)
hold on
plot(csc.tvec, SG_FG_ratio)


% c_ord = linspecer(length(on_ts));
Artch_green = [0.530 0.820 0.645];
for ii = 1:length(laser_iv.tstart)
    
    rectangle('position', [laser_iv.tstart(ii), min_csc, laser_iv.tend(ii) - laser_iv.tstart(ii), 100], 'facecolor', [Artch_green .2], 'edgecolor', [Artch_green .2])
    
end
%% Creating a Phase-Amplitud plot for the times when the laser was on 
%theta_amp;
%theta_phi;
phi_bins = -pi:pi/18:pi; 
%2. create a composite of the phase_theta, amplitud_SG
[~,theta_idx]= histc(theta_phi, phi_bins); %This creates a vector with the  cooreponding bin of each phase
%3. The phases are binned according to which phase they belong to and their
%mean is calculated
phase_means= zeros(1,length(unique(theta_idx)));
for ii= 1: length(unique(theta_idx));
    phase_means(ii)= nanmean(SG_amp(theta_idx == ii));
end
%4.Each mean amplitude is then divided by the sum over the bins 
phase_means_norm= phase_means/sum(phase_means);
%5 Plot the phase amplitude coupling
figure
clf
ax(1) = subplot(3,1,1:2);
deg_bins= (phi_bins(2:end)+pi)*(180/pi); %Here I am adding pi radinas to correct for the hilber transformation and converting to bins in degrees
bar(deg_bins,phase_means_norm,'FaceColor',Oxford_blue, 'EdgeColor',Oxford_blue )
%xlabel('Theta phase(Deg)');
ylabel('SG amplitude');
title('Theta-phase-SG-amplitude', 'Fontsize', 20);
custom_tick=[90:90:360];
xticks(custom_tick)
subplot (3,1,3);
Fs1= 1000;
F1= 1;
twin=[0 1];
f_tvec= twin(1):1/Fs1:twin(2);
fake_phase=[0:360/Fs1:360];
fake_signal=sin(2*pi*F1*f_tvec);
plot(fake_phase,fake_signal,'Color',Burt_orange )
xlim([0 fake_phase(end)]);
ylabel('Theta amplitude');
xlabel('Theta phase(Deg)','Fontsize', 14);
xticks(custom_tick);
%Restriuct the signals at the time when the laser was on 

%Extract those 


