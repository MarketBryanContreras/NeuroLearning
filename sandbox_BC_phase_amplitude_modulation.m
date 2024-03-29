%% sandbox  BC laser pulse start-end
% Initialize the meta structure
Phase_amp_meta= struct('SG_Mod_idx', [],'FG_Mod_idx', [] ,'SG_NS_Phase_bins', [],'SG_S_Phase_bins', [],'FG_NS_Phase_bins', [],'FG_S_Phase_bins', [],'SG_z_scr',[],'FG_z_scr',[]);
%% Load Meta data structure
data_dir_inter = 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter';
cd (data_dir_inter); 
load('Phase_amp_meta.mat')
%% Go to some data
%BC011
% data_dir = "C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC011_09_11_2023_D4_INHIBITION";
% cd (data_dir); %csc3,
% i=4;goodcsc={'CSC3.ncs'};
%BC1602
% data_dir = "C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC1602_07_12_2023_D4_INHIBITION";
% cd (data_dir); %csc2,4,5
% i=3;goodcsc={'CSC4.ncs'};
% BC1807
data_dir = 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC1807_07_12_2023_D4_INHIBITION';
cd (data_dir); %csc3,4
i=2;goodcsc={'CSC4.ncs'};
% %BC051
% data_dir = 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC051_09_11_2023_D4_INHIBITION';
% cd (data_dir); %csc2,csc5
% i=1;goodcsc={'CSC2.ncs'};
% %BC053
% data_dir =  'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC053_30_10_2023_D4_INHIBITION';
% cd (data_dir); %csc2,csc5
% i=1;goodcsc={'CSC2.ncs'};
%% Assigning color pallete
Archt_green = [0.530 0.820 0.645];
Oxford_blue = [0.039 0.137 0.259];
Powder_blue = [0.682 0.772 0.921];
Burnt_orange= [0.758 0.348 0.249];
Red_crayola= [0.937 0.176 0.337];
Web_orange= [1.00 0.678 0.020];

Archt_green_alpha = [0.530 0.820 0.645 0.3];
Oxford_blue_alpha= [0.039, 0.137, 0.259, 0.300];
Powder_blue_alpha= [0.682 0.772 0.921 0.3];
Burnt_orange_alpha= [0.933 0.423 0.302 0.3];
Web_orange_alpha= [1.00 0.678 0.020 0.3];

%% Standard time that my mice are placed into the maze
time_maze_start=30;
%% Load Events and CSC
%Events
evts = LoadEvents([]); 
%CSC
cfg =[];
cfg.fc = goodcsc;
csc = MS_LoadCSC(cfg);
fs = csc.cfg.hdr{1}.SamplingFrequency;
sSampCsc= time_maze_start*fs; %This is the idx of the samplimg at sec 30
restrictIvCsc=iv(csc.tvec(sSampCsc),csc.tvec(end)); %Creates and IV from the time the mouse was placed in the maze to the end of the recording
csc= restrict(csc,restrictIvCsc);

%EMG
cfg_csc = [];
cfg_csc.fc = {'CSC1.ncs'};
emg = MS_LoadCSC(cfg_csc);
% corrected_time=csc.tvec-min(csc.tvec); % creates n array of time corrected for the time that neuralynx started recording
% Correction of time stamps
strt=csc.tvec(1);
csc.tvec= csc.tvec-strt; %correct time csc
nn=size(evts.t,2);
evts.t{1,nn-1}=evts.t{1,nn-1}-strt % Coreccting laser events start
evts.t{1,nn}=evts.t{1,nn}-strt % Coreccting laser events end
%% Get events of interest: inhb, running, noInhb
% Start with off events
pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).'; 
% pattern = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).'; 

iv_inhb = MS_get_evts_off(evts, pattern); 
% Conitnue with the running events
pos = MS_DLC2TSD(cd, [], [7.2 7.2]);
sSamp=time_maze_start*30;
restrictIvPos=iv(pos.tvec(sSamp),pos.tvec(end));
pos=restrict(pos,restrictIvPos);
pos.tvec= pos.tvec-pos.tvec(1);
x_data=pos.data(3,:);
%assign your running area 
max_tresh=75;
min_tresh=20;
x_running_index= find(x_data>=min_tresh & x_data<=max_tresh);
x_dif=diff(x_running_index);
z=find(x_dif>1);
y=z+1;
jumps_end=x_running_index(z);
jumps_end=[jumps_end,x_running_index(end)];
jumps_start=x_running_index(y);
jumps_start=[x_running_index(1), jumps_start];
% f_jumps=[jumps_start,jumps_end];
% sort(f_jumps);
% %plot thhe data
% plot(x_data)
% hold on
% scatter(jumps_end,x_data(jumps_end),'filled')
% scatter(jumps_start,x_data(jumps_start),'filled')
%Convert this indices to time epochs
jump_start_time=pos.tvec(jumps_start);
jump_end_time=pos.tvec(jumps_end);
iv_running=iv(jump_start_time, jump_end_time); % This is the interval where the mouse is running
% ----To do---- 
%1.Remove those intervals where the mouse is to close to the
%previous one
%
%2.Check for velocity
%---

% Substarct the iv where the mosue is inhibited from those where it is
% running
cfg01=[];
cfg01.verbose=0;
iv_noInhb=DifferenceIV(cfg01,iv_running,iv_inhb);
%% Plot the x_pos by time and plot the int to collaborate that you got them right
clf;
fig=figure;
x = pos.tvec;
y = x_data;
z = zeros(size(x));
speed_data=pos.data(5,:);
col = speed_data;  % This is the color, it varies with x in this case.
ax1=subplot(2,1,1)
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
c=colorbar;
c.Label.String='Velocity (cm/s)';
c.Label.FontSize=12;
c.Position= [0.910 0.5900 0.0100 0.330] 
xlim([pos.tvec(1) pos.tvec(end)]);
set(gca, 'TickDir', 'out');

ax1=subplot(2,1,1)

%xlabel('Time (s)');
%ylabel('Position in the x axis of the maze (cm)');
%h03=LTplotIvBars(iv_running,x_data,Oxford_blue,0.3)
%corrected_time=csc.tvec-min(csc.tvec); % creates n array of time corrected for the time that neuralynx started recording
ax2=subplot(2,1,2)
plot(pos.tvec,x_data, 'Color', Oxford_blue);
h01=LTplotIvBars(iv_inhb,x_data,Archt_green,0.8)
h02=LTplotIvBars(iv_noInhb,x_data,Burnt_orange,0.4)
xlim([pos.tvec(1) pos.tvec(end)]);
ylim([0 100]);
% xlabel('Time (s)','Fontsize',14);
% ylabel('Position in the maze (cm)');
leg=legend([h01;h02],{'Laser 520nm','Running no inhibition'});
legendFontSize = 12; % Adjust the font size as needed
leg.Position=[0.850 0.3914 0.1047 0.0487];
set(leg, 'FontSize', legendFontSize);
box off;
legend boxoff ;
set(gca, 'TickDir', 'out');

%Give common xlabel, ylabel, and title
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
han.YLabel.Position=[-0.0300 0.5000 0];
ylabel(han,'Position in the x axis of the maze (cm)','FontSize', 16);
han.XLabel.Position=[0.5000 -0.060 0];
xlabel(han,'Time(s)','FontSize', 16);
title(han,'Example of displacement of mouse in the x axis during inhibition session','FontSize', 20);
fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1]; 
hold off;
% saveas(gcf, 'XPos_time.png');
% % The function fill to create rentangles use the arguments in the following fashion:fill([xmin xmax xmax xmin], [ymin ymin ymax ymax],Burt_orange, 'FaceAlpha', 0.3, 'LineStyle',"none")
%% Filtering in the required bands

% filter the LFP in the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [5 12]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 3; %type filter order
cfg_filt_t.display_filter = 0; % use this to see the fvtool
theta_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
% Filter the LFP in the slow gamma band
cfg_filt_t = [];
cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [30 58]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 4; %type filter order
cfg_filt_t.display_filter = 0; % use this to see the fvtool
SG_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
% Filter the LFP in the fast gamma band
cfg_filt_t = [];
cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [60 100]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 4; %type filter order
cfg_filt_t.display_filter = 0; % use this to see the fvtool
FG_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
% Filter the LFP in the 1-50 band
cfg_filt_t = [];
cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [1 50]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 4; %type filter order
cfg_filt_t.display_filter = 0; % use this to see the fvtool
ref_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using
%% Restricting data to intervals of inhb, noInhb
%inhb
csc_inhb = restrict(csc, iv_inhb);

theta_inhb=restrict(theta_csc, iv_inhb);
SG_inhb = restrict(SG_csc, iv_inhb);
FG_inhb = restrict(FG_csc, iv_inhb);
ref_inhb= restrict(ref_csc, iv_inhb);
%no-inhb
csc_noinhb = restrict(csc, iv_noInhb);

theta_noinhb = restrict(theta_csc, iv_noInhb);
SG_noinhb = restrict(SG_csc, iv_noInhb);
FG_noinhb = restrict(FG_csc, iv_noInhb);
ref_noinhb = restrict(ref_csc, iv_noInhb);
%running
theta_running = restrict(theta_csc, iv_running);
SG_running = restrict(SG_csc, iv_running);
FG_running = restrict(FG_csc, iv_running);
ref_running = restrict(ref_csc, iv_running);
%% Power analysis
%Obtaining power 
% Global
%   Inhb
ThetaInhbPower= BC_power(theta_inhb);
SGInhbPower= BC_power(SG_inhb);
FGInhbPower= BC_power(FG_inhb);
refInhbPower= BC_power(ref_inhb);
%    Normalizing to the power of reference
Thetainhb_norm= ThetaInhbPower./refInhbPower;
SGinhb_norm= SGInhbPower./refInhbPower;
FGinhb_norm= FGInhbPower./refInhbPower;
%   NoInhb
ThetaNoInhbPower= BC_power(theta_noinhb);
SGNoInhbPower= BC_power(SG_noinhb);
FGNoInhbPower= BC_power(FG_noinhb);
refNoInhbPower= BC_power(ref_noinhb);
%    Normalizing to the power of reference
ThetaNoinhb_norm= ThetaNoInhbPower./refNoInhbPower;
SGNoinhb_norm= SGNoInhbPower./refNoInhbPower;
FGNoinhb_norm= FGNoInhbPower./refNoInhbPower;

%    Obtaining the sum
%     Ihnb
sum_TInhbPower=sum(ThetaInhbPower);
sumSGInhbPower=sum(SGInhbPower);
sumFGInhbPower=sum(FGInhbPower);
%      Inhb normalized
sum_TInhbPowerNorm=sum(Thetainhb_norm);
sumSGInhbPowerNorm=sum(SGinhb_norm);
sumFGInhbPowerNorm=sum(FGinhb_norm);
%     NoIhnb
sum_TNoInhbPower=sum(ThetaNoInhbPower);
sumSGNoInhbPower=sum(SGNoInhbPower);
sumFGNoInhbPower=sum(FGNoInhbPower);
%      Inhb normalized
sum_TNoInhbPowerNorm=sum(ThetaNoinhb_norm);
sumSGNoInhbPowerNorm=sum(SGNoinhb_norm);
sumFGNoInhbPowerNorm=sum(FGNoinhb_norm);

%    Obtaining the mean
%     Ihnb
mean_TInhbPower=mean(ThetaInhbPower);
meanSGInhbPower=mean(SGInhbPower);
meanFGInhbPower=mean(FGInhbPower);
%      Inhb normalized
mean_TInhbPowerNorm=mean(Thetainhb_norm);
meanSGInhbPowerNorm=mean(SGinhb_norm);
meanFGInhbPowerNorm=mean(FGinhb_norm);
%     NoIhnb
mean_TNoInhbPower=mean(ThetaNoInhbPower);
meanSGNoInhbPower=mean(SGNoInhbPower);
meanFGNoInhbPower=mean(FGNoInhbPower);
%      Inhb normalized
mean_TNoInhbPowerNorm=mean(ThetaNoinhb_norm);
menSGNoInhbPowerNorm=mean(SGNoinhb_norm);
meanFGNoInhbPowerNorm=mean(FGNoinhb_norm);
%% Power analysis
%Obtaining power of each running epoch
% Per inhb epoch
sum_powers=[];
mean_powers=[];
for ie=1:length(iv_inhb.tstart)
    s=iv_inhb.tstart(ie);e=iv_inhb.tend(ie);
        if e-s<=9.9e-05
         sum_powers(ie,1)= NaN;
         sum_powers(ie,2)= NaN;
         sum_powers(ie,3)= NaN;
        continue
    end
    iv_momentary= iv(s,e);
    sum_powers(ie,1)= sum(BC_power(restrict(theta_csc, iv_momentary)));
    sum_powers(ie,2)= sum(BC_power(restrict(SG_csc, iv_momentary)));
    sum_powers(ie,3)= sum(BC_power(restrict(FG_csc, iv_momentary)));
end 
for ie=1:length(iv_noInhb.tstart)
    s=iv_noInhb.tstart(ie);e=iv_noInhb.tend(ie);
    iv_momentary= iv(s,e);
    if e-s<=9.9e-05
         sum_powers(ie,4)= NaN;
         sum_powers(ie,5)= NaN;
         sum_powers(ie,6)= NaN;
        continue
    end
    sum_powers(ie,4)= sum(BC_power(restrict(theta_csc, iv_momentary)));
    sum_powers(ie,5)= sum(BC_power(restrict(SG_csc, iv_momentary)));
    sum_powers(ie,6)= sum(BC_power(restrict(FG_csc, iv_momentary)));
end 
%%Converting 0 to NaN values
sum_powers(sum_powers==0)=NaN;
%%Determining significance
for ai=1:(size(sum_powers,2)/2)
    inhb=sum_powers(:,ai);noinhb=sum_powers(:,ai+3);
    bands={'theta','SG', 'FG'};
    if ttest(inhb, noinhb)
        fprintf('\n This mouse had a significance differnce in the %s band when looking at the sum\n', bands{ai});
    else
        fprintf('\n This mouse DID NOT had a significance differnce in the %s band when looking at the sum\n', bands{ai});
   end
end
for ie=1:length(iv_inhb.tstart)
    s=iv_inhb.tstart(ie);e=iv_inhb.tend(ie);
        if e-s<=9.9e-05
         mean_powers(ie,1)= NaN;
         mean_powers(ie,2)= NaN;
         mean_powers(ie,3)= NaN;
        continue
    end
    iv_momentary= iv(s,e);
    mean_powers(ie,1)= mean(BC_power(restrict(theta_csc, iv_momentary)));
    mean_powers(ie,2)= mean(BC_power(restrict(SG_csc, iv_momentary)));
    mean_powers(ie,3)= mean(BC_power(restrict(FG_csc, iv_momentary)));
end 
for ie=1:length(iv_noInhb.tstart)
    s=iv_noInhb.tstart(ie);e=iv_noInhb.tend(ie);
    iv_momentary= iv(s,e);
    if e-s<=9.9e-05
         mean_powers(ie,4)= NaN;
         mean_powers(ie,5)= NaN;
         mean_powers(ie,6)= NaN;
        continue
    end
    mean_powers(ie,4)= mean(BC_power(restrict(theta_csc, iv_momentary)));
    mean_powers(ie,5)= mean(BC_power(restrict(SG_csc, iv_momentary)));
    mean_powers(ie,6)= mean(BC_power(restrict(FG_csc, iv_momentary)));
end 
%%Converting 0 to NaN values
mean_powers(mean_powers==0)=NaN;
%%Determining significance
for ai=1:(size(mean_powers,2)/2)
    inhb2=mean_powers(:,ai);noinhb2=mean_powers(:,ai+3);
    bands={'theta','SG', 'FG'};
    if ttest(inhb2, noinhb2)
        fprintf('\n This mouse had a significance differnce in the %s band when looking at the mean\n', bands{ai});
    else
        fprintf('\n This mouse DID NOT had a significance differnce in the %s band when looking at the mean\n', bands{ai});
   end
end
%% Calculating the mod idx theta-slow gamma for intevals with inhibition and no inhibition
%Inhb
SGphiAmpNormInhb=BC_phase_amp_norm_bins(theta_inhb,SG_inhb);
SGInhb_modidx=MS_ModIdx(SGphiAmpNormInhb);
BC_plot_modidx(SGphiAmpNormInhb,Archt_green,SGInhb_modidx,'SG', 'silencing' )
%NoInhb
SGphiAmpNormNoInhb = BC_phase_amp_norm_bins(theta_noinhb,SG_noinhb);
SGNoInhb_modidx=MS_ModIdx(SGphiAmpNormNoInhb);
BC_plot_modidx(SGphiAmpNormNoInhb,Powder_blue,SGNoInhb_modidx, 'SG', 'no silencing')
% %Norestricted
% SGphiAmpNorm = BC_phase_amp_norm_bins(theta_csc,SG_csc);
% SGAll_modidx=MS_ModIdx(SGphiAmpNorm);
% BC_plot_modidx(SGphiAmpNorm)
% %All_running
% SGphiAmpNormRunning = BC_phase_amp_norm_bins(theta_running,SG_running);
% SGRunning_modidx=MS_ModIdx(SGphiAmpNormRunning);
% BC_plot_modidx(SGphiAmpNormRunning)
%Shifted
[SGshift_mean,SGshift_std]=LTshifted_meanModIdx(theta_running,SG_running);

% %% Lets make a plot with the mod idx 
% genmodidx=[ SGInhb_modidx,SGNoInhb_modidx];
% 
% h=bar(genmodidx, 'FaceColor','flat', 'EdgeColor', 'none');
% h.CData(1,:) = Archt_green;h.CData(2,:) = Oxford_blue;
% xticklabels({'Silencing', 'No silencing'});
% grid on
% ylabel('Mod idx')
% set(gca, 'TickDir', 'out');  % Move ticks outside the plot
% box off;                     % Turn off the box around the plot
% 
% % Add the results to the meta structure
% Phase_amp_meta.SG_Mod_idx(i,1)=SGInhb_modidx;Phase_amp_meta.SG_Mod_idx(i,2)=SGNoInhb_modidx;
% Phase_amp_meta.SG_S_Phase_bins(i,:)=SGphiAmpNormInhb;Phase_amp_meta.SG_NS_Phase_bins(i,:)=SGphiAmpNormNoInhb;
%% Calculating the mod idx theta-fast gamma for intevals with inhibition and no inhibition
%Inhb
FGphiAmpNormInhb=BC_phase_amp_norm_bins(theta_inhb,FG_inhb);
FGInhb_modidx=MS_ModIdx(FGphiAmpNormInhb);
BC_plot_modidx(FGphiAmpNormInhb)
%NoInhb
FGphiAmpNormNoInhb = BC_phase_amp_norm_bins(theta_noinhb,FG_noinhb);
FGNoInhb_modidx=MS_ModIdx(FGphiAmpNormNoInhb);
BC_plot_modidx(FGphiAmpNormNoInhb)
% %Norestricted
% FGphiAmpNorm = BC_phase_amp_norm_bins(theta_csc,FG_csc);
% FGAll_modidx=MS_ModIdx(FGphiAmpNorm);
% % BC_plot_modidx(FGphiAmpNorm)
% %All_running
% FGphiAmpNormRunning = BC_phase_amp_norm_bins(theta_running,FG_running);
% FGRunning_modidx=MS_ModIdx(FGphiAmpNormRunning);
% BC_plot_modidx(FGphiAmpNormRunning)
%Shifted
[FGshift_mean,FGshift_std]=LTshifted_meanModIdx(theta_running,FG_running);
% Add the results to the meta structure
% Phase_amp_meta.FG_Mod_idx(i,1)=FGInhb_modidx;Phase_amp_meta.FG_Mod_idx(i,2)=FGNoInhb_modidx;
% Phase_amp_meta.FG_S_Phase_bins(i,:)=FGphiAmpNormInhb;Phase_amp_meta.FG_NS_Phase_bins(i,:)=FGphiAmpNormNoInhb;
%% Plotting and t-test
%Lets check if our data past the t-test
[h, p] = ttest(Phase_amp_meta.SG_z_scr(1:3,1), Phase_amp_meta.SG_z_scr(1:3,2));
%Fast_gamma
x=[Phase_amp_meta.FG_z_scr(1:3,1),Phase_amp_meta.FG_z_scr(1:3,2)];
boxplot(x,'Labels',{'Silencing', 'No_silencing'});
title('Normalized FG ModIdx ArchT');
ylabel('Z-Scores');           

%Slow_gamma
y=[Phase_amp_meta.SG_z_scr(1:3,1),Phase_amp_meta.SG_z_scr(1:3,2)];
boxplot(y,'Labels',{'Silencing', 'No_silencing'})
title('Normalized SG ModIdx ArchT')
ylabel('Z-Scores');   

colors = [Burnt_orange;  % RGB values for Group 2
          Archt_green]  % RGB values for Group1;

h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for i = 1:numel(h)
    patch(get(h(i), 'XData'), get(h(i), 'YData'), colors(i, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot
grid on;
% Adjust figure properties
fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1];         % Set background color to white
% Adjust font properties
set(gca, 'FontName', 'Helvetica', 'FontSize', 12);  % Change font name and size
% Resize the figure (optional)
fig.Position = [100, 100, 800, 500];  % [x, y, width, height]
%saveas(gcf, 'FG_NormModIdx.png');
%% Normalizing ModIdx
%Slow gamma
SG_z_score_inhb= (SGInhb_modidx-SGshift_mean)/SGshift_std;
SG_z_score_NoInhb= (SGNoInhb_modidx-SGshift_mean)/SGshift_std;
Phase_amp_meta.SG_z_scr(i,1)=SG_z_score_inhb;Phase_amp_meta.SG_z_scr(i,2)=SG_z_score_NoInhb;
%Fast gamma
FG_z_score_inhb= (FGInhb_modidx-FGshift_mean)/FGshift_std;
FG_z_score_NoInhb= (FGNoInhb_modidx-FGshift_mean)/FGshift_std;
Phase_amp_meta.FG_z_scr(i,1)=FG_z_score_inhb;Phase_amp_meta.FG_z_scr(i,2)=FG_z_score_NoInhb;
%% Analizing the theta power 
%Create a spectogram to analize the power of different frequencies during
%the whole session
Fs= csc.cfg.hdr{1}.SamplingFrequency;
[S,F,T,P]=spectrogram(csc.data, hanning(512),256,1:200,Fs);
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
cor_iv_inhb=iv(iv_inhb.tstart-(csc.tvec(1)), iv_inhb.tend-(csc.tvec(1)));
cor_iv_NoInhb=iv(iv_noInhb.tstart-(csc.tvec(1)), iv_noInhb.tend-(csc.tvec(1)));
h01=LTplotIvBars(cor_iv_inhb,[1:200],Oxford_blue,0.2)
h02=LTplotIvBars(cor_iv_NoInhb,[1:200],Burnt_orange,0.2)
%Add soemthing that allows you to identify the periods of inh and no inhb

%Gather the average theta power of inhb and noinhb 

%% Saving values
zScores.data.FG.Inhb=[zScores.data.FG.Inhb,FG_z_score_inhb];
zScores.data.FG.NoInhb=[zScores.data.FG.NoInhb,FG_z_score_NoInhb];

zScores.data.SG.Inhb=[zScores.data.SG.Inhb,SG_z_score_inhb];
zScores.data.SG.NoInhb=[zScores.data.SG.NoInhb,SG_z_score_NoInhb];
%% Plotting
%Fast_gamma
x=[zScores.data.FG.Inhb',zScores.data.FG.NoInhb'];
boxplot(x,'Labels',zScores.label);
title('Z-Scores of FG ModIdx');
ylabel('Normalized ModIdx');           

%Slow_gamma
y=[zScores.data.SG.Inhb',zScores.data.SG.NoInhb']
boxplot(y,'Labels',zScores.label)
title('Z-Scores of SG ModIdx')
colors = [Burnt_orange;  % RGB values for Group 2
          Archt_green]  % RGB values for Group1;

h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for i = 1:numel(h)
    patch(get(h(i), 'XData'), get(h(i), 'YData'), colors(i, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot
grid on;
% Adjust figure properties
fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1];         % Set background color to white
% Adjust font properties
set(gca, 'FontName', 'Helvetica', 'FontSize', 12);  % Change font name and size
% Resize the figure (optional)
fig.Position = [100, 100, 800, 500];  % [x, y, width, height]
%saveas(gcf, 'FG_NormModIdx.png');
%% Theta power analysis 
win_size = 2^9; % works out to 512 samples. Out of some superstition I always use base 2 (b/c of bytes or something) when computing spectra. 
n_overlap = win_size/4; % just needs to be smaller than the window. 1/4 gives nice temporal resolution. 
fs = csc.cfg.hdr{1}.SamplingFrequency;
freq_range = 1:0.25:120; % range of frequencies. 

[~,F,T, P] = spectrogram(csc.data,win_size,n_overlap,freq_range,fs); 

T = T+csc.tvec(1); % just so that our time vectors are the same

subplot(3,1,2)
cla
imagesc(T, F, 10*log10(P))
set(gca, 'YDir', 'normal')
hold on

for ii = 1:length(laser_on)
    rectangle('position', [laser_on(ii), 0, laser_off(ii) - laser_on(ii), 120], 'FaceColor', 'none', 'EdgeColor', 'g')

end

%% 
% %on_ts2 = evts.t{(strcmp('TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).',evts.label))};
% on_idx = find(contains(evts.label, '(0x0002)'));  % find the correct label. 
% on_ts = evts.t{on_idx}; % grab those times
% off_idx = find(contains(evts.label, '(0x0000)'));  % find the corresponding off label
% off_ts = evts.t{off_idx}; 
% 
% % convert to 'iv' (inverval) format for simplicity
% laser_iv = iv(on_ts, off_ts); 
% %% plot some examples
% figure(101)
% clf
% hold on
% 
% max_csc = max(csc.data); 
% min_csc = min(csc.data); 
% offset = (max_csc - min_csc);%I use these lines to calculate the y limits (+- 10% of these value) 
% range_final= [min_csc-offset*0.1 max_csc+offset*0.1];
% height=(max_csc+offset*0.1)-(min_csc-offset*0.1);
% 
% for ii = 1:length(laser_iv.tstart)
%     
%     h=fill([(laser_iv.tstart(ii)-min(csc.tvec)) (laser_iv.tend(ii)-min(csc.tvec)) (laser_iv.tend(ii)-min(csc.tvec)) (laser_iv.tstart(ii)-min(csc.tvec))], [range_final(1) range_final(1) range_final(2) range_final(2)],Archt_green, 'FaceAlpha', 0.3, 'LineStyle',"none");
%    
% end
% 
% % for ii = 1:length(running_iv.tstart)
% %     
% %     h2=fill([(running_iv.tstart(ii)-min(csc.tvec)) (running_iv.tend(ii)-min(csc.tvec)) (running_iv.tend(ii)-min(csc.tvec)) (running_iv.tstart(ii)-min(csc.tvec))], [range_final(1) range_final(1) range_final(2) range_final(2)],Burnt_orange, 'FaceAlpha', 0.3, 'LineStyle',"none");
% %    
% % end
% % The function fill to create rentangles use the arguments in the following fashion:fill([xmin xmax xmax xmin], [ymin ymin ymax ymax],Burt_orange, 'FaceAlpha', 0.3, 'LineStyle',"none")
% %Adventage of fill over rentanlge is that it can be used in the label of the figure 
% % plot the data on top
% corrected_time=csc.tvec-min(csc.tvec); % creates n array of time corrected for the time that neuralynx started recording
% plot(corrected_time, csc.data, 'Color', Oxford_blue );
% xlim([corrected_time(1) corrected_time(end)]);
% ylim(range_final);
% xlabel('Time (s)');
% ylabel('Potential (�V)');
% title('Periods of inhibition in linear track', 'Fontsize', 20);
% lasleg=legend(h, 'Laser 520nm');
% 
% legendFontSize = 14; % Adjust the font size as needed
% set(lasleg, 'FontSize', legendFontSize);
% legend boxoff ;
% %print('-dpng', 'Inhibition_epochs.png', '-r300'); %code to save the figure
% 
% %with 300 dpi


% theta_amp = abs(hilbert(theta_csc.data)); % get the amplitude
% theta_phi  = angle(hilbert(theta_csc.data(1,:))); %get the phase
% %theta_csc.data = theta_csc.data(1,:); 
%% Filter on slow gamma


% SG_amp = abs(hilbert(SG_csc.data)); % get the amplitud
% 




% FG_amp = abs(hilbert(FG_csc.data)); % get the amplitude for Fast Gamm


%% Try some Phase-amp coupling

mod_th_g = MS_ModIdx_win(theta_csc, SG_csc, 30*theta_csc.cfg.hdr{1}.SamplingFrequency, iv_inhb);
title('Theta - SG Mod Index')
%print('-dpng', 'Thetha_slow_gamma_modulation_index.png', '-r300'); %code to save the figure

mod_th_fg = MS_ModIdx_win(theta_csc, FG_csc, 30*theta_csc.cfg.hdr{1}.SamplingFrequency, iv_inhb);
title('Theta - FG Mod Index')
%print('-dpng', 'Thetha_fast_gamma_modulation_index.png', '-r300'); %code to save the figure
%% Plot only the modulation indexes
figure(101)
clf
plot(theta_csc.tvec, mod_th_g, 'Color', Powder_blue); 
hold on
plot(theta_csc.tvec, mod_th_fg, 'Color', Red_crayola); 

h01=LTplotIvBars(iv_inhb,mod_th_fg,Archt_green,0.6)
h02=LTplotIvBars(iv_noInhb,mod_th_fg,Burnt_orange,0.4)
%fill([xmin xmax xmax xmin], [ymin ymin ymax ymax],Burt_orange, 'FaceAlpha', 0.3, 'LineStyle',"none")
ylabel('Mod Idx' )
xlabel('Time(s)')
legend({'Theta - SG', 'Theta - FG', 'Laser 520nm', 'No Inhb running'})
xlim([theta_csc.tvec(1) theta_csc.tvec(end)])
ylim([0 Mod_range_final(2)])
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
% 
% 
% SG_FG_ratio = smooth(SG_amp, 2*theta_csc.cfg.hdr{1}.SamplingFrequency)./smooth(FG_amp, 2*theta_csc.cfg.hdr{1}.SamplingFrequency);
% 
% figure(10)
% hold on
% plot(csc.tvec, SG_FG_ratio)
% 
% 
% % c_ord = linspecer(length(on_ts));
% Artch_green = [0.530 0.820 0.645];
% for ii = 1:length(laser_iv.tstart)
%     
%     rectangle('position', [laser_iv.tstart(ii), min_csc, laser_iv.tend(ii) - laser_iv.tstart(ii), 100], 'facecolor', [Artch_green .2], 'edgecolor', [Artch_green .2])
%     
% end
%% Creating a Phase-Amplitud plot for the times when the laser was on 
%theta_amp;
%theta_phi;
%1. Restrict the data to the intervals where you had the laser on 
csc_inhb = restrict(csc, iv_inhb);

theta_ihnb = restrict(theta_csc, iv_inhb);
SG_inhb = restrict(SG_csc, iv_inhb);
FG_inhb = restrict(FG_csc, iv_inhb);

thta_amp_laser= abs(hilbert(theta_laser.data));
theta_phi_laser= angle(hilbert(theta_laser.data(1,:)));

SG_amp_laser= abs(hilbert(SG_laser.data));
SG_phi_laser= angle(hilbert(SG_laser.data(1,:)));
%2. create a composite of the phase_theta, amplitud_SG
phi_bins = -pi:pi/18:pi; 
[~,theta_idx]= histc(theta_phi_laser, phi_bins); %This creates a vector with the  cooreponding bin of each phase
%3. The phases are binned according to which phase they belong to and their
%mean is calculated
phase_means= zeros(1,length(unique(theta_idx)));
for ii= 1: length(unique(theta_idx));
    phase_means(ii)= nanmean(SG_amp_laser(theta_idx == ii));
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
plot(fake_phase,fake_signal,'Color',Burnt_orange )
xlim([0 fake_phase(end)]);
ylabel('Theta amplitude');
xlabel('Theta phase(Deg)','Fontsize', 14);
xticks(custom_tick);

%% For the full session
phase_means_norm= BC_phase_amp_norm_bins(theta_csc,SG_csc);
%5. Calculating the ModIdx
% mod_idx_real= MS_ModIdx(phase_means_norm);
% %6. Shifting data
% nshuffles=500;
% for si=nshuffles:-1:1;
%         shift_phi= circshift(theta_phi,randi(length(theta_phi)-1));
%         shift_MI(si)= MS_Modidx(shift_phi,SG_amp);
% phase_means_shift= zeros(1,length(unique(theta_idx)));
% for ii= 1: length(unique(theta_idx));
%     phase_means_shift(ii)= nanmean(SG_shift_amp(theta_idx == ii));
% end
% phase_means_shift_norm= phase_means_shift/sum(phase_means_shift);
% mod_idx_shift= MS_ModIdx(phase_means_norm);

%5 Plot the phase amplitude coupling
% figure
% clf
% ax(1) = subplot(3,1,1:2);
% deg_bins= (phi_bins(2:end)+pi)*(180/pi); %Here I am adding pi radinas to correct for the hilber transformation and converting to bins in degrees
% bar(deg_bins,phase_means_norm,'FaceColor',Oxford_blue, 'EdgeColor',Oxford_blue )
% %xlabel('Theta phase(Deg)');
% ylabel('SG amplitude');
% title('Theta-phase-SG-amplitude', 'Fontsize', 20);
% custom_tick=[90:90:360];
% xticks(custom_tick)
% subplot (3,1,3);
% Fs1= 1000;
% F1= 1;
% twin=[0 1];
% f_tvec= twin(1):1/Fs1:twin(2);
% fake_phase=[0:360/Fs1:360];
% fake_signal=sin(2*pi*F1*f_tvec);
% plot(fake_phase,fake_signal,'Color',Burt_orange )
% xlim([0 fake_phase(end)]);
% ylabel('Theta amplitude');
% xlabel('Theta phase(Deg)','Fontsize', 14);
% xticks(custom_tick);
%Restriuct the signals at the time when the laser was on 
BC_plot_modidx(phase_means_norm);


%% Converting to a function
pos = MS_DLC2TSD(cd, [], [7.2 7.2]);
time_1st_end=30;
sSamp=time_1st_end*30;

pos=restrict(pos,sSamp,length(pos.tvec));
x_data=pos.data(3,:);
x_running_index= find(x_data>=20 & x_data<=75);
x_dif=diff(x_running_index);


z=find(x_dif>1);
y=z+1;
jumps_end=x_running_index(z);
jumps_end=[jumps_end,x_running_index(end)];
jumps_start=x_running_index(y);
jumps_start=[x_running_index(1), jumps_start];
f_jumps=[jumps_start,jumps_end];
sort(f_jumps);
% %plot thhe data
% plot(x_data)
% hold on
% scatter(jumps_end,x_data(jumps_end),'filled')
% scatter(jumps_start,x_data(jumps_start),'filled')

%Convert this indices to time epocs
jump_start_time=pos.tvec(jumps_start);
jump_end_time=pos.tvec(jumps_end);
running_iv=iv(jump_start_time, jump_end_time);

running_lfp= restrict(csc,running_iv.tstart,running_iv.tend)

