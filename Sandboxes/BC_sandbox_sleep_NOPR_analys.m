%% Sandbox for NOPR analsys

%% Initialize
% Go to the directory of some data
data_dir= '/Users/bryancontrerasmercado/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/CHRNA2_NOVEL_OBJECT/raw_data/NOPR/BC053_2023_11_16_D1_HAB_T2';
%data_dir= 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_NOVEL_OBJECT\raw_data\NOPR\BC053_2023_11_16_D1_HAB_T2'
cd(data_dir)
%% Parameters
emg_chan = 'CSC1.ncs';
lfp_chan = 'CSC2.ncs';
plot_flag = 01;

%% Loading some data

% Load the CSC guide
cfg = [];
cfg.fc = {lfp_chan}; 
csc = MS_LoadCSC(cfg); 

% Load the EMG
cfg_emg = [];
cfg_emg.fc = {emg_chan};
emg = MS_LoadCSC(cfg_emg); 

% load the events
evts = LoadEvents([]);

%% Restrict the data only to the sleep phase

%Extractinc the time stamps of the recording 
start_OF = evts.t{find(contains(evts.label, 'Starting Recording'))}(1); 
start_sleep = evts.t{find(contains(evts.label, 'Starting Recording'))}(2); 

end_OF = evts.t{find(contains(evts.label, 'Stopping Recording'))}(1); 
end_sleep = evts.t{find(contains(evts.label, 'Stopping Recording'))}(2); 

%Printing the duration of OF and sleeping 
fprintf('<strong>OF duration: %.2f mins</strong>\n', (end_OF - start_OF)/60); 
fprintf('<strong>Sleep duration: %.2f mins</strong>\n', (end_sleep - start_sleep)/60); 

%Restrict the data to just the sleep phase. 

csc_s = restrict(csc, start_sleep, end_sleep);
emg_s = restrict(emg, start_sleep, end_sleep);

%Correcting times

csc_s.tvec= csc_s.tvec-csc_s.tvec(1);
emg_s.tvec=emg_s.tvec- emg_s.tvec(1);

%% plot a bit of data for quality check
if plot_flag

    figure(1)
    clf

    ax(1) = subplot(2,1,1);
    plot(csc_s.tvec, csc_s.data);
    legend('HC LFP')

    ax(2) = subplot(2,1,2);
    plot(emg_s.tvec, emg_s.data - .001, 'r') % offset it a bit
    legend('emg')

    linkaxes(ax, 'x'); % locks the x axes of the subplots so if you zoom on one the other zooms in.

    % fix the annoying wide limits on the x axis
    xlim([csc_s.tvec(1) csc_s.tvec(end)])

end

%% sleep state

% specify known wake times or periods to ignore. 
wake_t = [0 4000 5015 5631 8824 9340 10605 10873 10618 108080 11531 11803 12450 13021 13618 13743 14188 14397]; 
%wake_t = [0 0]; 

wake_idx = nearest_idx(wake_t, csc_s.tvec);
wake_idx = reshape(wake_idx,2, length(wake_idx)/2)'; 
%% score the sleep. 
 
%% Getting the percentage of sleep sates 
figure(222)
clf
histogram(hypno.data)
