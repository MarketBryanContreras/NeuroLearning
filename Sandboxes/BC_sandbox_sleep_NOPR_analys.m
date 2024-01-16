%% Sandbox for NOPR analsys

%% Initialize
% Go to the directory of some data
%BC053 D1
%data_dir= '/Users/bryancontrerasmercado/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/CHRNA2_NOVEL_OBJECT/raw_data/NOPR/BC053_2023_11_16_D1_HAB_T2';
%BC053 D2
%data_dir= '/Users/bryancontrerasmercado/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/CHRNA2_NOVEL_OBJECT/raw_data/NOPR/BC053_2023_11_16_D1_HAB_T2';
%mac
data_dir= 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_NOVEL_OBJECT\raw_data\NOPR\BC053_2023_11_16_D1_HAB_T2';
%'windows
cd(data_dir)

%% Colloecting subject info
folder = pwd;
parts= split(folder, '/');
parts=parts{end};
parts= split(parts,'_');

info.subject=parts{1};
info.date=[ parts{2} '_' parts{3} '_' parts{4}];
info.session=parts{5};
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

% specify known wake times or periods to ignore. Note; The perios are in
% the interval format, so you want to write when do they start and when do
% they end in pairs
if strcmpi(info.subject,"BC053") %Specify the wake period according to subject and session
    wake_t = [0 3093 5015 5631 8824 9340 10605 10880 11531 11803 12450 13020 13618 13743 14188 14397]; 
end 
%wake_t = [0 0]; 
wake_idx = nearest_idx(wake_t, csc_s.tvec); %Converts time to sample idx
wake_idx = reshape(wake_idx,2, length(wake_idx)/2)'; 
%% score the sleep. 

[hypno, csc_out, emg_out] = dSub_Sleep_screener(csc_s, emg_s, wake_idx);  % can add in 'wake_idx' as the last input. 

%% Getting the percentage of sleep sates
figure(222)
clf
[y,x]=histcounts(hypno.data,[0.5:1:3.5]);
y_per=(y/sum(y))*100;
if plot_flag
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
    colororder reef
end
%to do add a label of the total time that the mice spent sleeping