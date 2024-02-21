%% Sandbox for NOPR analsys

%% Initialize

%% I have to create an automatic loader 
% Go to the directory of some data
%BC053 D1
%data_dir= '/Users/bryancontrerasmercado/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/CHRNA2_NOVEL_OBJECT/raw_data/NOPR/BC053_2023_11_16_D1_HAB_T2';
%BC053 D2
%data_dir= '/Users/bryancontrerasmercado/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/CHRNA2_NOVEL_OBJECT/raw_data/NOPR/BC053_2023_11_16_D1_HAB_T2';
data_dir='/Users/bryancontrerasmercado/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/CHRNA2_NOVEL_OBJECT/raw_data/Behavior/BC1807_2023_07_07_D2_NOPR';
%mac
%data_dir= 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_NOVEL_OBJECT\raw_data\NOPR\BC053_2023_11_16_D1_HAB_T2';
%data_dir= 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_NOVEL_OBJECT\raw_data\Behavior\BC1807_2023_07_07_D2_NOPR'
%'windows

cd(data_dir)

%% Colloecting subject info
% ----TO DO--- Adjust for MAc or WIN 
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
    new_colors= [...
        0.3467    0.5360    0.6907;
        0.9153    0.2816    0.2878;
        0.4416    0.7490    0.4322];
    colororder(newcolors);
    sgtitle(sprintf('Sleep phases %s', info.subject), 'fontweight', 'bold', 'fontsize', 16);
end
%to do add a label of the total time that the mice spent sleeping
%% Traacking movement and setting up intervals when it explored the objects
% Load the DLC data
pos = MS_DLC2TSD(cd, [], [4.5 4.5]); %Need to ask Eric how to solve for the size of the cage. For now we will assume its 4.5

%Correct time and remove time where the mouse was not in OF
pos.tvec=pos.tvec-pos.tvec(1);
%Correct positions of objects to the mean of the location for both objects
for ii= 5:1:8
    pos.data(ii,:)=mean(pos.data(ii,:));
end
%Calculate the distance between the objects and the mouse nose
%A. Calculate the distance between the nose and object A (Object A is the
%one on the left)
distances.data=[];
for iframe= length(pos.data(1,:)):-1:1
    x1=pos.data(1,iframe);
    y1=pos.data(2,iframe);
    x2=pos.data(5,iframe);
    y2=pos.data(6,iframe);
    d=sqrt(((x2-x1)^2)+((y2-y1)^2));
    distances.data(1,iframe)=d;
    if iframe==1
        clear iframe;clear x1;clear x2;clear y1;clear y2;clear d;
    end
end

%B.Same for object B (Object B is the one on the right)
for iframe= length(pos.data(1,:)):-1:1
    x1=pos.data(1,iframe);
    y1=pos.data(2,iframe);
    x2=pos.data(7,iframe);
    y2=pos.data(8,iframe);
    d=sqrt(((x2-x1)^2)+((y2-y1)^2));
    distances.data(2,iframe)=d;
end

%% Plot the position of the mouse
%parameters
video_flag=0;
%Initialize figure
figure (20)
%Actual position
subplot(8,1,1:4);
%plot the actual position 
time_frames=[1:length(pos.data)];
[aX,aY]=BC_Circle_plot(5.5,pos.data(5,1),pos.data(6,1), BC_color_genertor('Swamp_green'),plot_flag);
hold on
[bX,bY]=BC_Circle_plot(5.5,pos.data(7,1),pos.data(8,1), BC_color_genertor('Torment_blue'),plot_flag);

if plot_flag==1
scatter(pos.data(5,1),pos.data(6,1) , 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'none')
scatter(pos.data(7,1),pos.data(8,1) , 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
xlim([min(pos.data(1,:)) max(pos.data(1,:))]);
ylim([min(pos.data(2,:)) max(pos.data(2,:))]);

if video_flag==0;
    scatter(pos.data(1,:), pos.data(2,:), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.01);
else
    %%Trying the video
    for ii=1400:-1:900
        scatter(pos.data(1,ii), pos.data(2,ii), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.02);
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
%%% Rest of the plot
%plot X position of the mouse
subplot(8,1,5);
plot(time_frames,pos.data(1,:));
xlim([0,time_frames(end)])
%xlabel('Time(s)');
ylabel('X pos (cm)');
%plot Y position of the mouse
subplot(8,1,6);
plot(time_frames,pos.data(2,:), 'Color',BC_color_genertor('Oxford_blue'));
xlim([0,time_frames(end)])
%xlabel('Time(s)');
ylabel('Y pos (cm)');
%Plot the distances to object 1
subplot(8,1,7);
plot(time_frames, distances.data(1,:), 'm');
yline(5, 'Color', 'r')
xlim([0,time_frames(end)])
%xlabel('Time(s)');
ylabel('Distance to object A');
%Plot the distnaces to object 2
subplot(8,1,8);
plot(time_frames,distances.data(2,:),'r');
yline(5, 'Color', 'r')
xlim([0,time_frames(end)])
xlabel('Time(s)');
ylabel('Distance to object B');

hold off
%% Lets assign the intervals where te mice is inisde the radious
%Parameters
minFrames=5; % The minimum number of frames where the mouse is in radious

%Calculate the angle of the mouse head and if it is pointing towards the object

%Figure out intervals where the mouse nose was close to the object and the
%head direction is pointing towards the closest object
%Figure out whih data points are inside the radious
inA= inpolygon(pos.data(1,:),pos.data(2,:),aX,aY);
aidx=find(inA);
aInt=BC_object_interaction_intervals(aidx,minFrames);%Get periods where the mouse spends more than 5 frames (0.33 sec)
timeA=(sum(diff(aInt,1,2))/30);
inB= inpolygon(pos.data(1,:),pos.data(2,:),bX,bY);
bidx=find(inB);
bInt=BC_object_interaction_intervals(bidx,minFrames);
timeB=(sum(diff(bInt,1,2))/30);
fprintf('The mouse <strong>%s</strong> interacted <strong>%d</strong> times with the object 1 for a total of <strong>%.2f</strong> seconds and <strong>%d</strong> times with the object 2 for a total of <strong>%.2f</strong> seconds \n ', info.subject ,size(aInt,1),timeA,size(bInt,1), timeB)

%Lets verify that these idx correspond to times the mouse was with the objects
for jj=1:1:size(aInt,1)
    scatter(pos.data(1,aInt(jj,1):aInt(jj,2)), pos.data(2,aInt(jj,1):aInt(jj,2)), 'o', 'MarkerFaceColor', new_colors(jj,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
end
for jj=1:1:size(bInt,1)
    scatter(pos.data(1,bInt(jj,1):bInt(jj,2)), pos.data(2,bInt(jj,1):bInt(jj,2)), 'o', 'MarkerFaceColor', new_colors(jj,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
end
%% Automate collection of data for more individuals 
% Is 5 mc normally the case? how can I define the perimeter of my object ?
% Would this have to be done based on subject?

%Remove periods where the mouse is on top of the object (if any)