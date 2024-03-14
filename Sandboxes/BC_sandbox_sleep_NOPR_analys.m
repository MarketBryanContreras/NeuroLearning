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
video_flag=0;
save_output=1;
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

% get all sessions with 'BC'
inhib_dir = dir('*BC*');

%% Initializing outputs
    out=[];
%% Loop to load data from raw

for iS=1:length(inhib_dir)


    %% Colloecting subject info
    
    cd([inhib_dir(iS).folder filesep inhib_dir(iS).name])
    parts = pwd;
    parts= split(parts,filesep);
    parts=parts{end};
    parts= split(parts,'_');
    info.subject=parts{1};
    info.date=[ parts{2} '_' parts{3} '_' parts{4}];
    info.session=parts{5};
    
   
    %% Individual parameters
    if info.subject=="BC1807";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC3.ncs';
    elseif info.subject=="BC054";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC2.ncs';
    elseif info.subject=="BC053";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC2.ncs';
    elseif info.subject=="BC051";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC2.ncs';
    elseif info.subject=="BC014";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC2.ncs';
    elseif info.subject=="BC013";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC2.ncs';
    elseif info.subject=="BC011";
        emg_chan = 'CSC1.ncs';lfp_chan = 'CSC2.ncs';
    end
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

    if info.session=="D2"
        start_NOPR = evts.t{find(contains(evts.label, 'Starting Recording'))}(3);
        end_NOPR = evts.t{find(contains(evts.label, 'Stopping Recording'))}(3);
    end

    %Printing the duration of OF and sleeping
    fprintf('<strong>OF duration: %.2f mins</strong>\n', (end_OF - start_OF)/60);
    fprintf('<strong>Sleep duration: %.2f mins = %.2f hrs and %.2f min </strong>\n', (end_sleep - start_sleep)/60,floor(((end_sleep - start_sleep)/60)/60), ((end_sleep - start_sleep)/60) -(60*(floor(((end_sleep - start_sleep)/60)/60))));

    if info.session=="D2"
        fprintf('<strong>NOPR duration: %.2f mins</strong>\n', ((end_NOPR - start_NOPR)/60));
    end

    %Restrict the data to just the sleep phase.
    csc_s = restrict(csc, start_sleep, end_sleep);
    emg_s = restrict(emg, start_sleep, end_sleep);

    %Correcting times
    csc_s.tvec= csc_s.tvec-csc_s.tvec(1);
    emg_s.tvec=emg_s.tvec- emg_s.tvec(1);

    %% plot a bit of data for quality check end verify the awake states
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
    if info.session== "D1"
        if strcmpi(info.subject,"BC1807")
            wake_t = [0 3292 3657 5083 8332 10042 12118 13092];
        elseif strcmpi(info.subject,"BC053") %Specify the wake period according to subject and session
            wake_t = [0 3093 5015 5631 8824 9340 10605 10880 11531 11803 12450 13020 13618 13743 14188 14397];
        elseif strcmpi(info.subject,"BC011")
            wake_t = [0 1493 3455 4407 10405 13405 14654 14710];
        elseif strcmpi(info.subject,"BC013")
            wake_t = [0 4296 4735 5514 9823 10377 12838 13281];
        elseif strcmpi(info.subject,"BC014")
            wake_t = [0 4869 6952 9715 11557 11861 12753 13105 13301 13531];
        elseif strcmpi(info.subject,"BC051") % check this animal, apparenlty there was only 2 hrs of sleep 
            wake_t = [0 954 6238 8179];
        elseif strcmpi(info.subject,"BC054") 
            wake_t = [0 1423 2710 3859 4097 5355 6072 7153 9607 10098 11770 13206 13683 14402];
        end
    elseif info.session== "D2"
        if strcmpi(info.subject,"BC1807")
            wake_t = [0 4735 7179 9298 10043 11642 13919 14403];
        elseif strcmpi(info.subject,"BC013")
            wake_t = [0 688 1817 3130 4335 4783 6820 7511 10772 11917];
        elseif strcmpi(info.subject,"BC014")
            wake_t = [0 3730 ];
        elseif strcmpi(info.subject,"BC051") % check this animal, apparenlty there was only 2 hrs of sleep
            wake_t = [0 2453 2944 4320 4721 5119 6653 7569 10306 11233 13043 13639];
        elseif strcmpi(info.subject,"BC053") %Specify the wake period according to subject and session
            %wake_t = [0 3093 5015 5631 8824 9340 10605 10880 11531 11803 12450 13020 13618 13743 14188 14397];
        elseif strcmpi(info.subject,"BC054") % check this animal, apparenlty there was only 2 hrs of sleep
            wake_t = [0 2518 4597 5331 6521 8947 9211 7153 9607 10094 12110 12942];
        end
    end
    
    wake_idx = nearest_idx(wake_t, csc_s.tvec); %Converts time to sample idx
    wake_idx = reshape(wake_idx,2, length(wake_idx)/2)'; %reshape(columns, rows)
    %% score the sleep.

    [hypno, csc_out, emg_out] = dSub_Sleep_screener(plot_flag, csc_s, emg_s, wake_idx);  % can add in 'wake_idx' as the last input.
    
    %% Getting the percentage of sleep sates
[y,x]=histcounts(hypno.data,[0.5:1:3.5]);
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
    out.(info.subject).sleep_labels=hypno.labels;


    %% Tracking movement and setting up intervals when it explored the objects
    if info.session == "D2"
        % Load the DLC data
        distances=[];

        pos = MS_DLC2TSD_divide(cd, [], [4.5 4.5]); %Need to ask Eric how to solve for the size of the cage. For now we will assume its 4.5
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
        %% Plot the position of the mouse
        %%---To do--- Adapt this cell to the new structures
       
        
        minFrames=5; % The minimum number of frames where the mouse is in radious
        
        for iF=1:nfiles
            %Initialize figure
            figure (20+iF)
            %Actual position
            subplot(8,1,1:4);
            %plot the actual position
            time_frames=[1:length(pos.("File"+iF).data)];
            [aX,aY]=BC_Circle_plot(5.5,pos.("File"+iF).data(5,1),pos.("File"+iF).data(6,1), BC_color_genertor('Swamp_green'),plot_flag);
            hold on
            [bX,bY]=BC_Circle_plot(5.5,pos.("File"+iF).data(7,1),pos.("File"+iF).data(8,1), BC_color_genertor('Torment_blue'),plot_flag);
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


                scatter(pos.("File"+iF).data(5,1),pos.("File"+iF).data(6,1) , 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'none')
                scatter(pos.("File"+iF).data(7,1),pos.("File"+iF).data(8,1) , 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
                xlim([min(pos.("File"+iF).data(1,:)) max(pos.("File"+iF).data(1,:))]);
                ylim([min(pos.("File"+iF).data(2,:)) max(pos.("File"+iF).data(2,:))]);

                if video_flag==0;
                    scatter(pos.("File"+iF).data(1,:), pos.("File"+iF).data(2,:), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.01);
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
            %%% Rest of the plot
            %plot X position of the mouse
            subplot(8,1,5);
            plot(time_frames,pos.("File"+iF).data(1,:));
            xlim([0,time_frames(end)])
            %xlabel('Time(s)');
            ylabel('X pos (cm)');
            %plot Y position of the mouse
            subplot(8,1,6);
            plot(time_frames,pos.("File"+iF).data(2,:), 'Color',BC_color_genertor('Oxford_blue'));
            xlim([0,time_frames(end)])
            %xlabel('Time(s)');
            ylabel('Y pos (cm)');
            %Plot the distances to object 1
            subplot(8,1,7);
            plot(time_frames, distances.("File"+iF)(1,:), 'm');
            yline(5, 'Color', 'r')
            xlim([0,time_frames(end)])
            %xlabel('Time(s)');
            ylabel('Distance to object A');
            %Plot the distnaces to object 2
            subplot(8,1,8);
            plot(time_frames, distances.("File"+iF)(2,:),'r');
            yline(5, 'Color', 'r')
            xlim([0,time_frames(end)])
            xlabel('Time (s)');
            ylabel('Distance to object B');

            hold off

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
end

% Formating for saving output
if save_output
    cd(inter_dir)
    save(["out-" + date + ".mat"],'out')
    
end
%% Plot to compare time and # of interactions
%% Collect and make table
tbl = table(); 
stats = []; 
subject = []; 
trial_n = [];
opto = [];
cohort = []; 
SG_modidx = []; 
FG_modidx = []; 
t_bp=[];
fg_bp=[];
sg_bp=[];

all_files=fieldnames(out);
control_list= {'BC011' 'BC013'  'BC014'};
archT_list= {'BC051' 'BC054' 'BC1807'  'BC053'};
controlidx=[];
archTidx=[];

for idx = 1:numel(all_files)
    current_field = all_files{idx};
    if any(strcmp(current_field, control_list))
        controlidx = [controlidx, idx];
    end
    if any(strcmp(current_field, archT_list))
        archTidx = [archTidx, idx];
    end
end


for iSub = 1:length(all_files)
    sess_list = fieldnames(out.archT_list{iSub})); 
 
    n_inhib = length(out_Archt_07_nov_23.(archt_list{iSub}).D4.t_bp_inhib); 
    n_noinhib = length(out_Archt_07_nov_23.(archt_list{iSub}).D4.t_bp_noinhib);
    
    
     subject= [subject repmat(iSub,1, n_inhib + n_noinhib)]; 
     
     cohort = [cohort repmat(1,1, n_inhib + n_noinhib)];
     
     trial_n = [trial_n, [1:n_inhib, 1:n_noinhib]]; 
     
     opto = [opto, logical([repmat(1,1, n_inhib), repmat(0,1, n_noinhib)])]; 
     
     SG_modidx = [SG_modidx [out_Archt_07_nov_23.(archt_list{iSub}).D4.modidx_SG_inhib, out_Archt_07_nov_23.(archt_list{iSub}).D4.modidx_SG_noinhib]];
     
     FG_modidx = [FG_modidx [out_Archt_07_nov_23.(archt_list{iSub}).D4.modidx_FG_inhib, out_Archt_07_nov_23.(archt_list{iSub}).D4.modidx_FG_noinhib]];

     t_bp = [t_bp [out_Archt_07_nov_23.(archt_list{iSub}).D4.t_bp_inhib, out_Archt_07_nov_23.(archt_list{iSub}).D4.t_bp_noinhib]];
     
     sg_bp = [sg_bp [out_Archt_07_nov_23.(archt_list{iSub}).D4.sg_bp_inhib, out_Archt_07_nov_23.(archt_list{iSub}).D4.sg_bp_noinhib]];
     
     fg_bp = [fg_bp [out_Archt_07_nov_23.(archt_list{iSub}).D4.fg_bp_inhib, out_Archt_07_nov_23.(archt_list{iSub}).D4.fg_bp_noinhib]];

%      z_SGInhb_modidx=[z_SGInhb_modidx [out_Archt_07_nov_23.(archt_list{iSub}).D4.z_SGInhb_modidx]];
%      
%      z_FGInhb_modidx=[z_FGInhb_modidx [out_Archt_07_nov_23.(archt_list{iSub}).D4.z_FGInhb_modidx]];
%      
%      z_FGNoInhb_modidx=[z_FGNoInhb_modidx [out_Archt_07_nov_23.(archt_list{iSub}).D4.z_FGNoInhb_modidx]];
     
     if iSub==length(archt_list)
     archt_tbl= table(subject', cohort', opto', trial_n',SG_modidx',FG_modidx',t_bp',sg_bp',fg_bp','VariableNames', {'Subject', 'Cohort', 'Opto', 'Trial', 'SG_modidx', 'FG_modidx','Theta_bp','SG_bp','FG_bp'});
     end
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