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
%% Load data
cd(inter_dir);
load ("out-13-Jun-2024.mat"); %Selecte the data to be loaded
% out=rmfield(out,"BC051");
%% Collect and make table for D2 interactions
%Initialize columns
subject = []; 
cohort = []; 
pre_post = []; 
no_obj_a_int = []; 
no_obj_b_int=[];
obj_a_time= [];
obj_b_time= [];
stats=[];

all_files=fieldnames(out);%All the nems of the files or subjects
control_list= {'BC011' 'BC013'  'BC014'};
archT_list= {'BC051' 'BC054' 'BC1807'  'BC053'};
controlidx=[];
archTidx=[];
%Creates a list of the idx of either the control or the experimental 
for idx = 1:numel(all_files)
    current_field = all_files{idx};
    if any(strcmp(current_field, control_list))
        controlidx = [controlidx, idx];
    end
    if any(strcmp(current_field, archT_list))
        archTidx = [archTidx, idx];
    end
end
%Run until here to initialize any part of the analysis
%% Interaction analysis
for iSub = 1:length(all_files)
    current_subject = all_files{iSub};
    
    subject= [subject repmat(iSub,1,2)];
    
    if any(strcmp(current_subject, control_list))
        cohort = [cohort, repmat(1,1,2)];
    end
    if any(strcmp(current_subject, archT_list))
        cohort = [cohort, repmat(2,1,2)];
    end

    pre_post = [pre_post 1 2];
    no_obj_a_int= [no_obj_a_int out.(all_files{iSub}).D2.Behavior1.nObjAInteractions];
    no_obj_a_int= [no_obj_a_int out.(all_files{iSub}).D2.Behavior2.nObjAInteractions];

    no_obj_b_int= [no_obj_b_int out.(all_files{iSub}).D2.Behavior1.nObjBInteractions];
    no_obj_b_int= [no_obj_b_int out.(all_files{iSub}).D2.Behavior2.nObjBInteractions];

    obj_a_time= [obj_a_time out.(all_files{iSub}).D2.Behavior1.ObjATime];
    obj_a_time= [obj_a_time out.(all_files{iSub}).D2.Behavior2.ObjATime];

    obj_b_time= [obj_b_time out.(all_files{iSub}).D2.Behavior1.ObjBTime];
    obj_b_time= [obj_b_time out.(all_files{iSub}).D2.Behavior2.ObjBTime];

end
%% collect all in one table
D2Int= table(subject', cohort', pre_post', no_obj_a_int', no_obj_b_int', obj_a_time', obj_b_time','VariableNames', {'Subject', 'Cohort', 'Pre_post', 'No_obj_a_int', 'No_obj_b_int', 'ObjATime','ObjBTime'});
%% Calculating the percent of total investigation time  and discrimination index
%Investigation time percentage
post= D2Int(D2Int.Pre_post==2,:);
pre= D2Int(D2Int.Pre_post==1,:);
novel_invg_time= post.ObjATime;
novel_int= post.No_obj_a_int;
familiar_int=post.No_obj_b_int;
familiar_invg_time=post.ObjBTime;
pcnt_invg_time=(novel_invg_time./(novel_invg_time+familiar_invg_time))*100;
pcnt_invg_int=(novel_int)./(novel_int+familiar_int)*100;

%Discrimination index
disc_indx_time= (novel_invg_time+familiar_invg_time)./(novel_invg_time-familiar_invg_time);
disc_indx_int=(novel_int+familiar_int)./(novel_int-familiar_int);
%No. ob interaction post/ no. of interaction pre 
pcnt_objA_int= (post.No_obj_a_int./pre.No_obj_a_int)*100;
pcnt_objB_int= (post.No_obj_b_int./pre.No_obj_b_int)*100;

cohort=post.Cohort;
subject=post.Subject;
metrics= table(subject, cohort,pcnt_invg_time,pcnt_invg_int,disc_indx_time, disc_indx_int, pcnt_objA_int, pcnt_objB_int, 'VariableNames', {'subject' 'cohort' 'pcnt_inv_time' 'pcnt_inv_int' 'disc_indx_time' 'disc_indx_int' 'pcnt_objA_int_post' 'pcnt_objB_int_post'} );
metrics_ctrl=metrics(metrics.cohort==1,:);
metrics_archT=metrics(metrics.cohort==2,:);
%% Some stats in the discrimination index
statsMetrics.discIdx_time= fitlme(metrics, 'disc_indx_time ~ 1+cohort + (1|subject)');
statsMetrics.discIdx_int= fitlme(metrics, 'disc_indx_int ~ 1+cohort + (1|subject)');
statsMetrics.pcntInv_time= fitlme(metrics, 'pcnt_inv_time ~ 1+cohort + (1|subject)');
statsMetrics.pcntInv_int= fitlme(metrics, 'pcnt_inv_int ~ 1+cohort + (1|subject)');
statsMetrics.pcntObjApost= fitlme(metrics, 'pcnt_objA_int_post ~ 1+cohort + (1|subject)');
statsMetrics.pcntObjBpost= fitlme(metrics, 'pcnt_objB_int_post ~ 1+cohort + (1|subject)');

%% Intrecations Plotting
figure (1)
%Comparison of time with object A
sp1=subplot(2,8,5:6)

h = boxchart(D2Int.Pre_post, D2Int.ObjATime, 'GroupByColor', D2Int.Cohort);
ylabel('Time (s)'); xlabel(''); title('OBJECT A (moved)');

% Set box colors
h(2).BoxFaceColor = BC_color_genertor('ArchT_green');
h(1).BoxFaceColor = BC_color_genertor('oxford_blue');

% % Set x-tick positions and labels
xticks([1 2]); % Set x-tick positions for "Pre" and "Post"
xticklabels({'Pre', 'Post'}); % Set x-tick labels

hold on
xline(1.5, 'k--'); % Line to separate "Pre" and "Post"
%legend({'Control', 'ArchT', ''});
set(gca,'fontsize', 16)
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
set(gca,'box','off');

%Comparison of time with object B
sp2=subplot(2,8,7:8)

h = boxchart(D2Int.Pre_post,D2Int.ObjBTime,'GroupByColor',D2Int.Cohort)
%ylabel('Time (s)'); 
xlabel(''); title('OBJECT B');


% Set box colors
h(2).BoxFaceColor = BC_color_genertor('ArchT_green');
h(1).BoxFaceColor = BC_color_genertor('oxford_blue');

% % Set x-tick positions and labels
xticks([1 2]); % Set x-tick positions for "Pre" and "Post"
xticklabels({'Pre', 'Post'}); % Set x-tick labels
yticklabels({''});

hold on
xline(1.5, 'k--'); % Line to separate "Pre" and "Post"
legend({'Control', 'ArchT', ''});
set(gca,'fontsize', 16)
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
set(gca,'box','off');


%Comparison of no. of interactions with object A
sp3=subplot(2,8,13:14)

h = boxchart(D2Int.Pre_post, D2Int.No_obj_a_int, 'GroupByColor', D2Int.Cohort);
ylabel('No. of interactions '); xlabel('');

% Set box colors
h(2).BoxFaceColor = BC_color_genertor('ArchT_green');
h(1).BoxFaceColor = BC_color_genertor('oxford_blue');

% % Set x-tick positions and labels
xticks([1 2]); % Set x-tick positions for "Pre" and "Post"
xticklabels({'Pre', 'Post'}); % Set x-tick labels



hold on
xline(1.5, 'k--'); % Line to separate "Pre" and "Post"
%legend({'Control', 'ArchT', ''});
set(gca,'fontsize', 16)
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
set(gca,'box','off');

%Comparison of no. of interactions with object B
sp4=subplot(2,8,15:16)

h = boxchart(D2Int.Pre_post, D2Int.No_obj_b_int, 'GroupByColor', D2Int.Cohort);
%ylabel('No. of interactions');
xlabel(''); 

% Set box colors
h(2).BoxFaceColor = BC_color_genertor('ArchT_green');
h(1).BoxFaceColor = BC_color_genertor('oxford_blue');

% % Set x-tick positions and labels
xticks([1 2]); % Set x-tick positions for "Pre" and "Post"
xticklabels({'Pre', 'Post'}); % Set x-tick labels
yticklabels({''});
hold on
xline(1.5, 'k--'); % Line to separate "Pre" and "Post"
%legend({'Control', 'ArchT', ''});
set(gca,'fontsize', 16)
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
set(gca,'box','off');

linkaxes([sp1 sp2],'y')
linkaxes([sp3 sp4],'y')
%% Plotting the metrics

%Pcnt of inv time
ff=figure(900);
clf;
subplot(1,6,1)
boxplot(metrics.pcnt_inv_time, metrics.cohort,'Labels', {'Control', 'ArchT'})
xtickangle(45);
hold on

title('% of Inv. Time');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

% ctrl_val=(metrics.pcnt_inv_time(metrics.cohort==1));
% exp_val=(metrics.pcnt_inv_time(metrics.cohort==2));
% [h,p]=ttest2(ctrl_val,exp_val)

subplot(1,6,2)
boxplot(metrics.pcnt_inv_int,  metrics.cohort,'Labels', {'Control', 'ArchT'})
xtickangle(45);
hold on

title('% of Inv. Interactions');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

ctrl_val=(metrics.pcnt_inv_int(metrics.cohort==1));
exp_val=(metrics.pcnt_inv_int(metrics.cohort==2));
[h,p]=ttest2(ctrl_val,exp_val)

subplot(1,6,3)
boxplot(metrics.disc_indx_time,  metrics.cohort,'Labels', {'Control', 'ArchT'})
xtickangle(45);
hold on
title('Disc. idx time');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

% ctrl_val=(mmetrics.disc_indx_time(metrics.cohort==1))
% exp_val=(metrics.disc_indx_time(metrics.cohort==2))
%[h,p]=ttest2(ctrl_val,exp_val)

subplot(1,6,4)
boxplot(metrics.disc_indx_int,  metrics.cohort,'Labels', {'Control', 'ArchT'})
xtickangle(45);
hold on
title('Disc. idx interactions');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

% ctrl_val=(metrics.disc_indx_int(metrics.cohort==1))
% exp_val=(metrics.disc_indx_int(metrics.cohort==2))
% [h,p]=ttest2(ctrl_val,exp_val)


subplot(1,6,5)
boxplot(metrics.pcnt_objA_int_post,  metrics.cohort,'Labels', {'Control', 'ArchT'})
xtickangle(45);
hold on
title('% of interactions obj. A post');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

% ctrl_val=(metrics.pcnt_objA_int_post(metrics.cohort==1))
% exp_val=(metrics.pcnt_objA_int_post(metrics.cohort==2))
% [h,p]=ttest2(ctrl_val,exp_val)

subplot(1,6,6)
boxplot(metrics.pcnt_objB_int_post,  metrics.cohort,'Labels', {'Control', 'ArchT'})
xtickangle(45);
hold on
title('% of interactions obj. B post');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;          

% ctrl_val=(metrics.pcnt_objB_int_post(metrics.cohort==1))
% exp_val=(metrics.pcnt_objB_int_post(metrics.cohort==2))
% [h,p]=ttest2(ctrl_val,exp_val)

% Adjust figure properties
 %sgtitle('\bf Metrics');           % Add a general title above the entire figure
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 % Resize the figure (optional)
 fig.Position = [200, 200, 1000, 700];  % [x, y, width, height]
 %saveas(gcf, 'FG_NormModIdx.png');%%%% fill in nice plotting %%%%%%%

%% Stats
%Cheking if the archT pre-post has a significance on object A
archTd2Int = D2Int(D2Int.Cohort==2, :);
archtPreObjAInt=table2array(archTd2Int(archTd2Int.Pre_post==1,"ObjATime"));
archtPostObjAInt=table2array(archTd2Int(archTd2Int.Pre_post==2,"ObjATime"));
[h,p]=ttest(archtPreObjAInt,archtPostObjAInt);
stats.objectA_archT = fitlme(archTd2Int,'ObjATime ~ 1+ Pre_post +(1|Subject)');
%Cheking if the control pre-post has a significance on object A
controld2Int = D2Int(D2Int.Cohort==1, :);
stats.objectA_control = fitlme(controld2Int,'ObjATime ~ 1+ Pre_post +(1|Subject)');

% stats.t_bp = fitlme(tbl,'Theta_bp ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
% stats.sg_bp = fitlme(tbl,'SG_bp ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
% stats.fg_bp = fitlme(tbl,'FG_bp ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');

%% Collect data for the time on sleep, percentages,modulation index and bandpower
%Initialize columns
 
cohort = [];
day=[];
time_awk=[];
time_sws=[];
time_rem=[];
pcnt_awk=[];
pcnt_sws=[];
pcnt_rem=[];
FG_awk_modIdx=[];
SG_awk_modIdx=[];
FG_sws_modIdx=[];
SG_sws_modIdx=[];
FG_rem_modIdx=[];
SG_rem_modIdx=[];
Tta_awk_pwr=[];
Tta_sws_pwr=[];
Tta_rem_pwr=[];
SG_awk_pwr=[];
SG_sws_pwr=[];
SG_rem_pwr=[];
FG_awk_pwr=[];
FG_sws_pwr=[];
FG_rem_pwr=[];


for iSub = 1:length(all_files)
    current_subject = all_files{iSub};
    
    if any(strcmp(current_subject, control_list))
        cohort = [cohort, repmat(1,1,2)];
    end
    if any(strcmp(current_subject, archT_list))
        cohort = [cohort, repmat(2,1,2)];
    end

    day = [day 1 2];
    time_awk= [time_awk out.(all_files{iSub}).D1.sleep.times_sec(1)];
    time_awk= [time_awk out.(all_files{iSub}).D2.sleep.times_sec(1)];

    time_sws= [time_sws out.(all_files{iSub}).D1.sleep.times_sec(2)];
    time_sws= [time_sws out.(all_files{iSub}).D2.sleep.times_sec(2)];

    time_rem= [time_rem out.(all_files{iSub}).D1.sleep.times_sec(3)];
    time_rem= [time_rem out.(all_files{iSub}).D2.sleep.times_sec(3)];

    pcnt_awk= [pcnt_awk out.(all_files{iSub}).D1.sleep.percetages(1)];
    pcnt_awk= [pcnt_awk out.(all_files{iSub}).D2.sleep.percetages(1)];

    pcnt_sws= [pcnt_sws out.(all_files{iSub}).D1.sleep.percetages(2)];
    pcnt_sws= [pcnt_sws out.(all_files{iSub}).D2.sleep.percetages(2)];

    pcnt_rem= [pcnt_rem out.(all_files{iSub}).D1.sleep.percetages(3)];
    pcnt_rem= [pcnt_rem out.(all_files{iSub}).D2.sleep.percetages(3)];

    SG_awk_modIdx= [SG_awk_modIdx out.(all_files{iSub}).D1.sleep.mod_idx(1,1)];
    SG_awk_modIdx= [SG_awk_modIdx out.(all_files{iSub}).D2.sleep.mod_idx(1,1)];

    FG_awk_modIdx= [FG_awk_modIdx out.(all_files{iSub}).D1.sleep.mod_idx(1,2)];
    FG_awk_modIdx= [FG_awk_modIdx out.(all_files{iSub}).D2.sleep.mod_idx(1,2)];

    SG_sws_modIdx= [SG_sws_modIdx out.(all_files{iSub}).D1.sleep.mod_idx(2,1)];
    SG_sws_modIdx= [SG_sws_modIdx out.(all_files{iSub}).D2.sleep.mod_idx(2,1)];

    FG_sws_modIdx= [FG_sws_modIdx out.(all_files{iSub}).D1.sleep.mod_idx(2,2)];
    FG_sws_modIdx= [FG_sws_modIdx out.(all_files{iSub}).D2.sleep.mod_idx(2,2)];

    SG_rem_modIdx= [SG_rem_modIdx out.(all_files{iSub}).D1.sleep.mod_idx(3,1)];
    SG_rem_modIdx= [SG_rem_modIdx out.(all_files{iSub}).D2.sleep.mod_idx(3,1)];

    FG_rem_modIdx= [FG_rem_modIdx out.(all_files{iSub}).D1.sleep.mod_idx(3,2)];
    FG_rem_modIdx= [FG_rem_modIdx out.(all_files{iSub}).D2.sleep.mod_idx(3,2)];

    Tta_awk_pwr= [Tta_awk_pwr out.(all_files{iSub}).D1.sleep.powers.awk.Tta];
    Tta_awk_pwr= [Tta_awk_pwr out.(all_files{iSub}).D2.sleep.powers.awk.Tta];

    Tta_sws_pwr= [Tta_sws_pwr out.(all_files{iSub}).D1.sleep.powers.sws.Tta];
    Tta_sws_pwr= [Tta_sws_pwr out.(all_files{iSub}).D2.sleep.powers.sws.Tta];

    Tta_rem_pwr= [Tta_rem_pwr out.(all_files{iSub}).D1.sleep.powers.rem.Tta];
    Tta_rem_pwr= [Tta_rem_pwr out.(all_files{iSub}).D2.sleep.powers.rem.Tta];

    SG_awk_pwr= [SG_awk_pwr out.(all_files{iSub}).D1.sleep.powers.awk.SG];
    SG_awk_pwr= [SG_awk_pwr out.(all_files{iSub}).D2.sleep.powers.awk.SG];

    SG_sws_pwr= [SG_sws_pwr out.(all_files{iSub}).D1.sleep.powers.sws.SG];
    SG_sws_pwr= [SG_sws_pwr out.(all_files{iSub}).D2.sleep.powers.sws.SG];

    SG_rem_pwr= [SG_rem_pwr out.(all_files{iSub}).D1.sleep.powers.rem.SG];
    SG_rem_pwr= [SG_rem_pwr out.(all_files{iSub}).D2.sleep.powers.rem.SG];

    FG_awk_pwr= [FG_awk_pwr out.(all_files{iSub}).D1.sleep.powers.awk.FG];
    FG_awk_pwr= [FG_awk_pwr out.(all_files{iSub}).D2.sleep.powers.awk.FG];

    FG_sws_pwr= [FG_sws_pwr out.(all_files{iSub}).D1.sleep.powers.sws.FG];
    FG_sws_pwr= [FG_sws_pwr out.(all_files{iSub}).D2.sleep.powers.sws.FG];

    FG_rem_pwr= [FG_rem_pwr out.(all_files{iSub}).D1.sleep.powers.rem.FG];
    FG_rem_pwr= [FG_rem_pwr out.(all_files{iSub}).D2.sleep.powers.rem.FG];

end

%% Putting sleep data it in a table
sleep= table(subject', cohort', day', time_awk', time_sws', time_rem', pcnt_awk', pcnt_sws', pcnt_rem', SG_awk_modIdx',FG_awk_modIdx',SG_sws_modIdx', FG_sws_modIdx',SG_rem_modIdx', FG_rem_modIdx', Tta_awk_pwr', Tta_sws_pwr',Tta_rem_pwr',SG_awk_pwr',SG_sws_pwr',SG_rem_pwr', FG_awk_pwr', FG_sws_pwr',FG_rem_pwr','VariableNames', {'Subject', 'Cohort', 'Day', 'time_awk', 'time_sws', 'time_rem', 'pcnt_awk','pcnt_sws', 'pcnt_rem','SG_awk_modIdx','FG_awk_modIdx','SG_sws_modIdx', 'FG_sws_modIdx','SG_rem_modIdx', 'FG_rem_modIdx', 'Tta_awk_pwr', 'Tta_sws_pwr','Tta_rem_pwr','SG_awk_pwr','SG_sws_pwr','SG_rem_pwr', 'FG_awk_pwr', 'FG_sws_pwr','FG_rem_pwr'});
%% Stats for the modulation index
%Stats for the modulation index 
% stats=[];
% stats.SG_awk = fitlme(sleep, ...
%     'SG_awk_modIdx ~ 1 + Day*Cohort + (1|Subject)', ...
%     'FitMethod', 'REML');
% disp(stats.SG_awk);
% 
% stats.FG_awk = fitlme(sleep, ...
%     'FG_awk_modIdx ~ 1 + Day*Cohort + (1|Subject)', ...
%     'FitMethod', 'REML');
% disp(stats.FG_awk);
% 
% stats.SG_sws = fitlme(sleep, ...
%     'SG_sws_modIdx ~ 1 + Day*Cohort + (1|Subject)', ...
%     'FitMethod', 'REML');
% disp(stats.SG_sws);
% 
% stats.FG_sws = fitlme(sleep, ...
%     'FG_sws_modIdx ~ 1 + Day*Cohort + (1|Subject)', ...
%     'FitMethod', 'REML');
% disp(stats.FG_sws);
% 
% stats.SG_rem = fitlme(sleep, ...
%     'SG_rem_modIdx ~ 1 + Day*Cohort + (1|Subject)', ...
%     'FitMethod', 'REML');
% disp(stats.SG_rem);
% 
% stats.FG_rem = fitlme(sleep, ...
%     'FG_rem_modIdx ~ 1 + Day*Cohort + (1|Subject)', ...
%     'FitMethod', 'REML');
% disp(stats.FG_rem);
%% Stats for the times
%D1 intercohort
D1Slp=sleep(sleep.Day==1,:);
StatsSlpTime.D1AwkInterCohort= fitlme(D1Slp,'time_awk ~1+Cohort + (1|Subject)');
StatsSlpTime.D1SwsInterCohort= fitlme(D1Slp,'time_sws ~1+Cohort + (1|Subject)');
StatsSlpTime.D1RemInterCohort= fitlme(D1Slp,'time_rem ~1+Cohort + (1|Subject)');
%D2 intercohort
D2Slp=sleep(sleep.Day==2,:);
StatsSlpTime.D2AwkInterCohort= fitlme(D2Slp,'time_awk ~1+Cohort + (1|Subject)');
StatsSlpTime.D2SwsInterCohort= fitlme(D2Slp,'time_sws ~1+Cohort + (1|Subject)');
StatsSlpTime.D2RemInterCohort= fitlme(D2Slp,'time_rem ~1+Cohort + (1|Subject)');

%Control intracohort
C1Slp=sleep(sleep.Cohort==1,:);
StatsSlpTime.C1AwkIntraCohort= fitlme(C1Slp,'time_awk ~1+Day + (1|Subject)');
StatsSlpTime.C1SwsIntraCohort= fitlme(C1Slp,'time_sws ~1+Day + (1|Subject)');
StatsSlpTime.C1RemIntraCohort= fitlme(C1Slp,'time_rem ~1+Day + (1|Subject)');
%Experimental intracohort
C2Slp=sleep(sleep.Cohort==2,:);
StatsSlpTime.C2AwkIntraCohort= fitlme(C2Slp,'time_awk ~1+Day + (1|Subject)');
StatsSlpTime.C2SwsIntraCohort= fitlme(C2Slp,'time_sws ~1+Day + (1|Subject)');
StatsSlpTime.C2RemIntraCohort= fitlme(C2Slp,'time_rem ~1+Day + (1|Subject)');


%% Sleep data plotting dashboard
% Compariosn of time pcnt spent in each sleep phase
day1_cohort_1=sleep(sleep.Day==1 & sleep.Cohort==1, :);
day2_cohort_1=sleep(sleep.Day==2 & sleep.Cohort==1, :);
day1_cohort_2=sleep(sleep.Day==1 & sleep.Cohort==2, :);
day2_cohort_2=sleep(sleep.Day==2 & sleep.Cohort==2, :);

%Extracting sleep times by cohort
d1_c1_time_awk=day1_cohort_1.time_awk;
d1_c2_time_awk=day1_cohort_2.time_awk;
d2_c1_time_awk=day2_cohort_1.time_awk;
d2_c2_time_awk=day2_cohort_2.time_awk;

d1_c1_time_sws=day1_cohort_1.time_sws;
d1_c2_time_sws=day1_cohort_2.time_sws;
d2_c1_time_sws=day2_cohort_1.time_sws;
d2_c2_time_sws=day2_cohort_2.time_sws;

d1_c1_time_rem=day1_cohort_1.time_rem;
d1_c2_time_rem=day1_cohort_2.time_rem;
d2_c1_time_rem=day2_cohort_1.time_rem;
d2_c2_time_rem=day2_cohort_2.time_rem;

%Extracting sleep percentages by cohort

d1_c1_pcnt=[day1_cohort_1(:,4:6)];
d1_c2_pcnt=[day1_cohort_2(:,4:6)];
d2_c1_pcnt=[day2_cohort_1(:,4:6)];
d2_c2_pcnt=[day2_cohort_2(:,4:6)];

d1_c1_sum=table2array(sum(d1_c1_pcnt));
d1_c1_pcnt=d1_c1_sum/sum(d1_c1_sum);

d1_c2_sum=table2array(sum(d1_c2_pcnt));
d1_c2_pcnt=d1_c2_sum/sum(d1_c2_sum);

d2_c1_sum=table2array(sum(d2_c1_pcnt));
d2_c1_pcnt=d2_c1_sum/sum(d2_c1_sum);

d2_c2_sum=table2array(sum(d2_c2_pcnt));
d2_c2_pcnt=d2_c2_sum/sum(d2_c2_sum);



        

meanD1C1awk= mean(d1_c1_time_awk);
meanD1C1sws= mean(d1_c1_time_sws);
meanD1C1rem= mean(d1_c1_time_rem);

meanD2C1awk= mean(d2_c1_time_awk);
meanD2C1sws= mean(d2_c1_time_sws);
meanD2C1rem= mean(d2_c1_time_rem);

meanD1C2awk= mean(d1_c2_time_awk);
meanD1C2sws= mean(d1_c2_time_sws);
meanD1C2rem= mean(d1_c2_time_rem);

meanD2C2awk= mean(d2_c2_time_awk);
meanD2C2sws= mean(d2_c2_time_sws);
meanD2C2rem= mean(d2_c2_time_rem);


steD1C1awk= std(d1_c1_time_awk)/sqrt(size(d1_c1_time_awk,1));
steD1C1sws= std(d1_c1_time_sws)/sqrt(size(d1_c1_time_sws,1));
steD1C1rem= std(d1_c1_time_rem);sqrt(size(d1_c1_time_rem,1));

steD2C1awk= std(d2_c1_time_awk)/sqrt(size(d2_c1_time_awk,1));
steD2C1sws= std(d2_c1_time_sws)/sqrt(size(d2_c1_time_sws,1));
steD2C1rem= std(d2_c1_time_rem)/sqrt(size(d2_c1_time_rem,1));

steD1C2awk= std(d1_c2_time_awk)/sqrt(size(d1_c2_time_awk,1));
steD1C2sws= std(d1_c2_time_sws)/sqrt(size(d1_c2_time_sws,1));
steD1C2rem= std(d1_c2_time_rem)/sqrt(size(d1_c2_time_rem,1));

steD2C2awk= std(d2_c2_time_awk)/sqrt(size(d2_c2_time_awk,1));
steD2C2sws= std(d2_c2_time_sws)/sqrt(size(d2_c2_time_sws,1));
steD2C2rem= std(d2_c2_time_rem)/sqrt(size(d2_c2_time_rem,1));



xx= [1,2,3,4];
yy=[meanD1C1awk,meanD2C1awk,meanD1C2awk,meanD2C2awk;
    meanD1C1sws,meanD2C1sws,meanD1C2sws,meanD2C2sws;
    meanD1C1rem,meanD2C1rem,meanD1C2rem,meanD2C2rem];
ste=[steD1C1awk,steD2C1awk,steD1C2awk,steD2C2awk;
    steD1C1sws,steD2C1sws,steD1C2sws,steD2C2sws;
    steD1C1rem,steD2C1rem,steD1C2rem,steD2C2rem];
maxy=max(max(yy));
maxy=maxy+maxy*0.4;

subplot(2,8,[9:11])

num=4;
c=1:num;
control_x=[0.5 2.5 2.5 0.5];
control_y=[0 0 maxy maxy];
fill(control_x, control_y, BC_color_genertor('Oxford_blue'), 'FaceAlpha', 0.25 , 'LineStyle',"none")
hold on

archt_x=[2.5 5.0 5.0 2.5];
archt_y=[0 0 maxy maxy];
fill(archt_x, archt_y, BC_color_genertor('Archt_green'), 'FaceAlpha', 0.25 , 'LineStyle',"none")

for ii = 1:num
    bar(c(ii)-0.25,yy(1,ii),0.2,'EdgeColor','none');
    bar(c(ii),yy(2,ii),0.2,'EdgeColor','none');
    bar(c(ii)+0.25,yy(3,ii),0.2,'EdgeColor','none');
end
% Hold the plot for adding error bars


% Plot error bars
errH1 = errorbar(c-0.25,yy(1,:),ste(1,:),'.');
errH2 = errorbar(c,yy(2,:),ste(2,:),'.');
errH3 = errorbar(c+0.25,yy(3,:),ste(3,:),'.');
errH1.LineWidth = 1.25;
errH2.LineWidth = 1.25;
errH3.LineWidth = 1.25;

errH1.Color = [0 0 0];
errH2.Color = [0 0 0];
errH3.Color = [0 0 0];


ylim([0 maxy])
xlim([0.5 4.5])

xline(2.5,'k--');
legend({'Control','ArchT','Awake', 'SWS', 'REM', ''});
title('Time spent on corresponding sleep phase');
ylabel('Time (s)');
sleep_colors= [...
    0.3467    0.5360    0.6907;
    0.9153    0.2816    0.2878;
    0.4416    0.7490    0.4322];
colororder(sleep_colors);
set(gca,'fontsize', 16)
ax = gca;
ax.XTick = [1,2,3,4];
xticklabels({'D1', 'D2', 'D1', 'D2'}); % Set x-tick labels
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
set(gca,'box','off');

%Percentages D1 C1
subplot(4,8,1:2)
donutchart(d1_c1_pcnt, {'','',''});

%Percentages D1 C2
subplot(4,8,9:10)
donutchart(d1_c2_pcnt, {'','',''});

%Percentages D2 C1
subplot(4,8,3:4)
donutchart(d2_c1_pcnt, {'','',''});

%Percentages D2 C2
subplot(4,8,11:12)
donutchart(d2_c2_pcnt, {'','',''});

annotation('textbox', [0.1, 0.76, 0.1, 0.1], 'String', 'Control', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 22, 'FontWeight', 'bold');
annotation('textbox', [0.102, 0.54, 0.1, 0.1], 'String', 'ArhT', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 22, 'FontWeight', 'bold');

annotation('textbox', [0.195, 0.87 , 0.1, 0.1], 'String', 'Day 1', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 22, 'FontWeight', 'bold');
annotation('textbox', [0.405, 0.87 , 0.1, 0.1], 'String', 'Day 2', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 22, 'FontWeight', 'bold');



fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1];
%% [Not necesaary for plotting]Extracting SG Mod idx by cohort for the figure
 sleep_phases = {'Awake', 'SWS', 'REM'};
% cohorts = {'Control', 'ArchT'};
% days = {'D1', 'D2'};
% %Extracting SG Mod idx by cohort
% d1_c1_SG_awk_modIdx=day1_cohort_1.SG_awk_modIdx;
% d1_c2_SG_awk_modIdx=day1_cohort_2.SG_awk_modIdx;
% d1_c1_SG_sws_modIdx=day1_cohort_1.SG_sws_modIdx;
% d1_c2_SG_sws_modIdx=day1_cohort_2.SG_sws_modIdx;
% d1_c1_SG_rem_modIdx=day1_cohort_1.SG_rem_modIdx;
% d1_c2_SG_rem_modIdx=day1_cohort_2.SG_rem_modIdx;
% 
% d2_c1_SG_awk_modIdx=day2_cohort_1.SG_awk_modIdx;
% d2_c2_SG_awk_modIdx=day2_cohort_2.SG_awk_modIdx;
% d2_c1_SG_sws_modIdx=day2_cohort_1.SG_sws_modIdx;
% d2_c2_SG_sws_modIdx=day2_cohort_2.SG_sws_modIdx;
% d2_c1_SG_rem_modIdx=day2_cohort_1.SG_rem_modIdx;
% d2_c2_SG_rem_modIdx=day2_cohort_2.SG_rem_modIdx;
% 
% %Extracting FG Mod idx by cohort
% d1_c1_FG_awk_modIdx=day1_cohort_1.FG_awk_modIdx;
% d1_c2_FG_awk_modIdx=day1_cohort_2.FG_awk_modIdx;
% d1_c1_FG_sws_modIdx=day1_cohort_1.FG_sws_modIdx;
% d1_c2_FG_sws_modIdx=day1_cohort_2.FG_sws_modIdx;
% d1_c1_FG_rem_modIdx=day1_cohort_1.FG_rem_modIdx;
% d1_c2_FG_rem_modIdx=day1_cohort_2.FG_rem_modIdx;
% 
% d2_c1_FG_awk_modIdx=day2_cohort_1.FG_awk_modIdx;
% d2_c2_FG_awk_modIdx=day2_cohort_2.FG_awk_modIdx;
% d2_c1_FG_sws_modIdx=day2_cohort_1.FG_sws_modIdx;
% d2_c2_FG_sws_modIdx=day2_cohort_2.FG_sws_modIdx;
% d2_c1_FG_rem_modIdx=day2_cohort_1.FG_rem_modIdx;
% d2_c2_FG_rem_modIdx=day2_cohort_2.FG_rem_modIdx;

%% Plottting first figure for SG ModIDX
 sleep_phases = {'Awake', 'SWS', 'REM'};

% Combine data into one array for easier plotting
SGData = [sleep.SG_awk_modIdx sleep.SG_sws_modIdx sleep.SG_rem_modIdx];
FGData= [sleep.FG_awk_modIdx sleep.FG_sws_modIdx sleep.FG_rem_modIdx];
groupCohorts = [sleep.Cohort];
groupDays = [sleep.Day];
sg=figure(91);
clf;
for iSP=1:3
    %subplot(3,1,iSP)
    boxplot(SGData(:, iSP), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'}, 'LabelOrientation', 'inline');
    title([sleep_phases{1,iSP}]);
    ylabel('Mod Idx (AU)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',16);

end
ttl=annotation('textbox',[.35 0.9 0.1 0.1],'String','Slow Gamma PAC', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
sg.Position=[ 2000 0 400 1200]

% Figure of FG modIdx only REM
boxplot(FGData(:, 3), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'});
    title(' Fast Gamma during REM ');
    ylabel('Mod Idx (AU)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45)
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;

% Define colors
colors = [
     BC_color_genertor('Archt_green');
     BC_color_genertor('Oxford_blue');
     BC_color_genertor('Swamp_green');
    BC_color_genertor('Powder_blue'); ];

% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 
 % Figure of SG modIdx only REM
boxplot(SGData(:, 3), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'});
    title(' Fast Gamma during REM ');
    ylabel('Mod Idx (AU)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45)
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;

% Define colors
colors = [
     BC_color_genertor('Archt_green');
     BC_color_genertor('Oxford_blue');
     BC_color_genertor('Swamp_green');
    BC_color_genertor('Powder_blue'); ];

% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white



%% Plottting first figure for FG ModIDX

fg=figure(92);
clf;
for iSP=1:3
    subplot(3,1,iSP)
    boxplot(FGData(:, iSP), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'}, 'LabelOrientation', 'inline');
    title([sleep_phases{1,iSP}]);
    ylabel('Mod Idx (AU)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',16);

end
ttl=annotation('textbox',[.35 0.9 0.1 0.1],'String','Fast Gamma PAC', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
fg.Position=[ 2400 0 400 1200]

%% [Not necesaary for plotting]Extracting data for the power figure 
% %Extracting Tta pwr by cohort and day
% d1_c1_Tta_awk_pwr=day1_cohort_1.Tta_awk_pwr;
% d1_c2_Tta_awk_pwr=day1_cohort_2.Tta_awk_pwr;
% d1_c1_Tta_sws_pwr=day1_cohort_1.Tta_sws_pwr;
% d1_c2_Tta_sws_pwr=day1_cohort_2.Tta_sws_pwr;
% d1_c1_Tta_rem_pwr=day1_cohort_1.Tta_rem_pwr;
% d1_c2_Tta_rem_pwr=day1_cohort_2.Tta_rem_pwr;
% 
% d2_c1_Tta_awk_pwr=day2_cohort_1.Tta_awk_pwr;
% d2_c2_Tta_awk_pwr=day2_cohort_2.Tta_awk_pwr;
% d2_c1_Tta_sws_pwr=day2_cohort_1.Tta_sws_pwr;
% d2_c2_Tta_sws_pwr=day2_cohort_2.Tta_sws_pwr;
% d2_c1_Tta_rem_pwr=day2_cohort_1.Tta_rem_pwr;
% d2_c2_Tta_rem_pwr=day2_cohort_2.Tta_rem_pwr;
% 
% %Extracting SG pwr by cohort and day
% d1_c1_SG_awk_pwr=day1_cohort_1.SG_awk_pwr;
% d1_c2_SG_awk_pwr=day1_cohort_2.SG_awk_pwr;
% d1_c1_SG_sws_pwr=day1_cohort_1.SG_sws_pwr;
% d1_c2_SG_sws_pwr=day1_cohort_2.SG_sws_pwr;
% d1_c1_SG_rem_pwr=day1_cohort_1.SG_rem_pwr;
% d1_c2_SG_rem_pwr=day1_cohort_2.SG_rem_pwr;
% 
% d2_c1_SG_awk_pwr=day2_cohort_1.SG_awk_pwr;
% d2_c2_SG_awk_pwr=day2_cohort_2.SG_awk_pwr;
% d2_c1_SG_sws_pwr=day2_cohort_1.SG_sws_pwr;
% d2_c2_SG_sws_pwr=day2_cohort_2.SG_sws_pwr;
% d2_c1_SG_rem_pwr=day2_cohort_1.SG_rem_pwr;
% d2_c2_SG_rem_pwr=day2_cohort_2.SG_rem_pwr;
% 
% %Extracting FG pwr by cohort and day
% d1_c1_FG_awk_pwr=day1_cohort_1.SG_awk_pwr;
% d1_c2_FG_awk_pwr=day1_cohort_2.SG_awk_pwr;
% d1_c1_FG_sws_pwr=day1_cohort_1.SG_sws_pwr;
% d1_c2_FG_sws_pwr=day1_cohort_2.SG_sws_pwr;
% d1_c1_FG_rem_pwr=day1_cohort_1.SG_rem_pwr;
% d1_c2_FG_rem_pwr=day1_cohort_2.SG_rem_pwr;
% 
% d2_c1_FG_awk_pwr=day2_cohort_1.FG_awk_pwr;
% d2_c2_FG_awk_pwr=day2_cohort_2.FG_awk_pwr;
% d2_c1_FG_sws_pwr=day2_cohort_1.FG_sws_pwr;
% d2_c2_FG_sws_pwr=day2_cohort_2.FG_sws_pwr;
% d2_c1_FG_rem_pwr=day2_cohort_1.FG_rem_pwr;
% d2_c2_FG_rem_pwr=day2_cohort_2.FG_rem_pwr;
%% Plottting figure for Theta PWR

% Combine data into one array for easier plotting
Tta_pwr_Data = [sleep.Tta_awk_pwr sleep.Tta_sws_pwr sleep.Tta_rem_pwr];
SG_pwr_Data= [sleep.SG_awk_pwr sleep.SG_sws_pwr sleep.SG_rem_pwr];
FG_pwr_Data= [sleep.FG_awk_pwr sleep.FG_sws_pwr sleep.FG_rem_pwr];

groupCohorts = [sleep.Cohort];
groupDays = [sleep.Day];
sg=figure(101);
clf;
for iSP=1:3
    subplot(3,1,iSP)
    boxplot(Tta_pwr_Data(:, iSP), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'}, 'LabelOrientation', 'inline');
    title([sleep_phases{1,iSP}]);
    ylabel('Power (mW)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',16);set(gca, 'Box','off')
end



ttl=annotation('textbox',[.35 0.9 0.1 0.1],'String','Theta Band Power', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
sg.Position=[ 2000 0 400 1200]

% Figure of Theta pwr only REM
sg=figure(101);
boxplot(Tta_pwr_Data(:, 3), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'});
    title('Theta Power during REM');
    ylabel('Power (mW)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45)
% Customize line and whisker colors

set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;

% Define colors
colors = [
     BC_color_genertor('Archt_green');
     BC_color_genertor('Oxford_blue');
     BC_color_genertor('Swamp_green');
    BC_color_genertor('Powder_blue'); ];

% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 sg.Position=[ 2500 100 900 1000]
%% Plottting figure for SG PWR
sg=figure(102);
clf;
for iSP=1:3
    subplot(3,1,iSP)
    bp=boxplot(SG_pwr_Data(:, iSP), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'}, 'LabelOrientation','inline');
    title([sleep_phases{1,iSP}]);
    ylabel('Power (mW)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',16);set(gca, 'Box','off')
end

ttl=annotation('textbox',[.35 0.9 0.1 0.1],'String','SG Band Power', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
sg.Position=[ 2000 0 400 1200]

% Figure of SG pwr only REM
sg=figure(102);
boxplot(SG_pwr_Data(:, 3), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'});
    title('Slow Gamma Power during REM');
    ylabel('Power (mW)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45)
% Customize line and whisker colors

set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;

% Define colors
colors = [
     BC_color_genertor('Archt_green');
     BC_color_genertor('Oxford_blue');
     BC_color_genertor('Swamp_green');
    BC_color_genertor('Powder_blue'); ];

% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 sg.Position=[ 2500 100 900 1000]
%% Plottting figure for FG PWR
sg=figure(103);
clf;
for iSP=1:3
    subplot(3,1,iSP)
    bp=boxplot(FG_pwr_Data(:, iSP), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'}, 'LabelOrientation','inline');
    title([sleep_phases{1,iSP}]);
    ylabel('Power (mW)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',16);set(gca, 'Box','off')
end

ttl=annotation('textbox',[.35 0.9 0.1 0.1],'String','FG Band Power', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
sg.Position=[ 2000 0 400 1200]

% Figure of FG pwr only REM
sg=figure(102);
boxplot(FG_pwr_Data(:, 3), {groupCohorts, groupDays}, 'FactorSeparator', 1, ...
        'Labels', {'Control-D1', 'Control-D2', 'ArchT-D1', 'ArchT-D2'});
    title('Fast Gamma Power during REM');
    ylabel('Power (mW)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45)
% Customize line and whisker colors

set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;

% Define colors
colors = [
     BC_color_genertor('Archt_green');
     BC_color_genertor('Oxford_blue');
     BC_color_genertor('Swamp_green');
    BC_color_genertor('Powder_blue'); ];

% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 sg.Position=[ 2500 100 900 1000]

%% Stats for the power
stats_pwr=[];
%Theta
stats_pwr.Tta_awk = fitlme(sleep, ...
    'Tta_awk_pwr ~ 1 + Day*Cohort + (1|Subject)', ...
    'FitMethod', 'REML');
disp(stats_pwr.Tta_awk);

stats_pwr.Tta_sws = fitlme(sleep, ...
    'Tta_sws_pwr ~ 1 + Day*Cohort + (1|Subject)', ...
    'FitMethod', 'REML');
disp(stats_pwr.Tta_sws);

stats_pwr.Tta_rem = fitlme(sleep, ...
    'Tta_rem_pwr ~ 1 + Day*Cohort + (1|Subject)', ...
    'FitMethod', 'REML');
disp(stats_pwr.Tta_rem);
 
%SG
stats_pwr.SG_awk = fitlme(sleep, ...
    'SG_awk_pwr ~ 1 + Day*Cohort + (1|Subject)', ...
    'FitMethod', 'REML');
disp(stats_pwr.SG_awk);
 
stats_pwr.SG_sws = fitlme(sleep, ...
    'SG_sws_pwr ~ 1 + Day*Cohort + (1|Subject)', ...
    'FitMethod', 'REML');
disp(stats_pwr.SG_sws);

stats_pwr.SG_rem = fitlme(sleep, ...
    'SG_rem_pwr ~ 1 + Day*Cohort + (1|Subject)', ...
    'FitMethod', 'REML');
disp(stats_pwr.SG_rem);
%FG
stats_pwr.FG_awk = fitlme(sleep, ...
    'FG_awk_pwr ~ 1 + Day*Cohort + (1|Subject)', ...
    'FitMethod', 'REML');
disp(stats_pwr.FG_awk);

stats_pwr.FG_sws = fitlme(sleep, ...
    'FG_sws_pwr ~ 1 + Day*Cohort + (1|Subject)', ...
    'FitMethod', 'REML');
disp(stats_pwr.FG_sws);

stats_pwr.FG_rem = fitlme(sleep, ...
    'FG_rem_pwr ~ 1 + Day*Cohort + (1|Subject)', ...
    'FitMethod', 'REML');
disp(stats_pwr.FG_rem);
%% Calculating the Percentage of change in power for a more fair comaparison
%Pull the data from sleep into two different tables, 1 for D1 and another
%for D2

%Adjust according to the intervals with pwr
int1=1:2;
int2=(find(contains(fieldnames(sleep),'pwr')));
int3=(find(contains(fieldnames(sleep),'modIdx')));

colInt= [int1, int2, int3]; %Columns interval
Pwr_d1= sleep(sleep.Day==1,int2);
Pwr_d2= sleep(sleep.Day==2,int2);

Mod_d1=sleep(sleep.Day==1,:);
Mod_d2=sleep(sleep.Day==2,int3);

meta=sleep(sleep.Day==1,int1);
%Divide the coressponding columns of D2/D1 and multiply by 100 in pwr
pcntChange=(Pwr_d1./Pwr_d2).*100;
pcntChange= [meta, pcntChange];

%Divide the coressponding columns of D2/D1 and multiply by 100 in pwr
pcntModChange=(Mod_d1./Mod_d2).*100;
pcntModChange= [meta, pcntModChange];
% Define colors
colors = [
     BC_color_genertor('Archt_green');
     BC_color_genertor('Oxford_blue');
     ];

% Figure of percentage change only in rem
pcnt_rem_Data= [pcntChange.Tta_rem_pwr];
groupCohorts=[pcntChange.Cohort];
sg=figure(102);
boxplot(pcnt_rem_Data, {groupCohorts}, 'FactorSeparator', 1, ...
        'Labels', {'Control',  'ArchT'});
    title('Percentage change D1/D2 Theta Power');
    ylabel('Percentage (%)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45)
ylim([40 180])

% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;

% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 sg.Position=[ 2500 100 400 600]
 
% Figure of SG percentage change only in rem
pcnt_rem_Data= [pcntChange.SG_rem_pwr];
groupCohorts=[pcntChange.Cohort];
sg=figure(103);
boxplot(pcnt_rem_Data, {groupCohorts}, 'FactorSeparator', 1, ...
        'Labels', {'Control',  'ArchT'});
    title('Percentage change D1/D2 SG Power');
    ylabel('Percentage (%)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45)
ylim([40 180])

% Customize line and whisker colors

set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;
% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 sg.Position=[ 2500 100 400 600]

% Figure of FG percentage change only in rem
pcnt_rem_Data= [pcntChange.FG_rem_pwr];
groupCohorts=[pcntChange.Cohort];
sg=figure(104);
boxplot(pcnt_rem_Data, {groupCohorts}, 'FactorSeparator', 1, ...
        'Labels', {'Control',  'ArchT'});
    title('Percentage change D1/D2 FG Power');
    ylabel('Percentage (%)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45)
ylim([40 180])
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;
% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 sg.Position=[ 2500 100 400 600]

%Figure of MODIDX SG percentage change only in rem
pcnt_rem_Data= [pcntModChange.SG_rem_modIdx];
groupCohorts=[pcntChange.Cohort];
sg=figure(105);
boxplot(pcnt_rem_Data, {groupCohorts}, 'FactorSeparator', 1, ...
        'Labels', {'Control',  'ArchT'});
    title('Percentage change D1/D2 SG MOD IDX Power');
    ylabel('Percentage (%)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45);
ylim([30 210]);
yticks([30:45:210]);
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;
% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 sg.Position=[ 2500 100 400 600]

 %Figure of MODIDX SG percentage change only in rem
pcnt_rem_Data= [pcntModChange.FG_rem_modIdx];
groupCohorts=[pcntChange.Cohort];
sg=figure(106);
boxplot(pcnt_rem_Data, {groupCohorts}, 'FactorSeparator', 1, ...
        'Labels', {'Control',  'ArchT'});
    title('Percentage change D1/D2 FG MOD IDX Power');
    ylabel('Percentage (%)');set(gca, 'FontWeight','bold');set(gca, 'FontSize',18);

boxes = findobj(gca, 'Tag', 'Box');
xtickangle(45);
ylim([30 210]);
yticks([30:45:210]);
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;
% Apply colors to each box
for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), colors(i, :),'FaceAlpha', .8);
end
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 sg.Position=[ 2500 100 400 600]
%% Stats for the percentages changes in bandpowers
statBPPcntChng.remTta=fitlme(pcntChange,'Tta_rem_pwr ~1+Cohort+(1|Subject)');
statBPPcntChng.remSG=fitlme(pcntChange,'SG_rem_pwr ~1+Cohort+(1|Subject)')
statBPPcntChng.remFG=fitlme(pcntChange,'FG_rem_pwr ~1+Cohort+(1|Subject)')

%ModIdx stats
statBPPcntChng.SGModIdx=fitlme(pcntModChange,'SG_rem_modIdx ~1+Cohort+(1|Subject)');
statBPPcntChng.FGModIdx=fitlme(pcntModChange,'FG_rem_modIdx ~1+Cohort+(1|Subject)');

%% Try to put the data in a single box chart
% % Combine data for both objects
% combinedData = [D2Int.ObjATime, D2Int.ObjBTime];
% objectLabels = categorical(repelem({'Object A', 'Object B'}, size(D2Int, 1)));
% 
% % Expand pre_post and cohort variables for both objects
% pre_post_combined = repmat([D2Int.Pre_post], [2,1]);
% cohort_combined = repelem(D2Int.Cohort, 2);
% 
% % Plotting
% figure
% clf
% h = boxchart(pre_post_combined, combinedData(:),'GroupByColor', objectLabels);
% ylabel('Time (s)'); xlabel(''); title('Time spent with Objects A and B');
% legend({'Control', 'ArchT'}, 'Location', 'Best');
% 
% % Set box colors
% h(1).BoxFaceColor = BC_color_genertor('oxford_blue');
% h(2).BoxFaceColor = BC_color_genertor('oxford_blue');
% h(3).BoxFaceColor = BC_color_genertor('ArchT_green');
% h(4).BoxFaceColor = BC_color_genertor('ArchT_green');
% 
% % Set x-tick positions and labels
% xticks([1 2]); % Set x-tick positions for "Pre" and "Post"
% xticklabels({'Pre', 'Post'}); % Set x-tick labels
% 
% hold on
% xline(1.5, 'k--'); % Line to separate "Pre" and "Post"
% set(gca,'fontsize', 16);

%% Collect and make table for sleeping
% stats = []; 
% 
% subject = []; 
% cohort = []; 
% pre_post = []; 
% no_obj_a_int = []; 
% no_obj_b_int=[];
% obj_a_time= [];
% obj_b_time= [];
% 
% all_files=fieldnames(out);
% control_list= {'BC011' 'BC013'  'BC014'};
% archT_list= {'BC051' 'BC054' 'BC1807'  'BC053'};
% controlidx=[];
% archTidx=[];
% 
% for idx = 1:numel(all_files)
%     current_field = all_files{idx};
%     if any(strcmp(current_field, control_list))
%         controlidx = [controlidx, idx];
%     end
%     if any(strcmp(current_field, archT_list))
%         archTidx = [archTidx, idx];
%     end
% end
% 
% 
% for iSub = 1:length(all_files)
%     current_subject = all_files{iSub};
% 
%     subject= [subject repmat(iSub,1,2)];
% 
%     if any(strcmp(current_subject, control_list))
%         cohort = [cohort, repmat(1,1,2)];
%     end
%     if any(strcmp(current_subject, archT_list))
%         cohort = [cohort, repmat(2,1,2)];
%     end
% 
%     pre_post = [pre_post 1 2];
%     no_obj_a_int= [no_obj_a_int out.(all_files{iSub}).D2.Behavior1.nObjAInteractions];
%     no_obj_a_int= [no_obj_a_int out.(all_files{iSub}).D2.Behavior2.nObjAInteractions];
% 
%     no_obj_b_int= [no_obj_b_int out.(all_files{iSub}).D2.Behavior1.nObjBInteractions];
%     no_obj_b_int= [no_obj_b_int out.(all_files{iSub}).D2.Behavior2.nObjBInteractions];
% 
%     obj_a_time= [obj_a_time out.(all_files{iSub}).D2.Behavior1.ObjATime];
%     obj_a_time= [obj_a_time out.(all_files{iSub}).D2.Behavior2.ObjATime];
% 
%     obj_b_time= [obj_b_time out.(all_files{iSub}).D2.Behavior1.ObjBTime];
%     obj_b_time= [obj_b_time out.(all_files{iSub}).D2.Behavior2.ObjBTime];
% 
% end
% 
% 
% %% Plotting
% figure(1)
% boxchart(D2Int.Pre_post, D2Int.ObjATime, 'GroupByColor', D2Int.Cohort, 'ColorGroup', {'red', 'blue'});
% ylabel('Time (s)'); xlabel(''); title('Time spent with object A');
% legend({'Control', 'ArchT'});

%% Comodulation analysis
% Initializing the structure containing the comodulograms
CoMoGen=[];
ArchTList={};
CtrList={};
%Creates a list of the control and experimental subjects that are present in the dataset
for iSub = 1:length(all_files)
    current_subject = all_files{iSub};
    
    if any(strcmp(current_subject,archT_list ))
        ArchTList = [ArchTList; current_subject];
    elseif any(strcmp(current_subject, control_list))
        CtrList = [CtrList; current_subject];
    end
end
%This creates the CoMo structures dividing into control and D1 and D2 for
%the sleep analysis
for iArchT =1:length(ArchTList);
    
    %The four lines below are for out-24-may-2024 due to an error in the
    %output
    % CoMoGen.ArchT.D1.rem(:,:,iArchT)= out.(ArchTList{iArchT}).D1.sleep.CoMo.CoMoSws;
    % CoMoGen.ArchT.D1.awk(:,:,iArchT)= out.(ArchTList{iArchT}).D1.sleep.CoMo.CoMoAwk;
    % CoMoGen.ArchT.D2.rem(:,:,iArchT)= out.(ArchTList{iArchT}).D2.sleep.CoMo.CoMoSws;
    % CoMoGen.ArchT.D2.awk(:,:,iArchT)= out.(ArchTList{iArchT}).D2.sleep.CoMo.CoMoAwk;

    %These other lines are for any other utput 
    CoMoGen.ArchT.D1.awk(:,:,iArchT)= out.(ArchTList{iArchT}).D1.sleep.CoMo.CoMoAwk;
    CoMoGen.ArchT.D1.sws(:,:,iArchT)= out.(ArchTList{iArchT}).D1.sleep.CoMo.CoMoSws;
    CoMoGen.ArchT.D1.rem(:,:,iArchT)= out.(ArchTList{iArchT}).D1.sleep.CoMo.CoMoRem;
    CoMoGen.ArchT.D2.awk(:,:,iArchT)= out.(ArchTList{iArchT}).D2.sleep.CoMo.CoMoAwk;
    CoMoGen.ArchT.D2.sws(:,:,iArchT)= out.(ArchTList{iArchT}).D2.sleep.CoMo.CoMoSws;
    CoMoGen.ArchT.D2.rem(:,:,iArchT)= out.(ArchTList{iArchT}).D2.sleep.CoMo.CoMoRem;
end
for iCtr =1:length(CtrList);
    %The four lines below are for out-24-may-2024 due to an error in the
    %output
    % CoMoGen.Ctr.D1.rem(:,:,iArchT)= out.(CtrList{iArchT}).D1.sleep.CoMo.CoMoSws;
    % CoMoGen.Ctr.D1.awk(:,:,iArchT)= out.(CtrList{iArchT}).D1.sleep.CoMo.CoMoAwk;
    % CoMoGen.Ctr.D2.rem(:,:,iArchT)= out.(CtrList{iArchT}).D2.sleep.CoMo.CoMoSws;
    % CoMoGen.Ctr.D2.awk(:,:,iArchT)= out.(CtrList{iArchT}).D2.sleep.CoMo.CoMoAwk;

    %These other lines are for any other utput
    CoMoGen.Ctr.D1.awk(:,:,iArchT)= out.(CtrList{iArchT}).D1.sleep.CoMo.CoMoAwk;
    CoMoGen.Ctr.D1.sws(:,:,iArchT)= out.(CtrList{iArchT}).D1.sleep.CoMo.CoMoSws;
    CoMoGen.Ctr.D1.rem(:,:,iArchT)= out.(CtrList{iArchT}).D1.sleep.CoMo.CoMoRem;
    CoMoGen.Ctr.D2.awk(:,:,iArchT)= out.(CtrList{iArchT}).D2.sleep.CoMo.CoMoAwk;
    CoMoGen.Ctr.D2.sws(:,:,iArchT)= out.(CtrList{iArchT}).D2.sleep.CoMo.CoMoSws;
    CoMoGen.Ctr.D2.rem(:,:,iArchT)= out.(CtrList{iArchT}).D2.sleep.CoMo.CoMoRem;
end
%Calculates the avg of the CoMo for all the subjects
ArchTCoMoAvg=[];
CtrCoMoAvg=[];
sleepPhsList= fieldnames(CoMoGen.ArchT.D1);
for iSP= 1:length(sleepPhsList);
    ArchTCoMoAvg.D1.(sleepPhsList{iSP})=nanmean(CoMoGen.ArchT.D1.(sleepPhsList{iSP}),3);
    ArchTCoMoAvg.D2.(sleepPhsList{iSP})=nanmean(CoMoGen.ArchT.D2.(sleepPhsList{iSP}),3);
end
for iSP= 1:length(sleepPhsList);
    CtrCoMoAvg.D1.(sleepPhsList{iSP})=nanmean(CoMoGen.Ctr.D1.(sleepPhsList{iSP}),3);
    CtrCoMoAvg.D2.(sleepPhsList{iSP})=nanmean(CoMoGen.Ctr.D2.(sleepPhsList{iSP}),3);
end
%% Plotting the Comodulograms averages
%ArchT
fg=figure(1001);
clf;
phi_f=[4 12];
amp_f=[30 100];
phi_step=0.5;
amp_step=2;
phi_f=phi_f(1):phi_step:phi_f(2);
amp_f=amp_f(1):amp_step:amp_f(2);
canva=0;
%Loop over the days
for iD=1:2;
    %Loop over the sleep phases
    for iSP=1:length(sleepPhsList)
        canva=canva+1;
        ThisCoMo=ArchTCoMoAvg.("D"+iD).(sleepPhsList{iSP});
        pnls(canva)=subplot(2,length(sleepPhsList),canva);cla;
        imagesc(phi_f, amp_f, ThisCoMo');
        set(gca, 'ydir', 'normal')
        title(sleep_phases{iSP});
        xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)');
        colorbar('Location', 'southoutside');
    end
end

ttl=annotation('textbox',[.45 0.9 0.1 0.1],'String','Average PAC ArchT ', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
D1_lbl=annotation('textbox',[.05 0.72 0.1 0.1],'String','D1', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
D2_lbl=annotation('textbox',[.05 0.25 0.1 0.1],'String','D2', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
fg.Position=[ 100 100 1000 700]

set(pnls(1),'CLim',[0 3*10^-4]);
set(pnls(4),'CLim',[0 3*10^-4]);
set(pnls(2),'CLim',[0 3.5*10^-4]);
set(pnls(5),'CLim',[0 3.5*10^-4]);
set(pnls(3),'CLim',[0 10*10^-4]);
set(pnls(6),'CLim',[0 10*10^-4]);

%Initial configuration used
% set(pnls(1),'CLim',[0 4*10^-4]);
% set(pnls(4),'CLim',[0 4*10^-4]);
% set(pnls(2),'CLim',[0 2.5*10^-3]);
% set(pnls(5),'CLim',[0 2.5*10^-3]);
% set(pnls(3),'CLim',[0 11*10^-4]);
% set(pnls(6),'CLim',[0 11*10^-4]);

%Control
fg=figure(1002)
clf;
phi_f=[4 12];
amp_f=[30 100];
phi_step=0.5;
amp_step=2;
phi_f=phi_f(1):phi_step:phi_f(2);
amp_f=amp_f(1):amp_step:amp_f(2);
canva=0;
%Loop over the days
for iD=1:2;
    %Loop over the sleep phases
    for iSP=1:length(sleepPhsList)
        canva=canva+1;
        ThisCoMo=CtrCoMoAvg.("D"+iD).(sleepPhsList{iSP});
        pnls(canva)=subplot(2,length(sleepPhsList),canva);cla;
        imagesc(phi_f, amp_f, ThisCoMo');
        set(gca, 'ydir', 'normal')
        title(sleep_phases{iSP});
        xlabel('Phase Freq (Hz)'); ylabel('Amp Freq (Hz)');
        colorbar('Location', 'southoutside');
    end
end

ttl=annotation('textbox',[.45 0.9 0.1 0.1],'String','Average PAC control', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
D1_lbl=annotation('textbox',[.05 0.72 0.1 0.1],'String','D1', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
D2_lbl=annotation('textbox',[.05 0.25 0.1 0.1],'String','D2', 'FontSize', 16, 'FontWeight', 'bold', 'LineStyle','none')
fg.Position=[ 100 100 1000 700]
set(pnls(1),'CLim',[0 2.5*10^-4]);
set(pnls(4),'CLim',[0 2.5*10^-4]);
set(pnls(2),'CLim',[0 3*10^-4]);
set(pnls(5),'CLim',[0 3*10^-4]);
set(pnls(3),'CLim',[0 9*10^-4]);
set(pnls(6),'CLim',[0 9*10^-4]);
