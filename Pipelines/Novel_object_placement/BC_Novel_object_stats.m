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
%% Plot to compare time and # of interactions
cd(inter_dir);
load ("out-20-Mar-2024.mat");
% Collect and make table for D2 interactions


subject = []; 
cohort = []; 
pre_post = []; 
no_obj_a_int = []; 
no_obj_b_int=[];
obj_a_time= [];
obj_b_time= [];

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

%% Collect data for the time on sleep 
%Initialize columns
 
cohort = [];
day=[];
time_awk=[];
time_sws=[];
time_rem=[];
pcnt_awk=[];
pcnt_sws=[];
pcnt_rem=[];


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

end

%% Putting it in a table
sleep= table(subject', cohort', day', time_awk', time_sws', time_rem', pcnt_awk', pcnt_sws', pcnt_rem', 'VariableNames', {'Subject', 'Cohort', 'Day', 'time_awk', 'time_sws', 'time_rem', 'pcnt_awk','pcnt_sws', 'pcnt_rem'});
%% Stats

%% Plotting
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
%% Possible solution to error bars
% c = categorical({'CH','VC','GC','OC','BC','SC'});
% c = reordercats(c,{'CH','VC','GC','OC','BC','SC'});
% y = [707, 599; 464 444; 522 475; 566 346; 1329 1384; 459 498];
% std_dev = [321 271; 233 91; 202 132; 0 173; 410  850; 179 122];
% title('Title'); xlabel('x-label'); ylabel('y-label');
% figure 
% hold on
% bar(c,y)
% errorbar(y,std_dev,'.')
% 
% %Solution from stack flow
% close all
% clear all
% clc
% y = [707, 599; 464 444; 522 475; 566 346; 1329 1384; 459 498];
% std_dev = [321 271; 233 91; 202 132; 0 173; 410  850; 179 122];
% num = 6; %number of different subcategories
% c = 1:num;
% %%Figure
% figH = figure;
% axes1 = axes;
% title('Title'); xlabel('x-label'); ylabel('y-label');
% hold on
% %%Bar(s)
% %You can not color differently the same bar.
% for i = 1:num
%     bar(c(i)-0.15,y(i,1),0.2);
%     bar(c(i)+0.15,y(i,2),0.2);
% end
% %%Errorbar
% errH1 = errorbar(c-0.15,y(:,1),std_dev(:,1),'.','Color','b');
% errH2 = errorbar(c+0.15,y(:,2),std_dev(:,2),'.','Color','r');
% errH1.LineWidth = 1.5;
% errH2.LineWidth = 1.5;
% errH1.Color = [1 0.5 0];
% errH2.Color = [1 0.3 1];
% %%Set x-ticks
% set(axes1,'Xlim',[0.5 5.5]);
% set(axes1,'XTick',[1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6],'XTickLabel',...
%     {'CH',' ','VC',' ','GC',' ','OC',' ','BC',' ','SC'});
%   end
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

%% % Collect and make table for sleeping
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