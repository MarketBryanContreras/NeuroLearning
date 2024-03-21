%% Plot to compare time and # of interactions
cd(inter_dir);
load ("out-20-Mar-2024.mat");
% Collect and make table for D2 interactions

stats = []; 

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


%% Plotting
%Comparison of interactions with object A
figure(1)
clf
h = boxchart(D2Int.Pre_post, D2Int.ObjATime, 'GroupByColor', D2Int.Cohort);
ylabel('Time (s)'); xlabel(''); title('Time spent with object A (moved)');
legend({'Control', 'ArchT'});

% Set box colors
h(2).BoxFaceColor = BC_color_genertor('ArchT_green');
h(1).BoxFaceColor = BC_color_genertor('oxford_blue');

% % Set x-tick positions and labels
xticks([1 2]); % Set x-tick positions for "Pre" and "Post"
xticklabels({'Pre', 'Post'}); % Set x-tick labels

hold on
xline(1.5, 'k--'); % Line to separate "Pre" and "Post"
legend({'Control', 'ArchT', ''});
set(gca,'fontsize', 16) 

%Comparison of interactions with object B
figure (2)
clf
h = boxchart(D2Int.Pre_post,D2Int.ObjBTime,'GroupByColor',D2Int.Cohort)
ylabel('Time (s)'); xlabel(''); title('Time spent with object B');
legend({'Control', 'ArchT'});

% Set box colors
h(2).BoxFaceColor = BC_color_genertor('ArchT_green');
h(1).BoxFaceColor = BC_color_genertor('oxford_blue');

% % Set x-tick positions and labels
xticks([1 2]); % Set x-tick positions for "Pre" and "Post"
xticklabels({'Pre', 'Post'}); % Set x-tick labels

hold on
xline(1.5, 'k--'); % Line to separate "Pre" and "Post"
legend({'Control', 'ArchT', ''});
set(gca,'fontsize', 16) 



% %% Collect and make table for sleeping
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