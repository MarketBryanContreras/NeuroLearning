%% SUPER_BC_stats

% load(['C:\Users\ecarm\Downloads' filesep 'out_Archt_06_nov_23.mat'])
% load(['C:\Users\ecarm\Downloads' filesep 'out_eyfp_06_nov_23.mat'])
%% Adjust the name of the files to load
ArchT_file_name='out_ArchT_07_nov_23.mat';
eyfp_file_name='out_eyfp_07_nov_23.mat';
%% Load
load(['C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter' filesep ,ArchT_file_name])
load(['C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter' filesep ,eyfp_file_name])
%%
% parts = strsplit(ArchT_file_name, '.');
% archt_var= eval(parts{1})
%% collect and make table
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

archt_list = fieldnames(out_Archt_07_nov_23);
eyfp_list = fieldnames(out_eyfp_07_nov_23); 


for iSub = 1:length(archt_list)
    sess_list = fieldnames(out_Archt_07_nov_23.(archt_list{iSub})); 
 
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

%      z_SGInhb_modidx=[z_SGInhb_modidx [out_eyfp_07_nov_23.(archt_list{iSub}).D4.z_SGInhb_modidx]];
%      
%      z_FGInhb_modidx=[z_FGInhb_modidx [out_eyfp_07_nov_23.(archt_list{iSub}).D4.z_FGInhb_modidx]];
%      
%      z_FGNoInhb_modidx=[z_FGNoInhb_modidx [out_eyfp_07_nov_23.(archt_list{iSub}).D4.z_FGNoInhb_modidx]];
     
     if iSub==length(archt_list)
     archt_tbl= table(subject', cohort', opto', trial_n',SG_modidx',FG_modidx',t_bp',sg_bp',fg_bp','VariableNames', {'Subject', 'Cohort', 'Opto', 'Trial', 'SG_modidx', 'FG_modidx','Theta_bp','SG_bp','FG_bp'});
     end
end

% loop for eyfp

for iSub = 1:length(eyfp_list)
    sess_list = fieldnames(out_eyfp_07_nov_23.(eyfp_list{iSub})); 
 
    n_inhib = length(out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.t_bp_inhib); 
    n_noinhib = length(out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.t_bp_noinhib);
    
    
     subject= [subject repmat(iSub+length(archt_list),1, n_inhib + n_noinhib)]; 
     
     cohort = [cohort repmat(2,1, n_inhib + n_noinhib)];
     
     trial_n = [trial_n, [1:n_inhib, 1:n_noinhib]]; 
     
     opto = [opto, logical([repmat(1,1, n_inhib), repmat(0,1, n_noinhib)])]; 
     
     SG_modidx = [SG_modidx [out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.modidx_SG_inhib, out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.modidx_SG_noinhib]];
     
     FG_modidx = [FG_modidx [out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.modidx_FG_inhib, out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.modidx_FG_noinhib]];
     
     t_bp = [t_bp [out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.t_bp_inhib, out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.t_bp_noinhib]];
     
     sg_bp = [sg_bp [out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.sg_bp_inhib, out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.sg_bp_noinhib]];
     
     fg_bp = [fg_bp [out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.fg_bp_inhib, out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.fg_bp_noinhib]];
     
     if iSub==length(eyfp_list)
         eyfp_tbl= table(subject', cohort', opto', trial_n',SG_modidx',FG_modidx',t_bp',sg_bp',fg_bp','VariableNames', {'Subject', 'Cohort', 'Opto', 'Trial', 'SG_modidx', 'FG_modidx','Theta_bp','SG_bp','FG_bp'});
         eyfp_tbl= eyfp_tbl((size(archt_tbl,1)+1):size(eyfp_tbl,1),:);
     end
end

% opto_cells  = [];
% opto_cells(opto == 1) = 'Inhib'; 
% opto_cells(opto == 0) = 'Inhib'; 
opto_cells  =cell(size(opto));
opto_cells(opto == 1) = {'Light'}; 
opto_cells(opto == 0) = {'No Light'}; 

%% collect all in one table
tbl= table(subject', cohort', opto', trial_n',SG_modidx',FG_modidx',t_bp',sg_bp',fg_bp','VariableNames', {'Subject', 'Cohort', 'Opto', 'Trial', 'SG_modidx', 'FG_modidx','Theta_bp','SG_bp','FG_bp'});
%% simple stats
stats.SG_mod = fitlme(tbl,'SG_modidx ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
stats.FG_mod = fitlme(tbl,'FG_modidx ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
stats.t_bp = fitlme(tbl,'Theta_bp ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
stats.sg_bp = fitlme(tbl,'SG_bp ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
stats.fg_bp = fitlme(tbl,'FG_bp ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');

archt_stats.SG_mod = fitlme(archt_tbl,'SG_modidx ~ 1+ Opto + (1|Subject) + (1|Trial)');
archt_stats.FG_mod = fitlme(archt_tbl,'FG_modidx ~ 1+ Opto + (1|Subject) + (1|Trial)');
archt_stats.t_bp =   fitlme(archt_tbl,'Theta_bp ~ 1+ Opto + (1|Subject) + (1|Trial)');
archt_stats.sg_bp =   fitlme(archt_tbl,'SG_bp ~ 1+ Opto + (1|Subject)+ (1|Trial)');
archt_stats.fg_bp =   fitlme(archt_tbl,'FG_bp ~ 1+ Opto + (1|Subject)+ (1|Trial)');

eyfp_stats.SG_mod = fitlme(eyfp_tbl,'SG_modidx ~ 1+ Opto + (1|Subject) + (1|Trial)');
eyfp_stats.FG_mod = fitlme(eyfp_tbl,'FG_modidx ~ 1+ Opto + (1|Subject) + (1|Trial)');
eyfp_stats.t_bp =   fitlme(eyfp_tbl,'Theta_bp ~ 1+ Opto + (1|Subject) + (1|Trial)');
eyfp_stats.sg_bp =  fitlme(eyfp_tbl,'SG_bp ~ 1+ Opto + (1|Subject)+ (1|Trial)');
eyfp_stats.fg_bp =   fitlme(eyfp_tbl,'FG_bp ~ 1+ Opto + (1|Subject)+ (1|Trial)');

stats.SG_mod = fitglme(tbl,'SG_modidx ~ 1+ Opto + (1|Subject)');
%% Plotting
figure(102)
subplot(1,2,1)
%boxchart(tbl.Opto,tbl.SG_modidx,'GroupByColor',tbl.Opto)
boxplot( tbl.SG_modidx,tbl.Opto)

subplot(1,2,2)
%boxchart(tbl.Subject,tbl.FG_modidx,'GroupByColor',tbl.Opto)
boxplot( tbl.FG_modidx,tbl.Opto )



