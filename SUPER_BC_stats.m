%% SUPER_BC_stats

load(['C:\Users\ecarm\Downloads' filesep 'out_Archt_06_nov_23.mat'])
load(['C:\Users\ecarm\Downloads' filesep 'out_eyfp_06_nov_23.mat'])


%% collect and make table
tbl = table(); 
stats = []; 
subject = []; 
trial_n = [];
opto = [];
cohort = []; 
SG_modidx = []; 
FG_modidx = []; 

archt_list = fieldnames(out_Archt_06_nov_23);
eyfp_list = fieldnames(out_eyfp_06_nov_23); 


for iSub = 1:length(archt_list)
    sess_list = fieldnames(out_Archt_06_nov_23.(archt_list{iSub})); 
 
    n_inhib = length(out_Archt_06_nov_23.(archt_list{iSub}).D4.t_bp_inhib); 
    n_noinhib = length(out_Archt_06_nov_23.(archt_list{iSub}).D4.t_bp_noinhib);
    
    
     subject= [subject repmat(iSub,1, n_inhib + n_noinhib)]; 
     
     cohort = [cohort repmat(1,1, n_inhib + n_noinhib)];
     
     trial_n = [trial_n, [1:n_inhib, 1:n_noinhib]]; 
     
     opto = [opto, logical([repmat(1,1, n_inhib), repmat(0,1, n_noinhib)])]; 
     
     SG_modidx = [SG_modidx [out_Archt_06_nov_23.(archt_list{iSub}).D4.modidx_SG_inhib, out_Archt_06_nov_23.(archt_list{iSub}).D4.modidx_SG_noinhib]];
     
     FG_modidx = [FG_modidx [out_Archt_06_nov_23.(archt_list{iSub}).D4.modidx_FG_inhib, out_Archt_06_nov_23.(archt_list{iSub}).D4.modidx_FG_noinhib]];

end

% loop for eyfp

for iSub = 1:length(eyfp_list)
    sess_list = fieldnames(out_eyfp_06_nov_23.(eyfp_list{iSub})); 
 
    n_inhib = length(out_eyfp_06_nov_23.(eyfp_list{iSub}).D4.t_bp_inhib); 
    n_noinhib = length(out_eyfp_06_nov_23.(eyfp_list{iSub}).D4.t_bp_noinhib);
    
    
     subject= [subject repmat(iSub+length(archt_list),1, n_inhib + n_noinhib)]; 
     
     cohort = [cohort repmat(2,1, n_inhib + n_noinhib)];
     
     trial_n = [trial_n, [1:n_inhib, 1:n_noinhib]]; 
     
     opto = [opto, logical([repmat(1,1, n_inhib), repmat(0,1, n_noinhib)])]; 
     
     SG_modidx = [SG_modidx [out_eyfp_06_nov_23.(eyfp_list{iSub}).D4.modidx_SG_inhib, out_eyfp_06_nov_23.(eyfp_list{iSub}).D4.modidx_SG_noinhib]];
     
     FG_modidx = [FG_modidx [out_eyfp_06_nov_23.(eyfp_list{iSub}).D4.modidx_FG_inhib, out_eyfp_06_nov_23.(eyfp_list{iSub}).D4.modidx_FG_noinhib]];

end


opto_cells  =cell(size(opto));
opto_cells{opto == 1} = 'Light'; 
opto_cells{opto == 0} = 'No Light'; 

%% collect all in one table
tbl= table(subject', cohort', opto', trial_n',SG_modidx',FG_modidx', 'VariableNames', {'Subject', 'Cohort', 'Opto', 'Trial', 'SG_modidx', 'FG_modidx'});

figure(102)
subplot(1,2,1)
boxchart(tbl.opto,tbl.SG_modidx,'GroupByColor',tbl.Opto)


subplot(1,2,2)
boxchart(tbl.Subject,tbl.FG_modidx,'GroupByColor',tbl.Opto)

% simple stats
stats.SG_mod = fitlme(tbl,'SG_modidx ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
stats.SG_mod = fitglme(tbl,'SG_modidx ~ 1+ Opto + (1|Subject)');

% stats.SG_mod_yfp = fitlme(tbl,'SG_modidx ~ 1+ Opto + (1|Subject) + (1|Trial)');

stats.SG_mod_coh = fitlme(tbl,'SG_modidx ~ 1+ Opto*Cohort + (1|Subject) + (1|Trial)');



stats.FG_mod = fitlme(tbl,'FG_modidx ~ 1+ Opto + (1|Subject) + (1|Trial)');




