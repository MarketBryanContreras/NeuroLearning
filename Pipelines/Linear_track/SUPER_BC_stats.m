%% SUPER_BC_stats

% load(['C:\Users\ecarm\Downloads' filesep 'out_Archt_06_nov_23.mat'])
% load(['C:\Users\ecarm\Downloads' filesep 'out_eyfp_06_nov_23.mat'])
%% Adjust the name of the files to load
ArchT_file_name='out_Arch_07_nov_23.mat';
eyfp_file_name='out_eyfp_07_nov_23.mat';
%% Dynamic loader
sys=computer;
if contains(sys,'PCWIN')
    % Dynamic based on computer user windows.
    out_archt=load([getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter' filesep ,ArchT_file_name])
    out_eyfp=load([getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter' filesep ,eyfp_file_name])
elseif contains(sys,'MAC')
    % Dynamic based on computer user for mac
    out_archt=load([getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'inter' filesep ,ArchT_file_name]);
    out_eyfp=load([getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'inter' filesep ,eyfp_file_name]);
    inter_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'inter'];
else disp("This system is not supported on the dynamic loader yet, pleae add it or contact BC for it")
end
out_archt=out_arch.out;
out_eyfp=out_eyfp.out;
% out_archt=out_archt.out_Arch_07_nov_23;
% out_eyfp=out_eyfp.out_eyfp_07_nov_23;
archt_list = fieldnames(out_archt);
eyfp_list = fieldnames(out_eyfp); 
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

for iSub = 1:length(archt_list)
    sess_list = fieldnames(out_archt.(archt_list{iSub})); 
 
    n_inhib = length(out_archt.(archt_list{iSub}).D4.t_bp_inhib); 
    n_noinhib = length(out_archt.(archt_list{iSub}).D4.t_bp_noinhib);
    
    
     subject= [subject repmat(iSub,1, n_inhib + n_noinhib)]; 
     
     cohort = [cohort repmat(1,1, n_inhib + n_noinhib)];
     
     trial_n = [trial_n, [1:n_inhib, 1:n_noinhib]]; 
     
     opto = [opto, logical([repmat(1,1, n_inhib), repmat(0,1, n_noinhib)])]; 
     
     SG_modidx = [SG_modidx [out_archt.(archt_list{iSub}).D4.modidx_SG_inhib, out_archt.(archt_list{iSub}).D4.modidx_SG_noinhib]];
     
     FG_modidx = [FG_modidx [out_archt.(archt_list{iSub}).D4.modidx_FG_inhib, out_archt.(archt_list{iSub}).D4.modidx_FG_noinhib]];

     t_bp = [t_bp [out_archt.(archt_list{iSub}).D4.t_bp_inhib, out_archt.(archt_list{iSub}).D4.t_bp_noinhib]];
     
     sg_bp = [sg_bp [out_archt.(archt_list{iSub}).D4.sg_bp_inhib, out_archt.(archt_list{iSub}).D4.sg_bp_noinhib]];
     
     fg_bp = [fg_bp [out_archt.(archt_list{iSub}).D4.fg_bp_inhib, out_archt.(archt_list{iSub}).D4.fg_bp_noinhib]];

%      z_SGInhb_modidx=[z_SGInhb_modidx [out_Archt_07_nov_23.(archt_list{iSub}).D4.z_SGInhb_modidx]];
%      
%      z_FGInhb_modidx=[z_FGInhb_modidx [out_Archt_07_nov_23.(archt_list{iSub}).D4.z_FGInhb_modidx]];
%      
%      z_FGNoInhb_modidx=[z_FGNoInhb_modidx [out_Archt_07_nov_23.(archt_list{iSub}).D4.z_FGNoInhb_modidx]];
     
     if iSub==length(archt_list)
     archt_tbl= table(subject', cohort', opto', trial_n',SG_modidx',FG_modidx',t_bp',sg_bp',fg_bp','VariableNames', {'Subject', 'Cohort', 'Opto', 'Trial', 'SG_modidx', 'FG_modidx','Theta_bp','SG_bp','FG_bp'});
     end
end

% loop for eyfp

for iSub = 1:length(eyfp_list)
    sess_list = fieldnames(out_eyfp.(eyfp_list{iSub})); 
 
    n_inhib = length(out_eyfp.(eyfp_list{iSub}).D4.t_bp_inhib); 
    n_noinhib = length(out_eyfp.(eyfp_list{iSub}).D4.t_bp_noinhib);
    
    
     subject= [subject repmat(iSub+length(archt_list),1, n_inhib + n_noinhib)]; 
     
     cohort = [cohort repmat(2,1, n_inhib + n_noinhib)];
     
     trial_n = [trial_n, [1:n_inhib, 1:n_noinhib]]; 
     
     opto = [opto, logical([repmat(1,1, n_inhib), repmat(0,1, n_noinhib)])]; 
     
     SG_modidx = [SG_modidx [out_eyfp.(eyfp_list{iSub}).D4.modidx_SG_inhib, out_eyfp.(eyfp_list{iSub}).D4.modidx_SG_noinhib]];
     
     FG_modidx = [FG_modidx [out_eyfp.(eyfp_list{iSub}).D4.modidx_FG_inhib, out_eyfp.(eyfp_list{iSub}).D4.modidx_FG_noinhib]];
     
     t_bp = [t_bp [out_eyfp.(eyfp_list{iSub}).D4.t_bp_inhib, out_eyfp.(eyfp_list{iSub}).D4.t_bp_noinhib]];
     
     sg_bp = [sg_bp [out_eyfp.(eyfp_list{iSub}).D4.sg_bp_inhib, out_eyfp.(eyfp_list{iSub}).D4.sg_bp_noinhib]];
     
     fg_bp = [fg_bp [out_eyfp.(eyfp_list{iSub}).D4.fg_bp_inhib, out_eyfp.(eyfp_list{iSub}).D4.fg_bp_noinhib]];
     
     if iSub==length(eyfp_list)
         eyfp_tbl= table(subject', cohort', opto', trial_n',SG_modidx',FG_modidx',t_bp',sg_bp',fg_bp','VariableNames', {'Subject', 'Cohort', 'Opto', 'Trial', 'SG_modidx', 'FG_modidx','Theta_bp','SG_bp','FG_bp'});
         eyfp_tbl= eyfp_tbl((size(archt_tbl,1)+1):size(eyfp_tbl,1),:);
     end
end


% opto_cells  = [];
% opto_cells(opto == 1) = 'Inhib'; 
% opto_cells(opto == 0) = 'Inhib'; 
% opto_cells  =cell(size(opto));
% opto_cells(opto == 1) = {'Light'}; 
% opto_cells(opto == 0) = {'No Light'}; 


%% collect all in one table
tbl= table(subject', cohort', opto', trial_n',SG_modidx',FG_modidx',t_bp',sg_bp',fg_bp','VariableNames', {'Subject', 'Cohort', 'Opto', 'Trial', 'SG_modidx', 'FG_modidx','Theta_bp','SG_bp','FG_bp'});
%% simple stats
stats.SG_mod = fitlme(tbl,'SG_modidx ~ 1+ Opto*Cohort + (1|Subject) + (1|Trial)');
stats.FG_mod = fitlme(tbl,'FG_modidx ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
stats.t_bp = fitlme(tbl,'Theta_bp ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
stats.sg_bp = fitlme(tbl,'SG_bp ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');
stats.fg_bp = fitlme(tbl,'FG_bp ~ 1+ Opto + (1|Cohort) + (1|Subject) + (1|Trial)');

archt_stats.SG_mod = fitlme(archt_tbl,'SG_modidx ~ 1+ Opto+ (1|Subject) + (1|Trial)');
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

boxchart(tbl.Subject,tbl.SG_modidx,'GroupByColor',tbl.Opto)
boxchart(tbl.Subject,tbl.FG_modidx,'GroupByColor',tbl.Opto)
boxchart(tbl.Cohort,tbl.SG_modidx,'GroupByColor',tbl.Opto )
%% ArchT
% Filtering SG Data
archt_sg_modidx_values_opto_1 = tbl.SG_modidx(tbl.Opto == 1 & tbl.Cohort==1);
archt_sg_modidx_values_opto_0 = tbl.SG_modidx(tbl.Opto == 0 & tbl.Cohort==1);

archt_sg_xx=[archt_sg_modidx_values_opto_1' archt_sg_modidx_values_opto_0'];

archt_sg_aa=repmat(1,1,length(archt_sg_modidx_values_opto_1));
archt_sg_bb=repmat(0,1,length(archt_sg_modidx_values_opto_0));

archt_sg_yy=[archt_sg_aa archt_sg_bb];
% Filtering FG Data
archt_fg_modidx_values_opto_1 = tbl.FG_modidx(tbl.Opto == 1 & tbl.Cohort==1);
archt_fg_modidx_values_opto_0 = tbl.FG_modidx(tbl.Opto == 0 & tbl.Cohort==1);

archt_fg_xx=[archt_fg_modidx_values_opto_1' archt_fg_modidx_values_opto_0'];

archt_fg_aa=repmat(1,1,length(archt_fg_modidx_values_opto_1));
archt_fg_bb=repmat(0,1,length(archt_fg_modidx_values_opto_0));
archt_fg_yy=[archt_fg_aa archt_fg_bb];

% Filtering Theta_bp Data
archt_theta_bp_values_opto_1 = tbl.Theta_bp(tbl.Opto == 1 & tbl.Cohort==1);
archt_theta_bp_values_opto_0 = tbl.Theta_bp(tbl.Opto == 0 & tbl.Cohort==1);

archt_theta_bp_xx=[archt_theta_bp_values_opto_1' archt_theta_bp_values_opto_0'];

archt_theta_bp_aa=repmat(1,1,length(archt_theta_bp_values_opto_1));
archt_theta_bp_bb=repmat(0,1,length(archt_theta_bp_values_opto_0));

archt_theta_bp_yy=[archt_theta_bp_aa archt_theta_bp_bb];

% Filtering SG_bp Data
archt_sg_bp_values_opto_1 = tbl.SG_bp(tbl.Opto == 1 & tbl.Cohort==1);
archt_sg_bp_values_opto_0 = tbl.SG_bp(tbl.Opto == 0 & tbl.Cohort==1);

archt_sg_bp_xx=[archt_sg_bp_values_opto_1' archt_sg_bp_values_opto_0'];

archt_sg_bp_aa=repmat(1,1,length(archt_sg_bp_values_opto_1));
archt_sg_bp_bb=repmat(0,1,length(archt_sg_bp_values_opto_0));

archt_sg_bp_yy=[archt_sg_bp_aa archt_sg_bp_bb];

% Filtering FG_bp Data
archt_fg_bp_values_opto_1 = tbl.FG_bp(tbl.Opto == 1 & tbl.Cohort==1);
archt_fg_bp_values_opto_0 = tbl.FG_bp(tbl.Opto == 0 & tbl.Cohort==1);

archt_fg_bp_xx=[archt_fg_bp_values_opto_1' archt_fg_bp_values_opto_0'];

archt_fg_bp_aa=repmat(1,1,length(archt_fg_bp_values_opto_1));
archt_fg_bp_bb=repmat(0,1,length(archt_fg_bp_values_opto_0));

archt_fg_bp_yy=[archt_fg_bp_aa archt_fg_bp_bb];

%% eyfp
% Filtering SG Data
eyfp_sg_modidx_values_opto_1 = tbl.SG_modidx(tbl.Opto == 1 & tbl.Cohort==2);
eyfp_sg_modidx_values_opto_0 = tbl.SG_modidx(tbl.Opto == 0 & tbl.Cohort==2);

eyfp_sg_xx=[eyfp_sg_modidx_values_opto_1' eyfp_sg_modidx_values_opto_0'];

eyfp_sg_aa=repmat(1,1,length(eyfp_sg_modidx_values_opto_1));
eyfp_sg_bb=repmat(0,1,length(eyfp_sg_modidx_values_opto_0));

eyfp_sg_yy=[eyfp_sg_aa eyfp_sg_bb];
% Filtering FG Data
eyfp_fg_modidx_values_opto_1 = tbl.FG_modidx(tbl.Opto == 1 & tbl.Cohort==2);
eyfp_fg_modidx_values_opto_0 = tbl.FG_modidx(tbl.Opto == 0 & tbl.Cohort==2);

eyfp_fg_xx=[eyfp_fg_modidx_values_opto_1' eyfp_fg_modidx_values_opto_0'];

eyfp_fg_aa=repmat(1,1,length(eyfp_fg_modidx_values_opto_1));
eyfp_fg_bb=repmat(0,1,length(eyfp_fg_modidx_values_opto_0));
eyfp_fg_yy=[eyfp_fg_aa eyfp_fg_bb];

% Filtering Theta_bp Data
eyfp_theta_bp_values_opto_1 = tbl.Theta_bp(tbl.Opto == 1 & tbl.Cohort==2);
eyfp_theta_bp_values_opto_0 = tbl.Theta_bp(tbl.Opto == 0 & tbl.Cohort==2);

eyfp_theta_bp_xx=[eyfp_theta_bp_values_opto_1' eyfp_theta_bp_values_opto_0'];

eyfp_theta_bp_aa=repmat(1,1,length(eyfp_theta_bp_values_opto_1));
eyfp_theta_bp_bb=repmat(0,1,length(eyfp_theta_bp_values_opto_0));

eyfp_theta_bp_yy=[eyfp_theta_bp_aa eyfp_theta_bp_bb];

% Filtering SG_bp Data
eyfp_sg_bp_values_opto_1 = tbl.SG_bp(tbl.Opto == 1 & tbl.Cohort==2);
eyfp_sg_bp_values_opto_0 = tbl.SG_bp(tbl.Opto == 0 & tbl.Cohort==2);

eyfp_sg_bp_xx=[eyfp_sg_bp_values_opto_1' eyfp_sg_bp_values_opto_0'];

eyfp_sg_bp_aa=repmat(1,1,length(eyfp_sg_bp_values_opto_1));
eyfp_sg_bp_bb=repmat(0,1,length(eyfp_sg_bp_values_opto_0));

eyfp_sg_bp_yy=[eyfp_sg_bp_aa eyfp_sg_bp_bb];

% Filtering FG_bp Data
eyfp_fg_bp_values_opto_1 = tbl.FG_bp(tbl.Opto == 1 & tbl.Cohort==2);
eyfp_fg_bp_values_opto_0 = tbl.FG_bp(tbl.Opto == 0 & tbl.Cohort==2);

eyfp_fg_bp_xx=[eyfp_fg_bp_values_opto_1' eyfp_fg_bp_values_opto_0'];

eyfp_fg_bp_aa=repmat(1,1,length(eyfp_fg_bp_values_opto_1));
eyfp_fg_bp_bb=repmat(0,1,length(eyfp_fg_bp_values_opto_0));

eyfp_fg_bp_yy=[eyfp_fg_bp_aa eyfp_fg_bp_bb];

%% Plotting Mod IDX stats
%Archt
figure(102)
clf
subplot(2,2,1)
boxplot(archt_sg_xx,archt_sg_yy, 'Labels', {'No silencing', 'Silencing'});
%Significance plot
yt = get(gca, 'YTick');
set(gca, 'Xtick', 1:3);
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*yt(end), '-k',mean(xt([1 2])),yt(end)*1.1, '*k', mean(xt([1 2]))+0.05*mean(xt([1 2])),yt(end)*1.1, '*k',mean(xt([1 2]))-0.05*mean(xt([1 2])),yt(end)*1.1, '*k') %Three asteriks
%plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2]))-0.025*mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.025*mean(xt([1 2])),yt(5)*1.1,'*k') 2 asteriks

title('SG ArchT');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('Powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(2, 2, 2);
boxplot(archt_fg_xx,archt_fg_yy, 'Labels', {'No silencing', 'Silencing'});
title('FG  ArchT');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('Powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

%eyfp
subplot(2,2,3)
boxplot(eyfp_sg_xx,eyfp_sg_yy, 'Labels', {'No silencing', 'Silencing'});
title('SG control');
colors = [BC_color_genertor('swamp_green');  % RGB values for Group 2
    BC_color_genertor('torment_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(2, 2, 4);
boxplot(eyfp_fg_xx,eyfp_fg_yy, 'Labels', {'No silencing', 'Silencing'});

title('FG control');
colors = [BC_color_genertor('swamp_green');  % RGB values for Group 2
    BC_color_genertor('torment_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

%sgtitle('General Title');


 % Adjust figure properties
 t=title(" Mod IDX");           % Add a general title above the entire figure")
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 % Resize the figure (optional)
 fig.Position = [200, 200, 1000, 700];  % [x, y, width, height]
 %saveas(gcf, 'FG_NormModIdx.png');%%%% fill in nice plotting %%%%%%%
%% Plotting Band power stats
fig103=figure(103)
clf


subplot(2,3,1)
boxplot(archt_theta_bp_xx,archt_theta_bp_yy, 'Labels', {'No silencing', 'Silencing'});
ylabel('Power (mW)');
%Significance plot
yt = get(gca, 'YTick');
ylim([0 yt(end)*1.2])
% set(gca, 'Xtick', 1:3);
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*yt(end)*1.1, '-k',mean(xt([1 2])),yt(end)*1.15, '*k', mean(xt([1 2]))+0.05*mean(xt([1 2])),yt(end)*1.15, '*k',mean(xt([1 2]))-0.05*mean(xt([1 2])),yt(end)*1.15, '*k') %Three asteriks
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2]))-0.025*mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.025*mean(xt([1 2])),yt(5)*1.1,'*k') %2 asteriks

title('Theta ArchT');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('Powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(2, 3, 2);
boxplot(archt_sg_bp_xx,archt_sg_bp_yy, 'Labels', {'No silencing', 'Silencing'});
%%Significance plot
% yt = get(gca, 'YTick');
% set(gca, 'Xtick', 1:3);
% xt = get(gca, 'XTick');
% hold on
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.05*mean(xt([1 2])),yt(5)*1.1, '*k',mean(xt([1 2]))-0.05*mean(xt([1 2])),yt(5)*1.1, '*k') %Three asteriks
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2]))-0.025*mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.025*mean(xt([1 2])),yt(5)*1.1,'*k') 2 asteriks

title('SG ArchT');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('Powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(2, 3, 3);
boxplot(archt_fg_bp_xx,archt_fg_bp_yy, 'Labels', {'No silencing', 'Silencing'});
%Significance plot
yt = get(gca, 'YTick');
ylim([0 yt(end)*1.5])
% set(gca, 'Xtick', 1:3);
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*yt(end)*1.3, '-k',mean(xt([1 2])),yt(end)*1.35, '*k', mean(xt([1 2]))+0.05*mean(xt([1 2])),yt(end)*1.35, '*k',mean(xt([1 2]))-0.05*mean(xt([1 2])),yt(end)*1.35, '*k') %Three asteriks
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2]))-0.025*mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.025*mean(xt([1 2])),yt(5)*1.1,'*k') %2 asteriks

title('FG ArchT');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('Powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

%EYFP

subplot(2,3,4)
boxplot(eyfp_theta_bp_xx,eyfp_theta_bp_yy, 'Labels', {'No silencing', 'Silencing'});
ylabel('Power (mW)');
%%Significance plot
% yt = get(gca, 'YTick');
% set(gca, 'Xtick', 1:3);
% xt = get(gca, 'XTick');
% hold on
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.05*mean(xt([1 2])),yt(5)*1.1, '*k',mean(xt([1 2]))-0.05*mean(xt([1 2])),yt(5)*1.1, '*k') %Three asteriks
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2]))-0.025*mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.025*mean(xt([1 2])),yt(5)*1.1,'*k') 2 asteriks

title('Theta control');
colors = [BC_color_genertor('swamp_green');  % RGB values for Group 2
    BC_color_genertor('torment_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(2, 3, 5);
boxplot(eyfp_sg_bp_xx,eyfp_sg_bp_yy, 'Labels', {'No silencing', 'Silencing'});
%%Significance plot
% yt = get(gca, 'YTick');
% set(gca, 'Xtick', 1:3);
% xt = get(gca, 'XTick');
% hold on
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.05*mean(xt([1 2])),yt(5)*1.1, '*k',mean(xt([1 2]))-0.05*mean(xt([1 2])),yt(5)*1.1, '*k') %Three asteriks
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2]))-0.025*mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.025*mean(xt([1 2])),yt(5)*1.1,'*k') 2 asteriks

title('SG control');
colors = [BC_color_genertor('swamp_green');  % RGB values for Group 2
    BC_color_genertor('torment_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(2, 3, 6);
boxplot(eyfp_fg_bp_xx,eyfp_fg_bp_yy, 'Labels', {'No silencing', 'Silencing'});
%%Significance plot
% yt = get(gca, 'YTick');
% set(gca, 'Xtick', 1:3);
% xt = get(gca, 'XTick');
% hold on
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.05*mean(xt([1 2])),yt(5)*1.1, '*k',mean(xt([1 2]))-0.05*mean(xt([1 2])),yt(5)*1.1, '*k') %Three asteriks
% plot(xt([1 2]), [1 1]*yt(5), '-k',mean(xt([1 2]))-0.025*mean(xt([1 2])),yt(5)*1.1, '*k', mean(xt([1 2]))+0.025*mean(xt([1 2])),yt(5)*1.1,'*k') 2 asteriks

title('FG control');
colors = [BC_color_genertor('swamp_green');  % RGB values for Group 2
    BC_color_genertor('torment_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

% Adjust figure properties

 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 % Resize the figure (optional)
 fig.Position = [200, 200, 1000, 700];  % [x, y, width, height]
 %saveas(gcf, 'FG_NormModIdx.png');%%%% fill in nice plotting %%%%%%%




%% average value tabl

tbl_avg = table(); 
stats_avg = []; 
subject = []; 
opto = [];
cohort = []; 
SG_modidx_z = []; 
FG_modidx_z = []; 


archt_list = fieldnames(out_archt);
eyfp_list = fieldnames(out_eyfp); 


for iSub = 1:length(archt_list)
    
subject = [subject iSub iSub];

cohort = [cohort 1 1]; 
opto = [opto 1 0]; 

SG_modidx_z = [SG_modidx_z out_archt.(archt_list{iSub}).D4.z_SGInhb_modidx out_archt.(archt_list{iSub}).D4.z_SGNoInhb_modidx]; % has both inhib and no_inhib

FG_modidx_z = [FG_modidx_z out_archt.(archt_list{iSub}).D4.z_FGInhb_modidx out_archt.(archt_list{iSub}).D4.z_FGNoInhb_modidx]; % has both inhib and no_inhib

end
tbl_avg_archT= table(subject', cohort', opto',SG_modidx_z', FG_modidx_z' ,'VariableNames', {'Subject', 'Cohort', 'Opto', 'SG_modidx_z', 'FG_modidx_z'});



for iSub = 1:length(eyfp_list)
    
subject = [subject iSub+length(archt_list) iSub+length(archt_list)];

cohort = [cohort 0 0]; 
opto = [opto 1 0]; 

SG_modidx_z = [SG_modidx_z out_eyfp.(eyfp_list{iSub}).D4.z_SGInhb_modidx out_eyfp.(eyfp_list{iSub}).D4.z_SGNoInhb_modidx]; % has both inhib and no_inhib

FG_modidx_z = [FG_modidx_z out_eyfp.(eyfp_list{iSub}).D4.z_FGInhb_modidx out_eyfp.(eyfp_list{iSub}).D4.z_FGNoInhb_modidx]; % has both inhib and no_inhib

end

tbl_avg= table(subject', cohort', logical(opto'),SG_modidx_z', FG_modidx_z' ,'VariableNames', {'Subject', 'Cohort', 'Opto', 'SG_modidx_z', 'FG_modidx_z'});
archTavgtbl=tbl_avg(tbl_avg.Cohort==1,:);
ctrlavgtbl=tbl_avg(tbl_avg.Cohort==0,:);
%% stats for avg Z mod
zstats.fg=fitlme(tbl_avg,'FG_modidx_z ~ 1+ Opto + (1|Cohort) + (1|Subject)');
zstats.sg=fitlme(tbl_avg,'SG_modidx_z ~ 1+ Opto + (1|Cohort) + (1|Subject)');
zstats.archtFG=fitlme(archTavgtbl,'FG_modidx_z ~ 1+ Opto  + (1|Subject)');
zstats.archtSG=fitlme(archTavgtbl,'SG_modidx_z ~ 1+ Opto  + (1|Subject)');
zstats.ctrlFG=fitlme(ctrlavgtbl,'FG_modidx_z ~ 1+ Opto  + (1|Subject)');
zstats.ctrlSG=fitlme(ctrlavgtbl,'SG_modidx_z ~ 1+ Opto  + (1|Subject)');

%% Preparation for the graph below
% subplot(1,2,1)
% boxchart(tbl.Opto,tbl.SG_modidx,,tbl.Opto)
% boxchart(tbl_avg.Subject,tbl_avg.SG_modidx_z,'GroupByColor',tbl_avg.Opto )
% 
% subplot(1,2,2)
% boxchart(tbl.Subject,tbl.FG_modidx,'GroupByColor',tbl.Opto)
% boxchart(tbl_avg.Subject,tbl_avg.FG_modidx_z,'GroupByColor',tbl_avg.Opto )
archt_gv_opto=[];
eyfp_gv_opto=[];
archt_silencing_SG_modidx_z=(tbl_avg.SG_modidx_z(tbl_avg.Opto==1 & tbl_avg.Cohort==1));
archt_gv_opto=[archt_gv_opto (tbl_avg.Opto(tbl_avg.Opto==1 & tbl_avg.Cohort==1))'];
archt_nosilencing_SG_modidx_z=(tbl_avg.SG_modidx_z(tbl_avg.Opto==0 & tbl_avg.Cohort==1));
archt_gv_opto=[archt_gv_opto (tbl_avg.Opto(tbl_avg.Opto==0 & tbl_avg.Cohort==1))'];

archt_SG_modidx_z=[archt_silencing_SG_modidx_z' archt_nosilencing_SG_modidx_z'];

archt_silencing_FG_modidx_z=(tbl_avg.FG_modidx_z(tbl_avg.Opto==1 & tbl_avg.Cohort==1));
archt_nosilencing_FG_modidx_z=(tbl_avg.FG_modidx_z(tbl_avg.Opto==0 & tbl_avg.Cohort==1));

archt_FG_modidx_z=[archt_silencing_FG_modidx_z' archt_nosilencing_FG_modidx_z'];

eyfp_silencing_SG_modidx_z=(tbl_avg.SG_modidx_z(tbl_avg.Opto==1 & tbl_avg.Cohort==0));
eyfp_gv_opto=[eyfp_gv_opto (tbl_avg.Opto(tbl_avg.Opto==1 & tbl_avg.Cohort==0))'];
eyfp_nosilencing_SG_modidx_z=(tbl_avg.SG_modidx_z(tbl_avg.Opto==0 & tbl_avg.Cohort==0));
eyfp_gv_opto=[eyfp_gv_opto (tbl_avg.Opto(tbl_avg.Opto==0 & tbl_avg.Cohort==0))'];

eyfp_SG_modidx_z=[eyfp_silencing_SG_modidx_z' eyfp_nosilencing_SG_modidx_z'];

eyfp_silencing_FG_modidx_z=(tbl_avg.FG_modidx_z(tbl_avg.Opto==1 & tbl_avg.Cohort==0));
eyfp_nosilencing_FG_modidx_z=(tbl_avg.FG_modidx_z(tbl_avg.Opto==0 & tbl_avg.Cohort==0));

eyfp_FG_modidx_z=[eyfp_silencing_FG_modidx_z' eyfp_nosilencing_FG_modidx_z'];

%% Plotting z_scores ***This is where you figured out the points connecting  pairwise comparison***
ff=figure(104);
clf;
fontsz=20;
ff.Position = [200, 200, 1000, 700];  % [x, y, width, height]


subplot(1,4,3)
boxplot(archt_SG_modidx_z, archt_gv_opto,'Labels', {'No silencing', 'Silencing'})
xtickangle(45);
hold on
%Add the lines
for ii=1:length(archt_silencing_SG_modidx_z)
    x=archt_silencing_SG_modidx_z(ii);
    y=archt_nosilencing_SG_modidx_z(ii);
    plot([1, 2], [y, x], '-o', 'Color', BC_color_genertor('oxford_blue',0.3), 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor',BC_color_genertor('oxford_blue') , 'MarkerEdgeColor', BC_color_genertor('oxford_blue'), 'Marker', 'o');  % Traces a line to the first point  
end
%Significance plot
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
plot(xt([1 2]), [1 1]*yt(end)*1.07, '-k',mean(xt([1 2])),yt(end)*1.08, '*k') %one asteriks
ylim([0 100])
title('SG');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
%Add colors
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
yticks([0:20:100]);
set(gca,'FontSize',fontsz);
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(1,4,4)
boxplot(archt_FG_modidx_z, archt_gv_opto,'Labels', {'No silencing', 'Silencing'})
xtickangle(45);
hold on
for ii=1:length(archt_silencing_FG_modidx_z)
    
    x=archt_silencing_FG_modidx_z(ii);
    y=archt_nosilencing_FG_modidx_z(ii);
    plot([1, 2], [y, x], '-o', 'Color', BC_color_genertor('oxford_blue',0.3), 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor',BC_color_genertor('oxford_blue') , 'MarkerEdgeColor', BC_color_genertor('oxford_blue'), 'Marker', 'o');  % Traces a line to the first point
    
end
title('FG');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
yticks([0:20:100]);
set(gca,'FontSize',fontsz);
ylim([0 100]);
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

%EYFP

subplot(1,4,1)
boxplot(eyfp_SG_modidx_z, eyfp_gv_opto,'Labels', {'No silencing', 'Silencing'})
xtickangle(45);
hold on
for ii=1:length(eyfp_silencing_SG_modidx_z)
    
    x=eyfp_silencing_SG_modidx_z(ii);
    y=eyfp_nosilencing_SG_modidx_z(ii);
    plot([1, 2], [y, x], '-o', 'Color', BC_color_genertor('oxford_blue',0.1), 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor',BC_color_genertor('oxford_blue') , 'MarkerEdgeColor', BC_color_genertor('oxford_blue'), 'Marker', 'o');  % Traces a line to the first point
    
end
title('SG');
colors = [BC_color_genertor('swamp_green');  % RGB values for Group 2
    BC_color_genertor('torment_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
ylabel('Z-Score (AU)');
yticks([0:20:100]);
set(gca,'FontSize',fontsz);
ylim([0 100]);
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(1,4,2)
boxplot(eyfp_FG_modidx_z, eyfp_gv_opto,'Labels', {'No silencing', 'Silencing'})
xtickangle(45);
hold on
for ii=1:length(eyfp_silencing_FG_modidx_z)
    
    x=eyfp_silencing_FG_modidx_z(ii);
    y=eyfp_nosilencing_FG_modidx_z(ii);
    plot([1, 2], [y, x], '-o', 'Color', BC_color_genertor('oxford_blue',0.3), 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor',BC_color_genertor('oxford_blue') , 'MarkerEdgeColor', BC_color_genertor('oxford_blue'), 'Marker', 'o');  % Traces a line to the first point
    
end
title('FG');
colors = [BC_color_genertor('swamp_green');  % RGB values for Group 2
    BC_color_genertor('torment_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
yticks([0:20:100]);
set(gca,'FontSize',fontsz);
ylim([0 100]);
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

% Adjust figure properties
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 % Resize the figure (optional)
% fig.Position = [200, 200, 1000, 700];  % [x, y, width, height]
 %saveas(gcf, 'FG_NormModIdx.png');%%%% fill in nice plotting %%%%%%%%% avg lme
%% Collecting data for the amplitud_velocity analysis


stats_avg.SG_mod = fitlme(tbl_avg, 'SG_modidx_z ~ 1+ Opto*Cohort + (1|Subject)')


stats_Arch_avg.SG_mod = fitlme(tbl_avg_archT, 'SG_modidx_z ~ 1+ Opto + (1|Subject)')

%% Putting all the the data for speed- amplitude analysis together

ArchtInhbSpdPwrTbl=BC_SpdPwrTblGenerator(out_archt,'inhb');
ArchtNoInhbSpdPwrTbl=BC_SpdPwrTblGenerator(out_archt,'noInhb');
ArchtRunningSpdPwrTbl=BC_SpdPwrTblGenerator(out_archt,'running');

CtrlInhbSpdPwrTbl=BC_SpdPwrTblGenerator(out_eyfp,'inhb');
CtrlNoInhbSpdPwrTbl=BC_SpdPwrTblGenerator(out_eyfp,'noInhb');
CtrlRunningSpdPwrTbl=BC_SpdPwrTblGenerator(out_eyfp,'running');
%% New version
ArchtInhbSpdPwrTbl=BC_SpdPwrDcmtdTblExtractor(out_archt,'Inhb');
ArchtNoInhbSpdPwrTbl=BC_SpdPwrDcmtdTblExtractor(out_archt,'NoInhb');
ArchtRunningSpdPwrTbl=BC_SpdPwrDcmtdTblExtractor(out_archt,'Running');

CtrlInhbSpdPwrTbl=BC_SpdPwrDcmtdTblExtractor(out_eyfp,'Inhb');
CtrlNoInhbSpdPwrTbl=BC_SpdPwrDcmtdTblExtractor(out_eyfp,'NoInhb');
CtrlRunningSpdPwrTbl=BC_SpdPwrDcmtdTblExtractor(out_eyfp,'Running');

%% stats
%Preparing data format
%Incorporating the opto column to archt
opto_archtInhb= logical(repmat(1,height(ArchtInhbSpdPwrTbl),1));
ArchtInhbSpdPwrTbl= addvars(ArchtInhbSpdPwrTbl,opto_archtInhb,'NewVariableNames','Opto' );
opto_archtNoInhb= logical(repmat(0,height(ArchtNoInhbSpdPwrTbl),1));
ArchtNoInhbSpdPwrTbl=addvars(ArchtNoInhbSpdPwrTbl,opto_archtNoInhb,'NewVariableNames','Opto' );
%Incorporating the opto column to control
opto_ctrlInhb= logical(repmat(1,height(CtrlInhbSpdPwrTbl),1));
CtrlInhbSpdPwrTbl= addvars(CtrlInhbSpdPwrTbl,opto_ctrlInhb,'NewVariableNames','Opto' );
opto_ctrlNoInhb= logical(repmat(0,height(CtrlNoInhbSpdPwrTbl),1));
CtrlNoInhbSpdPwrTbl=addvars(CtrlNoInhbSpdPwrTbl,opto_ctrlNoInhb,'NewVariableNames','Opto' );

%Combining archt tables into a single one and adding the cohort column
archTtbl= [ArchtInhbSpdPwrTbl;ArchtNoInhbSpdPwrTbl];
cohortArcht=repmat(1,height(archTtbl),1);
archTtbl=addvars(archTtbl,cohortArcht, 'NewVariableNames','Cohort');
%Combining ctrl tables into a single one and adding the cohort column
ctrlTtbl=[CtrlInhbSpdPwrTbl;CtrlNoInhbSpdPwrTbl];
cohortCtrl=repmat(0,height(ctrlTtbl),1);
ctrlTtbl=addvars(ctrlTtbl,cohortCtrl, 'NewVariableNames','Cohort');
ctrlTtbl.Subject=ctrlTtbl.Subject+max(unique(archTtbl.Subject)); %Adding the number of subjetcs from archt to ctrl table
%Making a big table
totalTtbl=[archTtbl;ctrlTtbl];
%Running the mixed effect model
stats.speed.tresh0=fitlme(totalTtbl,'Spd ~ 1+ Opto + (1|Cohort) + (1|Subject) ');
stats.archtSpeed.tresh0=fitlme(archTtbl,'Spd ~ 1+ Opto + (1|Cohort) + (1|Subject)');
stats.ctrlSpeed.tresh0=fitlme(ctrlTtbl,'Spd ~ 1+ Opto + (1|Cohort) + (1|Subject) ');

%trying to remove velocities below tresh
totalTblTresh=[];
archtTblTresh=[];
ctrlTblTresh=[];
for tresh=5:5:25;
    totalTblTresh.(['tresh' num2str(tresh)])=totalTtbl;
    totalTblTresh.(['tresh' num2str(tresh)])(totalTtbl.Spd<=tresh,:)=[];
    archtTblTresh.(['tresh' num2str(tresh)])=archTtbl;
    archtTblTresh.(['tresh' num2str(tresh)])(archTtbl.Spd<=tresh,:)=[];
    ctrlTblTresh.(['tresh' num2str(tresh)])=ctrlTtbl;
    ctrlTblTresh.(['tresh' num2str(tresh)])(ctrlTtbl.Spd<=tresh,:)=[];
end

%Running the mixed effect model with treshholds
for tresh=5:5:25;
    stats.speed.(['tresh' num2str(tresh)])=fitlme(totalTblTresh.(['tresh' num2str(tresh)]),'Spd ~ 1+ Opto + (1|Cohort) + (1|Subject)');
    stats.archtSpeed.(['tresh' num2str(tresh)])=fitlme(archtTblTresh.(['tresh' num2str(tresh)]),'Spd ~ 1+ Opto + (1|Cohort) + (1|Subject)');
    stats.ctrlSpeed.(['tresh' num2str(tresh)])=fitlme(ctrlTblTresh.(['tresh' num2str(tresh)]),'Spd ~ 1+ Opto + (1|Cohort) + (1|Subject)');
end
%% Graphing speed-theta amp
%Speed Vs Theta Amp
PcntTresh=75;

ArchtInhbSpd=ArchtInhbSpdPwrTbl.Spd;
ArchtInhbAmp=ArchtInhbSpdPwrTbl.ttaAmp;
ArchtInhbTresh= prctile(ArchtInhbSpd,PcntTresh)
outArchtIdxInhb=ArchtInhbSpd>ArchtInhbTresh;
ArchtInhbSpd=ArchtInhbSpd(outArchtIdxInhb);
ArchtInhbAmp=ArchtInhbAmp(outArchtIdxInhb);

ArchtNoInhbSpd=ArchtNoInhbSpdPwrTbl.Spd;
ArchtNoInhbAmp=ArchtNoInhbSpdPwrTbl.ttaAmp;
ArchtNoInhbTresh= prctile(ArchtNoInhbSpd,PcntTresh)
outArchtIdxNoInhb=ArchtNoInhbSpd>ArchtNoInhbTresh;
ArchtNoInhbSpd=ArchtNoInhbSpd(outArchtIdxNoInhb);
ArchtNoInhbAmp=ArchtNoInhbAmp(outArchtIdxNoInhb);

CtrlInhbSpd=CtrlInhbSpdPwrTbl.Spd;
CtrlInhbAmp=CtrlInhbSpdPwrTbl.ttaAmp;
CtrlInhbTresh= prctile(CtrlInhbSpd,PcntTresh)
outCtrlIdxInhb=CtrlInhbSpd>CtrlInhbTresh;
CtrlInhbSpd=CtrlInhbSpd(outCtrlIdxInhb);
CtrlInhbAmp=CtrlInhbAmp(outCtrlIdxInhb);

CtrlNoInhbSpd=CtrlNoInhbSpdPwrTbl.Spd;
CtrlNoInhbAmp=CtrlNoInhbSpdPwrTbl.ttaAmp;
CtrlNoInhbTresh= prctile(CtrlNoInhbSpd,PcntTresh)
outCtrlIdxNoInhb=CtrlNoInhbSpd>CtrlNoInhbTresh;
CtrlNoInhbSpd=CtrlNoInhbSpd(outCtrlIdxNoInhb);
CtrlNoInhbAmp=CtrlNoInhbAmp(outCtrlIdxNoInhb);

xmin=min([ArchtInhbTresh,ArchtNoInhbTresh,CtrlInhbTresh,CtrlNoInhbTresh]);

figure(101)
s1=subplot(2,2,2)
scatter(ArchtInhbSpd,ArchtInhbAmp, 'MarkerEdgeColor',BC_color_genertor('Archt_green'),'SizeData',11)
lsq1=lsline(s1);lsq1.Color='g',lsq1.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]);
xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
yticks([0 0.5 1]);
RAI=corr2(ArchtInhbSpd,ArchtInhbAmp);
title(sprintf('ArchT Silencing R=%f',RAI))
set(gca, 'TickDir', 'out');



s2=subplot(2,2,1)
scatter(ArchtNoInhbSpd,ArchtNoInhbAmp,'MarkerEdgeColor',BC_color_genertor('Oxford_blue'),'SizeData',11)
lsq2=lsline(s2);lsq2.Color='r',lsq2.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]); 
xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
yticks([0 0.5 1]);
RANI=corr2(ArchtNoInhbSpd,ArchtNoInhbAmp);
title(sprintf('ArchT No Silencing R=%f',RANI))
set(gca, 'TickDir', 'out');



s3=subplot(2,2,4)
scatter(CtrlInhbSpd,CtrlInhbAmp, 'MarkerEdgeColor',BC_color_genertor('Swamp_green'),'SizeData',11)
lsq3=lsline(s3);lsq3.Color='r',lsq3.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]); 
xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
yticks([0 0.5 1]);
RCI=corr2(CtrlInhbSpd,CtrlInhbAmp);
title(sprintf('Control Silencing R=%f',RCI))
set(gca, 'TickDir', 'out');


s4=subplot(2,2,3)
scatter(CtrlNoInhbSpd,CtrlNoInhbAmp,'MarkerEdgeColor',BC_color_genertor('Powder_blue'),'SizeData',11)
lsq4=lsline(s4);lsq4.Color='r',lsq4.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]); 
xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
yticks([0 0.5 1]);
RCNI=corr2(CtrlNoInhbSpd,CtrlNoInhbAmp);
title(sprintf('Control No Silencing R=%f',RCNI))
set(gca, 'TickDir', 'out');

box off;
fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1];
fig.Color = [1 1 1];         % Set background color to white
hold off;
%% Graphing theta amp vs sg amp

%Theta Amp Vs SG Amp

ArchtInhbAmpTta=ArchtInhbSpdPwrTbl.ttaAmp;
fillmissing(ArchtInhbAmpTta,'spline');
ArchtInhbAmpSG=ArchtInhbSpdPwrTbl.sgAmp;
fillmissing(ArchtInhbAmpSG,'spline');


ArchtNoInhbAmpTta=ArchtNoInhbSpdPwrTbl.ttaAmp;
fillmissing(ArchtNoInhbAmpTta,'spline');
ArchtNoInhbAmpSG=ArchtNoInhbSpdPwrTbl.sgAmp;
fillmissing(ArchtNoInhbAmpSG,'spline');


CrtlInhbAmpTta=CtrlInhbSpdPwrTbl.ttaAmp;
fillmissing(CrtlInhbAmpTta,'spline');
CrtlInhbAmpSG=CtrlInhbSpdPwrTbl.sgAmp;
fillmissing(CrtlInhbAmpSG,'spline');

CrtlNoInhbAmpTta=CtrlNoInhbSpdPwrTbl.ttaAmp;
fillmissing(CrtlNoInhbAmpTta,'spline');
CrtlNoInhbAmpSG=CtrlNoInhbSpdPwrTbl.sgAmp;
fillmissing(CrtlNoInhbAmpSG,'spline');


figure(102)
s1=subplot(2,2,2)
scatter(ArchtInhbAmpTta,ArchtInhbAmpSG, 'MarkerEdgeColor',BC_color_genertor('Archt_green'),'SizeData',11)
lsq1=lsline(s1);lsq1.Color='g',lsq1.LineWidth=2;
xlim([0 1]); ylim([0 1]);
xlabel('Theta Amplitude Normalized');ylabel('SG Amplitude Normalized')
yticks([0 0.5 1]);xticks([0 0.5 1]);
RAI=corr2(ArchtInhbAmpTta,ArchtInhbAmpSG)
title(sprintf('ArchT Silencing R=%f',RAI))
set(gca, 'TickDir', 'out');



s2=subplot(2,2,1)
scatter(ArchtNoInhbAmpTta,ArchtNoInhbAmpSG,'MarkerEdgeColor',BC_color_genertor('Oxford_blue'),'SizeData',11)
lsq2=lsline(s2);lsq2.Color='r',lsq2.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]); 
xlabel('Theta Amplitude Normalized');ylabel('SG Amplitude Normalized')
yticks([0 0.5 1]);xticks([0 0.5 1]);
RANI=corr2(ArchtNoInhbAmpTta,ArchtNoInhbAmpSG);
title(sprintf('ArchT No Silencing R=%f',RANI))
set(gca, 'TickDir', 'out');



s3=subplot(2,2,4)
scatter(CrtlInhbAmpTta,CrtlInhbAmpSG, 'MarkerEdgeColor',BC_color_genertor('Swamp_green'),'SizeData',11)
lsq3=lsline(s3);lsq3.Color='r',lsq3.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]); 
xlabel('Theta Amplitude Normalized');ylabel('SG Amplitude Normalized')
yticks([0 0.5 1]);xticks([0 0.5 1]);
RCI=corr2(CrtlInhbAmpTta,CrtlInhbAmpSG);
title(sprintf('Control Silencing R=%f',RCI))
set(gca, 'TickDir', 'out');


s4=subplot(2,2,3)
scatter(CrtlNoInhbAmpTta,CrtlNoInhbAmpSG,'MarkerEdgeColor',BC_color_genertor('Powder_blue'),'SizeData',11)
lsq4=lsline(s4);lsq4.Color='r',lsq4.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]); 
xlabel('Theta Amplitude Normalized');ylabel('SG Amplitude Normalized')
yticks([0 0.5 1]);xticks([0 0.5 1]);
RCNI=corr2(CrtlNoInhbAmpTta,CrtlNoInhbAmpSG);
title(sprintf('Control No Silencing R=%f',RCNI))
set(gca, 'TickDir', 'out');

box off;
fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1];
fig.Color = [1 1 1];         % Set background color to white
%% Graphing theta amp vs fg amp

%Theta Amp Vs FG Amp

ArchtInhbAmpTta=ArchtInhbSpdPwrTbl.ttaAmp;
fillmissing(ArchtInhbAmpTta,'spline');
ArchtInhbAmpFG=ArchtInhbSpdPwrTbl.fgAmp;
fillmissing(ArchtInhbAmpFG,'spline');


ArchtNoInhbAmpTta=ArchtNoInhbSpdPwrTbl.ttaAmp;
fillmissing(ArchtNoInhbAmpTta,'spline');
ArchtNoInhbAmpFG=ArchtNoInhbSpdPwrTbl.fgAmp;
fillmissing(ArchtNoInhbAmpFG,'spline');


CrtlInhbAmpTta=CtrlInhbSpdPwrTbl.ttaAmp;
fillmissing(CrtlInhbAmpTta,'spline');
CrtlInhbAmpFG=CtrlInhbSpdPwrTbl.fgAmp;
fillmissing(CrtlInhbAmpFG,'spline');

CrtlNoInhbAmpTta=CtrlNoInhbSpdPwrTbl.ttaAmp;
fillmissing(CrtlNoInhbAmpTta,'spline');
CrtlNoInhbAmpFG=CtrlNoInhbSpdPwrTbl.fgAmp;
fillmissing(CrtlNoInhbAmpFG,'spline');


figure(103)
s1=subplot(2,2,2)
scatter(ArchtInhbAmpTta,ArchtInhbAmpFG, 'MarkerEdgeColor',BC_color_genertor('Archt_green'),'SizeData',11)
lsq1=lsline(s1);lsq1.Color='g',lsq1.LineWidth=2;
xlim([0 1]); ylim([0 1]);
xlabel('Theta Amplitude Normalized');ylabel('FG Amplitude Normalized')
yticks([0 0.5 1]);xticks([0 0.5 1]);
RAI=corr2(ArchtInhbAmpTta,ArchtInhbAmpFG)
title(sprintf('ArchT Silencing R=%f',RAI))
set(gca, 'TickDir', 'out');



s2=subplot(2,2,1)
scatter(ArchtNoInhbAmpTta,ArchtNoInhbAmpFG,'MarkerEdgeColor',BC_color_genertor('Oxford_blue'),'SizeData',11)
lsq2=lsline(s2);lsq2.Color='r',lsq2.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]); 
xlabel('Theta Amplitude Normalized');ylabel('FG Amplitude Normalized')
yticks([0 0.5 1]);xticks([0 0.5 1]);
RANI=corr2(ArchtNoInhbAmpTta,ArchtNoInhbAmpFG);
title(sprintf('ArchT No Silencing R=%f',RANI))
set(gca, 'TickDir', 'out');



s3=subplot(2,2,4)
scatter(CrtlInhbAmpTta,CrtlInhbAmpFG, 'MarkerEdgeColor',BC_color_genertor('Swamp_green'),'SizeData',11)
lsq3=lsline(s3);lsq3.Color='r',lsq3.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]); 
xlabel('Theta Amplitude Normalized');ylabel('FG Amplitude Normalized')
yticks([0 0.5 1]);xticks([0 0.5 1]);
RCI=corr2(CrtlInhbAmpTta,CrtlInhbAmpFG);
title(sprintf('Control Silencing R=%f',RCI))
set(gca, 'TickDir', 'out');


s4=subplot(2,2,3)
scatter(CrtlNoInhbAmpTta,CrtlNoInhbAmpFG,'MarkerEdgeColor',BC_color_genertor('Powder_blue'),'SizeData',11)
lsq4=lsline(s4);lsq4.Color='r',lsq4.LineWidth=2;
%xlim([xmin 40]); ylim([0 1]); 
xlabel('Theta Amplitude Normalized');ylabel('FG Amplitude Normalized')
yticks([0 0.5 1]);xticks([0 0.5 1]);
RCNI=corr2(CrtlNoInhbAmpTta,CrtlNoInhbAmpFG);
title(sprintf('Control No Silencing R=%f',RCNI))
set(gca, 'TickDir', 'out');

box off;
fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1];
fig.Color = [1 1 1];         % Set background color to white

%% Distribution of the speed
tresh=5;
vlnCtrlNoInhb=CtrlNoInhbSpdPwrTbl.Spd;
vlnCtrlNoInhb=vlnCtrlNoInhb(vlnCtrlNoInhb>tresh);
vlnCtrlInhb=CtrlInhbSpdPwrTbl.Spd;
vlnCtrlInhb=vlnCtrlInhb(vlnCtrlInhb>tresh);
vlnArchtNoInhb=ArchtNoInhbSpdPwrTbl.Spd;
vlnArchtNoInhb=vlnArchtNoInhb(vlnArchtNoInhb>tresh);
vlnArchtInhb=ArchtInhbSpdPwrTbl.Spd;
vlnArchtInhb=vlnArchtInhb(vlnArchtInhb>tresh);
%Stats
[h1,p1]=ttest2(vlnCtrlNoInhb',vlnCtrlInhb');
[h2,p2]=ttest2(vlnArchtNoInhb,vlnArchtInhb);
[h3,p3]=ttest2(vlnArchtNoInhb,vlnArchtInhb);

dataviolin=([vlnCtrlNoInhb;vlnCtrlInhb;vlnArchtNoInhb;vlnArchtInhb]);
condition_names={'Control No Silencing','Control Silencing','ArchT No Silencing','ArchT Silencing'};
groups=([ones(1,length(vlnCtrlNoInhb)),2.*ones(1,length(vlnCtrlInhb)),3.*ones(1,length(vlnArchtNoInhb)),4.*ones(1,length(vlnArchtInhb))]);
c=[BC_color_genertor('Powder_blue');...
    BC_color_genertor('Swamp_green');...
    BC_color_genertor('Oxford_blue');...
    BC_color_genertor('Archt_green')];
% adding jittered scattered data same color boxplots for 2x2 data
figure(332)
h = daviolinplot(dataviolin,'groups',groups,'outsymbol','k+',...
    'boxcolors','same','colors',c,'scatter',1,'jitter',0,'violinalpha',0.7,'xtlabels', condition_names);
ylabel('Speed (cm/s)');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
ylim([0 40]);
set(gca,'FontSize',14);
yticks([0 10 20 30 40]);
% (*)Aesthetics (*)
box off;
set(gca, 'TickDir', 'out');
fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1];
fig.Color = [1 1 1];         % Set background color to white
hold off;

%stats


%% Graphing speed-theta amp
%Speed Vs Theta Amp
figure(101)
s1=subplot(2,2,2)
scatter(ArchtInhbSpdPwrTbl.Spd,ArchtInhbSpdPwrTbl.ttaPwr, 'MarkerEdgeColor',BC_color_genertor('Archt_green'),'SizeData',11)
lsq1=lsline(s1);lsq1.Color='r',lsq1.LineWidth=2;
%xlim([0 50]); ylim([0 1]); 
xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
title('ArchT Silencing')

s2=subplot(2,2,1)
scatter(ArchtNoInhbSpdPwrTbl.Spd,ArchtNoInhbSpdPwrTbl.ttaPwr,'MarkerEdgeColor',BC_color_genertor('Oxford_blue'),'SizeData',11)
lsq2=lsline(s2);lsq2.Color='r',lsq2.LineWidth=2;
%xlim([0 50]); ylim([0 1]); 
xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
title('ArchT No Silencing')

s3=subplot(2,2,4)
scatter(CtrlInhbSpdPwrTbl.Spd,CtrlInhbSpdPwrTbl.ttaPwr, 'MarkerEdgeColor',BC_color_genertor('Swamp_green'),'SizeData',11)
lsq3=lsline(s3);lsq3.Color='r',lsq3.LineWidth=2;
%xlim([0 50]); ylim([0 1]); 
xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
title('Control Silencing')

s4=subplot(2,2,3)
scatter(CtrlNoInhbSpdPwrTbl.Spd,CtrlNoInhbSpdPwrTbl.ttaPwr,'MarkerEdgeColor',BC_color_genertor('Powder_blue'),'SizeData',11)
lsq4=lsline(s4);lsq4.Color='r',lsq4.LineWidth=2;
%xlim([0 50]); ylim([0 1]); 
xlabel('Speed (cm/s)');ylabel('Theta Amplitude Normalized')
title('Control No Silencing')
%%
%Speed Vs SGamma Amp
figure(102)
subplot(2,2,2)
s=scatter(ArchtInhbSpdPwrTbl.Spd,ArchtInhbSpdPwrTbl.sgAmp, 'MarkerEdgeColor',BC_color_genertor('Archt_green'),'SizeData',11);
xlim([0 50]); ylim([0 1]); xlabel('Speed (cm/s)');ylabel('SG Amplitude (V)')
title('ArchT Silencing')
subplot(2,2,1)
s=scatter(ArchtNoInhbSpdPwrTbl.Spd,ArchtNoInhbSpdPwrTbl.sgAmp,'MarkerEdgeColor',BC_color_genertor('Oxford_blue'),'SizeData',11);
xlim([0 50]); ylim([0 1]); xlabel('Speed (cm/s)');ylabel('SG Amplitude (V)')
title('ArchT No Silencing')
subplot(2,2,4)
s=scatter(CtrlInhbSpdPwrTbl.Spd,CtrlInhbSpdPwrTbl.sgAmp, 'MarkerEdgeColor',BC_color_genertor('Swamp_green'),'SizeData',11);
xlim([0 50]); ylim([0 1]); xlabel('Speed (cm/s)');ylabel('SG Amplitude (V)')
title('Control Silencing')
subplot(2,2,3)
s=scatter(CtrlNoInhbSpdPwrTbl.Spd,CtrlNoInhbSpdPwrTbl.sgAmp,'MarkerEdgeColor',BC_color_genertor('Powder_blue'),'SizeData',11);
xlim([0 50]); ylim([0 1]); xlabel('Speed (cm/s)');ylabel('SG Amplitude (V)')
title('Control No Silencing')


%Speed Vs FGamma Amp
figure(103)
subplot(2,2,2)
s=scatter(ArchtInhbSpdPwrTbl.Spd,ArchtInhbSpdPwrTbl.FGAmp, 'MarkerEdgeColor',BC_color_genertor('Archt_green'),'SizeData',11);
xlim([0 50]); ylim([0 3*10^-4]); xlabel('Speed (cm/s)');ylabel('FG Amplitude (V)')
title('ArchT Silencing')
subplot(2,2,1)
s=scatter(ArchtNoInhbSpdPwrTbl.Spd,ArchtNoInhbSpdPwrTbl.FGAmp,'MarkerEdgeColor',BC_color_genertor('Oxford_blue'),'SizeData',11);
xlim([0 50]); ylim([0 3*10^-4]); xlabel('Speed (cm/s)');ylabel('FG Amplitude (V)')
title('ArchT No Silencing')
subplot(2,2,4)
s=scatter(CtrlInhbSpdPwrTbl.Spd,CtrlInhbSpdPwrTbl.FGAmp, 'MarkerEdgeColor',BC_color_genertor('Swamp_green'),'SizeData',11);
xlim([0 50]); ylim([0 3*10^-4]); xlabel('Speed (cm/s)');ylabel('FG Amplitude (V)')
title('Control Silencing')
subplot(2,2,3)
s=scatter(CtrlNoInhbSpdPwrTbl.Spd,CtrlNoInhbSpdPwrTbl.FGAmp,'MarkerEdgeColor',BC_color_genertor('Powder_blue'),'SizeData',11);
xlim([0 50]); ylim([0 3*10^-4]); xlabel('Speed (cm/s)');ylabel('FG Amplitude (V)')
title('Control No Silencing')

%% Putting all the the data for amplitude-amplitude analysis together

ArchtInhbSpdPwrTblTotal=BC_SpdPwrTblTotalGenerator(out_archt,'inhb');
ArchtNoInhbSpdPwrTblTotal=BC_SpdPwrTblTotalGenerator(out_archt,'noInhb');
ArchtRunningSpdPwrTblTotal=BC_SpdPwrTblTotalGenerator(out_archt,'running');

CtrlInhbSpdPwrTblTotal=BC_SpdPwrTblTotalGenerator(out_eyfp,'inhb');
CtrlNoInhbSpdPwrTblTotal=BC_SpdPwrTblTotalGenerator(out_eyfp,'noInhb');
CtrlRunningSpdPwrTblTotal=BC_SpdPwrTblTotalGenerator(out_eyfp,'running');

%% Graphing theta-Sgamma amplitude

figure(102)
clf;
subplot(2,2,4)
s=scatter(ArchtInhbSpdPwrTblTotal.TtaAmp,ArchtInhbSpdPwrTblTotal.SGAmp,[],ArchtInhbSpdPwrTblTotal.Subject,'SizeData',5);
s.MarkerEdgeAlpha = 0.1;
xlim([0 6*10^-4]); ylim([0 3.5*10^-4]); xlabel('Theta Amplitude (V)');ylabel('SG Amplitude (V)')
title('ArchT Silencing')

subplot(2,2,3)
s=scatter(ArchtNoInhbSpdPwrTblTotal.TtaAmp,ArchtNoInhbSpdPwrTblTotal.SGAmp,[],ArchtNoInhbSpdPwrTblTotal.Subject,'SizeData',5);
s.MarkerEdgeAlpha = 0.1;
xlim([0 6*10^-4]); ylim([0 3.5*10^-4]); xlabel('Theta Amplitude (V)');ylabel('SG Amplitude (V)')
title('ArchT No Silencing')

subplot(2,2,2)
s=scatter(CtrlInhbSpdPwrTblTotal.TtaAmp,CtrlInhbSpdPwrTblTotal.SGAmp,[],CtrlInhbSpdPwrTblTotal.Subject,'SizeData',5);
s.MarkerEdgeAlpha = 0.1;
xlim([0 6*10^-4]); ylim([0 3.5*10^-4]); xlabel('Theta Amplitude (V)');ylabel('SG Amplitude (V)')
title('Control Silencing')

subplot(2,2,1)
s=scatter(CtrlNoInhbSpdPwrTblTotal.TtaAmp,CtrlNoInhbSpdPwrTblTotal.SGAmp,[],CtrlNoInhbSpdPwrTblTotal.Subject,'SizeData',5);
s.MarkerEdgeAlpha = 0.1;
xlim([0 6*10^-4]); ylim([0 3.5*10^-4]); xlabel('Theta Amplitude (V)');ylabel('SG Amplitude (V)')
title('Control No Silencing')
%Amplitude ThetaVs FG
figure(103)
clf;
subplot(2,2,4)
s=scatter(ArchtInhbSpdPwrTblTotal.TtaAmp,ArchtInhbSpdPwrTblTotal.FGAmp,[],ArchtInhbSpdPwrTblTotal.Subject,'SizeData',5);
s.MarkerEdgeAlpha = 0.1;
xlim([0 6*10^-4]); ylim([0 3.5*10^-4]); xlabel('Theta Amplitude (V)');ylabel('FG Amplitude (V)')
title('ArchT Silencing')

subplot(2,2,3)
s=scatter(ArchtNoInhbSpdPwrTblTotal.TtaAmp,ArchtNoInhbSpdPwrTblTotal.FGAmp,[],ArchtNoInhbSpdPwrTblTotal.Subject,'SizeData',5);
s.MarkerEdgeAlpha = 0.1;
xlim([0 6*10^-4]); ylim([0 3.5*10^-4]); xlabel('Theta Amplitude (V)');ylabel('FG Amplitude (V)')
title('ArchT No Silencing')

subplot(2,2,2)
s=scatter(CtrlInhbSpdPwrTblTotal.TtaAmp,CtrlInhbSpdPwrTblTotal.FGAmp,[],CtrlInhbSpdPwrTblTotal.Subject,'SizeData',5);
s.MarkerEdgeAlpha = 0.1;
xlim([0 6*10^-4]); ylim([0 3.5*10^-4]); xlabel('Theta Amplitude (V)');ylabel('FG Amplitude (V)')
title('Control Silencing')

subplot(2,2,1)
s=scatter(CtrlNoInhbSpdPwrTblTotal.TtaAmp,CtrlNoInhbSpdPwrTblTotal.FGAmp,[],CtrlNoInhbSpdPwrTblTotal.Subject,'SizeData',5);
s.MarkerEdgeAlpha = 0.1;
xlim([0 6*10^-4]); ylim([0 3.5*10^-4]); xlabel('Theta Amplitude (V)');ylabel('FG Amplitude (V)')
title('Control No Silencing')

