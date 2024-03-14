%% SUPER_BC_stats

% load(['C:\Users\ecarm\Downloads' filesep 'out_Archt_06_nov_23.mat'])
% load(['C:\Users\ecarm\Downloads' filesep 'out_eyfp_06_nov_23.mat'])
%% Adjust the name of the files to load
ArchT_file_name='out_ArchT_07_nov_23.mat';
eyfp_file_name='out_eyfp_07_nov_23.mat';
%% Dynamic loader
sys=computer;
if contains(sys,'PCWIN')
    % Dynamic based on computer user windows.
    load([getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter' filesep ,ArchT_file_name])
    load([getenv('USERPROFILE') '\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter' filesep ,eyfp_file_name])
elseif contains(sys,'MAC')
    % Dynamic based on computer user for mac
    load([getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'inter' filesep ,ArchT_file_name]);
    load([getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'inter' filesep ,eyfp_file_name]);
    inter_dir = [getenv('HOME')  filesep 'Williams Lab Dropbox' filesep 'Williams Lab Team Folder' filesep 'Bryan_DropBox' filesep 'CHRNA2_LINEAR_TRACK' filesep 'inter'];
else disp("This system is not supported on the dynamic loader yet, pleae add it or contact BC for it")
end

%% Prototype to load a file from setting above
% parts = strsplit(ArchT_file_name, '.');
% archt_var= parts{1};
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


archt_list = fieldnames(out_Archt_07_nov_23);
eyfp_list = fieldnames(out_eyfp_07_nov_23); 


for iSub = 1:length(archt_list)
    
subject = [subject iSub iSub];

cohort = [cohort 1 1]; 
opto = [opto 1 0]; 

SG_modidx_z = [SG_modidx_z out_Archt_07_nov_23.(archt_list{iSub}).D4.z_SGInhb_modidx out_Archt_07_nov_23.(archt_list{iSub}).D4.z_SGNoInhb_modidx]; % has both inhib and no_inhib

FG_modidx_z = [FG_modidx_z out_Archt_07_nov_23.(archt_list{iSub}).D4.z_FGInhb_modidx out_Archt_07_nov_23.(archt_list{iSub}).D4.z_FGNoInhb_modidx]; % has both inhib and no_inhib

end
tbl_avg_archT= table(subject', cohort', opto',SG_modidx_z', FG_modidx_z' ,'VariableNames', {'Subject', 'Cohort', 'Opto', 'SG_modidx_z', 'FG_modidx_z'});



for iSub = 1:length(eyfp_list)
    
subject = [subject iSub+length(archt_list) iSub+length(archt_list)];

cohort = [cohort 0 0]; 
opto = [opto 1 0]; 

SG_modidx_z = [SG_modidx_z out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.z_SGInhb_modidx out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.z_SGNoInhb_modidx]; % has both inhib and no_inhib

FG_modidx_z = [FG_modidx_z out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.z_FGInhb_modidx out_eyfp_07_nov_23.(eyfp_list{iSub}).D4.z_FGNoInhb_modidx]; % has both inhib and no_inhib

end

tbl_avg= table(subject', logical(cohort'), logical(opto'),SG_modidx_z', FG_modidx_z' ,'VariableNames', {'Subject', 'Cohort', 'Opto', 'SG_modidx_z', 'FG_modidx_z'});

%% stats for avg Z mod


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

%% Plotting z_scores
ff=figure(104);
clf;
subplot(1,4,1)
boxplot(archt_SG_modidx_z, archt_gv_opto,'Labels', {'No silencing', 'Silencing'})
xtickangle(45);
hold on
for ii=1:length(archt_silencing_SG_modidx_z)
    
    x=archt_silencing_SG_modidx_z(ii);
    y=archt_nosilencing_SG_modidx_z(ii);
    plot([1, 2], [y, x], '-o', 'Color', BC_color_genertor('oxford_blue',0.3), 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor',BC_color_genertor('oxford_blue') , 'MarkerEdgeColor', BC_color_genertor('oxford_blue'), 'Marker', 'o');  % Traces a line to the first point
    
end
title('SG ArchT');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(1,4,2)
boxplot(archt_FG_modidx_z, archt_gv_opto,'Labels', {'No silencing', 'Silencing'})
xtickangle(45);
hold on
for ii=1:length(archt_silencing_FG_modidx_z)
    
    x=archt_silencing_FG_modidx_z(ii);
    y=archt_nosilencing_FG_modidx_z(ii);
    plot([1, 2], [y, x], '-o', 'Color', BC_color_genertor('oxford_blue',0.3), 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor',BC_color_genertor('oxford_blue') , 'MarkerEdgeColor', BC_color_genertor('oxford_blue'), 'Marker', 'o');  % Traces a line to the first point
    
end
title('FG ArchT');
colors = [BC_color_genertor('archt_green');  % RGB values for Group 2
    BC_color_genertor('powder_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.5);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

%EFP

subplot(1,4,3)
boxplot(eyfp_SG_modidx_z, eyfp_gv_opto,'Labels', {'No silencing', 'Silencing'})
xtickangle(45);
hold on
for ii=1:length(eyfp_silencing_SG_modidx_z)
    
    x=eyfp_silencing_SG_modidx_z(ii);
    y=eyfp_nosilencing_SG_modidx_z(ii);
    plot([1, 2], [y, x], '-o', 'Color', BC_color_genertor('oxford_blue',0.1), 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor',BC_color_genertor('oxford_blue') , 'MarkerEdgeColor', BC_color_genertor('oxford_blue'), 'Marker', 'o');  % Traces a line to the first point
    
end
title('SG control');
colors = [BC_color_genertor('swamp_green');  % RGB values for Group 2
    BC_color_genertor('torment_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

subplot(1,4,4)
boxplot(eyfp_FG_modidx_z, eyfp_gv_opto,'Labels', {'No silencing', 'Silencing'})
xtickangle(45);
hold on
for ii=1:length(eyfp_silencing_FG_modidx_z)
    
    x=eyfp_silencing_FG_modidx_z(ii);
    y=eyfp_nosilencing_FG_modidx_z(ii);
    plot([1, 2], [y, x], '-o', 'Color', BC_color_genertor('oxford_blue',0.3), 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor',BC_color_genertor('oxford_blue') , 'MarkerEdgeColor', BC_color_genertor('oxford_blue'), 'Marker', 'o');  % Traces a line to the first point
    
end
title('FG control');
colors = [BC_color_genertor('swamp_green');  % RGB values for Group 2
    BC_color_genertor('torment_blue')]; % RGB values for Group1;
h = findobj(gca, 'Tag', 'Box');  % Get handles to the box objects
for aa = 1:numel(h)
    patch(get(h(aa), 'XData'), get(h(aa), 'YData'), colors(aa, :), 'FaceAlpha', 0.8);
end
% Customize line and whisker colors
set(findobj(gca,'Type','Line'),'Color',[0.2 0.2 0.2]);
% Adjust plot aesthetics
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

% Adjust figure properties
 suptitle('\bf MOD IDX Z scores');           % Add a general title above the entire figure
 fig = gcf;                   % Get current figure handle
 fig.Color = [1 1 1];         % Set background color to white
 % Resize the figure (optional)
 fig.Position = [200, 200, 1000, 700];  % [x, y, width, height]
 %saveas(gcf, 'FG_NormModIdx.png');%%%% fill in nice plotting %%%%%%%
%% avg lme


stats_avg.SG_mod = fitlme(tbl_avg, 'SG_modidx_z ~ 1+ Opto*Cohort + (1|Subject)')


stats_Arch_avg.SG_mod = fitlme(tbl_avg_archT, 'SG_modidx_z ~ 1+ Opto + (1|Subject)')

