function graph_phase_amp=BC_plot_modidx(phase_amp, color, mod_idx,phi_name, period);
%% BC_ModIdx: It bins the index as a continuous variable between two signals. Indended for theta - gamma modulation as per Tort et al. 2010/2009/2008
%
%    %    Inputs: 
%    - phi_amp: [1 x nSamples] the mean amplitude x phase bins
%
%
%
%    Outputs: 
%    - ModIdx: [double] modulation index from Tort et al. 
%
%    
%
%
%
%
% BC 2023-08-24   initial version 
%
%% Assigning color pallete

Oxford_blue = [0.039 0.137 0.259];

%% Plot the phase amplitude coupling
figure (1)
clf
%subplot(2,5,[1:4,6:9]);
div=(length(phase_amp))/2;
phi_bins = -pi:pi/div:pi; 
deg_bins= (phi_bins(2:end)+pi)*(180/pi); %Here I am adding pi radinas to correct for the hilbert transformation and converting to bins in degrees
deg_bins=[deg_bins (deg_bins+360)];
bar(deg_bins,[phase_amp fliplr(phase_amp)],'FaceColor',color, 'EdgeColor',color )
%xlabel('Theta phase(Deg)');
ylabel(['Normalized mean ',phi_name,' amplitude'],'Fontsize', 12);
%title('Theta-phase-SG-amplitude', 'Fontsize', 20);
custom_tick=[90:90:720];
xticks([0 360 720]);
yticks([0.05 0.055 0.06])
ylim([0.05 0.061]);
set(gca,'FontSize', 18)
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot
%--to do-- Put title as a globl title and not only on this subplot
%title(['Thetha-phase modulation of ',phi_name,' amplitude during ',period ],'Fontsize',14)
%Bar plot mod idx
% subplot(2,5,[5,10])
% bar(mod_idx,'FaceColor',Oxford_blue, 'EdgeColor',Oxford_blue);
% max_y= ceil(log10(mod_idx));
% ylim([0 10^max_y]);
% ylabel('Modulation index','Fontsize', 12);
% set(gca, 'TickDir', 'out');  % Move ticks outside the plot
% box off;                     % Turn off the box around the plot
%Guiding phase
% subplot (3,5,[11:14]);
% Fs1= 1000;
% F1= 1;
% twin=[0 1];
% f_tvec= twin(1):1/Fs1:twin(2);
% fake_phase=[0:360/Fs1:360];
% fake_phase=[fake_phase (fake_phase+360)];
% fake_signal=sin(2*pi*F1*f_tvec);
% fake_signal=[fake_signal fake_signal];
% plot(fake_phase,fake_signal,'Color',Oxford_blue);
% xlim([0 fake_phase(end)]);
% ylabel('Theta amplitude','Fontsize', 12);
% xlabel('Theta phase (Deg)','Fontsize', 12);
% xticks(custom_tick);
% % Adjust plot aesthetics
% set(gca, 'TickDir', 'out');  % Move ticks outside the plot
% box off;                     % Turn off the box around the plot
% % Adjust figure properties
fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1];         % Set background color to white
fig.Position = [100, 100, 600, 550];  % [x, y, width, height]
%% Plotting only the curve
figure (2)
x = 1:length(phase_amp);
xx = linspace(1, length(phase_amp), 100);
smoothed_data = spline(x, phase_amp, xx);
%plot([x x+max(x)], [phase_amp flip(phase_amp)], 'b-', 'DisplayName', 'Original');
%hold on;
plot([xx xx+max(xx)], [smoothed_data flip(smoothed_data)], 'Color',color, 'DisplayName', 'Control FG No Light ', 'LineWidth',2);
ylim([0.05 0.061])
xlim([1 36])
xticks([1 18 36]);
xticklabels({'0' '360' '720'})
yticks([0.05 0.055 0.06])
set(gca,'FontSize', 18)
set(gca, 'TickDir', 'out');  % Move ticks outside the plot
box off;                     % Turn off the box around the plot

legend;
hold on;
fig = gcf;                   % Get current figure handle
fig.Color = [1 1 1];         % Set background color to white
fig.Position = [600, 100, 600, 550];  % [x, y, width, height]
