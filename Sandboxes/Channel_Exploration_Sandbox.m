%Channel exploratory sandbox

%% Dynamic loader
BC_linearTrack_dynamicLoader('experimental'); %"experimental" for archT, "control" 
%% Parameters
plot_flag = 1; % switch to 0 if you want to supress verification figures.
store_outputs=1; %switch to 0 for not storing the figures
time_maze_start = 30;
min_trial_dur = 0.5;
mouse_group=1; %1 for ArchT and 2 for eYFP. This just modify color of the plots
%% Generating color pallet for the mouse group
if mouse_group==1
    gcolors = [BC_color_genertor('powder_blue');  % RGB values for Group 2
        BC_color_genertor('archt_green')]; % RGB values for Group1;
else
    gcolors = [BC_color_genertor('torment_blue');  % RGB values for Group 2
        BC_color_genertor('swamp_green')]; % RGB values for Group1;
end
%% Loop over sessions
% Removing non- desired recordings
inhib_dir = dir('*D4_INHIBITION*');% get all sessions with 'D4_INHIBITION'
indices_to_remove = contains({inhib_dir.name}, 'disconnected');% Find the indices of directories containing 'disconnected'
inhib_dir(indices_to_remove) = [];% Remove the directories containing 'disconnected'