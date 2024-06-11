%Channel exploratory sandbox

%% Parameters:
plot=00; % Activate to show the plots, not sure if is going to be necessary
save_output=00; %Activate if you wan to save the output figures in a folder
experment_type= 00; %00 for linear track, 01 for NOPR , and 02 for fear conditioning. 

%% Load all the CSC present in the folder
csc_list=dir('*ncs*');
csc_list=csc_list.name;
for iC= 1:length(csc_list);

%% Plot the raw traces from all the CSC. Filter in theta, SG and FG

%Linear track: Highlght the times of the running epochs (light and no light) of the mouse

%Open Field: Find the times when the mouse was in either of  the sleep states
%highlight them in the graphs
%% Linear track: Look for the theta, SG and FG powers comparison during light and non-light running epochs. Add PAC comparison

%% Open field: 

