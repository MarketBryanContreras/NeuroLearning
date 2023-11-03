function [mouseName, day] = LT_name_nDay(sessionFolder)
%% LT_name_day: Extracts the name and day of the recording
% Inputs: [string] Name of the folder where you wan tto extract name and
% day of session named with the convention
% name_month_day_year_dayOfTraining_inhibition
% 
% 
% 
% 
% 

%% Parse session folder name to extract mouse name and day information
    parts = strsplit(sessionFolder, '_');
    mouseName = parts{1};
    day = parts{5};

