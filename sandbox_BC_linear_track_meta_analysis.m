%% Sandbox Linear_track_meta_analysis
%Go to you new data folder. In the future this should contain a file with the structure of the data to be analyzed
%This should be able to run from my lienar track folder and go into every
%folder. Fot he moment I will jut switch to each subdirectory and collect
%the data as necessary

working_dir = 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\inter';
raw_data_dir = 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\testing';
cd(raw_data_dir);
% raw_data_dir= 'C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw'
% cd(raw_data_dir)
% cd(data_dir)
% data_dir=["C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC1807_07_12_2023_D4_INHIBITION";
% "C:\Users\bcont\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC1602_07_12_2023_D4_INHIBITION"   
% ];
% folder_names=["BC1807_07_12_2023_D4_INHIBITION";
% "BC1602_07_12_2023_D4_INHIBITION" ];
%% Little hack to collect you data into a structure that contains the data you require. This can be seen as your draft for future structuring you data
%This structure is intented to conatin the data session from different mice
s= struct();
%Determine from how many macie you have info session
% s.mice1= struct;
% s.mice2= struct;
% % Determine how many session you have from each mice. For now lets keep
% it to every mouse assuming only one session
% s.mice1.session1=struct;
% s.mice2.session1=struct;
 
%% Assign the required field to each structure by mice
function LN_structure=LN_createDataStruct(raw_data_dir)
mainS=struct();
%Get a list of the folders in the raw folder
sessionList = dir((raw_data_dir));
sessionList = sessionList([sessionList.isdir]);  % Keep only folders
%% Loop through the sessions
for si=length(sessionList):-1:3
    sessionFolder = sessionList(si).name;
    sessionPath = fullfile(raw_data_dir, sessionFolder);
    % Extract mouse and day information from the session folder name
    [mouseName, day] = LT_name_nDay(sessionFolder);
    % Load session data
    session_data= LN_loadSessionData(sessionPath);
    cd(sessionPath);
    %Decide wether to save it as new mouse or add it to an existent one
    if ~isfield(mainS, mouseName)
        mainS.(mouseName).(day)=struct;
    end
    %Store the session
    mainS.(mouseName).(day)=session_data;
end

%%
%% Loop to filter you data in the theta, SG and FG band and add it to the structure

%% Loop through the session and calculate the modidx of your data. Add it to the structure, caluclate the mean of your shift data and the std

%% Caluclate the Z score of your modidx for laser on and laser off ruuning

end
