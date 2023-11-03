 function sessionData = LT_loadSessionData(sessionPath)
%% LoadSessionData Creates structure with the data 
% Inputs:
% 
% Outputs:

%%
cd(sessionPath); %I will try to modify this line in the future so it does not has to change in every iteration
%Load the events of the sseion
sessionData=[];
sessionData.events=LoadEvents([]); 
%Load the csc of interest
cfg =[];
cfg.fc = {'CSC6.ncs'};% Need to modify to be able to select different CSC
sessionData.csc = MS_LoadCSC(cfg);
sessionData.corrected_time=sessionData.csc.tvec-min(sessionData.csc.tvec);
end