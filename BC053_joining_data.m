evts1= LoadEvents([])
Chan_to_use={};
cfg = [];
cfg.fc = Chan_to_use; 
cfg.desired_sampling_frequency = 2000; % helps with speed. 

csc = MS_LoadCSC(cfg);


csc1=[];
for ii=1:7
    lfp_chan=['CSC' num2str(ii) '.ncs']
    cfg = [];
    cfg.fc = {lfp_chan};
    csc1.num2str(ii) = MS_LoadCSC(cfg);
end