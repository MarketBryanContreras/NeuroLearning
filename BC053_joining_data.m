evts1= LoadEvents([])
Chan_to_use={'CSC1.ncs' 'CSC2.ncs' 'CSC3.ncs' 'CSC4.ncs' 'CSC5.ncs' 'CSC6.ncs' 'CSC7.ncs'};
cfg = [];
cfg.fc = Chan_to_use; 
cfg.desired_sampling_frequency = 2000; % helps with speed. 

csc1 = MS_LoadCSC(cfg);

evts2= LoadEvents([])
cfg = [];
cfg.fc = Chan_to_use; 
cfg.desired_sampling_frequency = 2000; % helps with speed. 

csc2 = MS_LoadCSC(cfg);
